#' Estimate the polyRAD overdispersion parameter efficiently
#'
#' Two-stage search: a coarse, roughly log-spaced grid followed by an integer
#' refinement between the neighbours of the coarse optimum. Rough genotype
#' estimates are computed once up front so that \code{TestOverdispersion}
#' does not recompute them on every call. Large values mean the read counts
#' are close to binomial (little overdispersion), so the search is capped at
#' \code{max_od}; if the optimum is at the cap, the cap is returned.
#'
#' @param polyRAD_obj RADdata object
#' @param max_od Upper bound for the overdispersion parameter
#'
#' @return numeric, the selected overdispersion parameter
#' @import polyRAD
#' @importFrom stats ppoints
#' @noRd
estimate_overdispersion <- function(polyRAD_obj, max_od = 1000) {
  # Rough genotype estimation under HWE, done once (same steps TestOverdispersion
  # would otherwise run internally on every call)
  od_obj <- polyRAD_obj
  if (is.null(od_obj$posteriorProb)) {
    od_obj <- AddAlleleFreqHWE(od_obj)
    od_obj <- AddGenotypePriorProb_HWE(od_obj)
    od_obj <- AddGenotypeLikelihood(od_obj)
    od_obj <- AddPloidyChiSq(od_obj)
    od_obj <- AddGenotypePosteriorProb(od_obj)
  }

  # Score each tested value with the same metric TestOverdispersion uses, so we
  # can pick the minimum ourselves (TestOverdispersion returns NA at the edges)
  score_od <- function(to_test) {
    res <- suppressMessages(utils::capture.output(
      out <- TestOverdispersion(od_obj, to_test = to_test)
    ))
    out$optimal <- NULL
    scores <- vapply(out, function(x) {
      x1 <- sort(x)
      sqrt(mean((x1 - ppoints(length(x1)))^2))
    }, numeric(1))
    stats::setNames(scores, to_test)
  }

  # Stage 1: coarse grid spanning 2 to max_od
  coarse <- unique(round(c(2, 4, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 80,
                           100, 125, 150, 200, 250, 300, 400, 500, 650, 800, max_od)))
  coarse <- coarse[coarse <= max_od]
  coarse_scores <- score_od(coarse)
  i <- which.min(coarse_scores)

  if (i == length(coarse)) {
    warning("Overdispersion optimum is at the upper bound (", max_od,
            "); read counts look close to binomial. Using ", max_od, ".")
    return(max_od)
  }

  # Stage 2: integer refinement between neighbouring coarse values (at most ~30 points)
  lo <- coarse[max(1, i - 1)]
  hi <- coarse[i + 1]
  fine <- unique(round(seq(lo, hi, length.out = min(hi - lo + 1, 30))))
  fine_scores <- score_od(fine)
  OD <- fine[which.min(fine_scores)]
  message("Optimal overdispersion value: ", OD)
  OD
}

#' Perform dosage calling with polyRAD
#'
#' Reads a VCF with read depth information, keeps only biallelic markers,
#' imports them into a polyRAD \code{RADdata} object, estimates the
#' overdispersion parameter and Hind/He statistics, and calls allele dosages
#' using the selected polyRAD model.
#'
#' @param vcf Path to the input VCF file (plain, gzipped or bgzipped). It must
#'   contain allele read depth information (\code{AD} FORMAT field) for each
#'   sample. Multiallelic markers are removed before the analysis. If
#'   \code{REF}/\code{ALT} are missing or coded as \code{A}/\code{B}, random
#'   distinct nucleotides are assigned temporarily so polyRAD can read the file.
#' @param ploidy Integer. Ploidy of the samples, passed to \code{taxaPloidy}
#'   in \code{\link[polyRAD]{VCF2RADdata}}.
#' @param model Character. polyRAD model used for dosage calling. One of
#'   \code{"IterateHWE"} (Hardy-Weinberg equilibrium, e.g. diversity panels),
#'   \code{"IteratePopStruct"} (population structure),
#'   \code{"IterateHWE_LD"} or \code{"IteratePopStruct_LD"} (as above, also
#'   using linkage disequilibrium with nearby loci to update priors), or
#'   \code{"PipelineMapping2Parents"} (biparental mapping populations). Any
#'   value other than the four \code{Iterate*} models runs
#'   \code{\link[polyRAD]{PipelineMapping2Parents}}.
#' @param p1 Character or \code{NULL}. Sample name of the donor parent; used
#'   only with \code{"PipelineMapping2Parents"}. Must match a sample name in
#'   the VCF.
#' @param p2 Character or \code{NULL}. Sample name of the recurrent parent; used
#'   only with \code{"PipelineMapping2Parents"}. Must match a sample name in
#'   the VCF.
#' @param backcross.gen Integer. Number of generations of backcrossing to the
#'   recurrent parent (\code{p2}); used only with
#'   \code{"PipelineMapping2Parents"}. Default 0.
#' @param intercross.gen Integer. Number of generations of random intermating
#'   within the population; used only with \code{"PipelineMapping2Parents"}.
#'   Default 0 (e.g. 1 for an F1 cross).
#' @param selfing.gen Integer. Number of generations of self-fertilization;
#'   used only with \code{"PipelineMapping2Parents"}. Default 0.
#' @param contamRate Numeric. Expected rate of sample cross-contamination,
#'   passed to \code{\link[polyRAD]{VCF2RADdata}}. Default 0.001.
#' @param min_ind_read Integer. Minimum number of samples with reads for a
#'   marker to be kept (\code{min.ind.with.reads} in
#'   \code{\link[polyRAD]{VCF2RADdata}}). Default 1.
#' @param min_ind_maf Integer. Minimum number of samples carrying reads of the
#'   minor allele for a marker to be kept (\code{min.ind.with.minor.allele} in
#'   \code{\link[polyRAD]{VCF2RADdata}}). Default 0.
#' @param tol Numeric. Convergence tolerance for the iterative polyRAD models,
#'   also passed to \code{\link[polyRAD]{VCF2RADdata}}. Default 1e-05.
#' @param session Shiny \code{session} object. Kept for compatibility with the
#'   Shiny module; not used inside the function.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{Genos}{Matrix of the most probable dosages, with markers in rows
#'     and samples in columns. Each row is named \code{<marker>_<allele>} and
#'     holds the dosage (0 to \code{ploidy}) of that allele; the other allele
#'     of each biallelic marker is omitted. Samples with no reads at a marker
#'     are \code{NA}.}
#'   \item{RADHindHe}{Matrix of Hind/He values (samples by loci) from
#'     \code{\link[polyRAD]{HindHe}}, used to flag loci and samples whose
#'     heterozygosity does not match the expected ploidy.}
#' }
#'
#' @examples
#' \dontrun{
#' vcf_file <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")
#' res <- polyRAD_dosage_call(vcf = vcf_file, ploidy = 2, model = "IterateHWE")
#' }
#'
#' @import polyRAD
#' @importFrom vcfR read.vcfR is.biallelic write.vcf
#' @importFrom R.utils gunzip
#'
#' @export 
polyRAD_dosage_call <- function(vcf, ploidy, model, p1 = NULL, p2 = NULL,
                                backcross.gen = 0, intercross.gen = 0, selfing.gen = 0,
                                contamRate = 0.001, min_ind_read = 1, min_ind_maf = 0,
                                tol = 1e-05, session) {
  
  # Variables
  vcf_path <- vcf
  n.bx <- backcross.gen
  n.inter <- intercross.gen
  n.self <- selfing.gen
  
  # Having some issues formatting the output when multiallelic SNPs is input, so excluding for now
  temp_vcf <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
  temp_vcf <- temp_vcf[is.biallelic(temp_vcf),]
  
  # Function to randomly assign non-matching bases for polyRAD
  assign_random_bases <- function(vcf) {
    # Get the number of rows in the VCF
    num_rows <- nrow(vcf@fix)
    
    # Generate random bases for REF and ALT
    bases <- c("A", "C", "G", "T")
    ref_bases <- sample(bases, num_rows, replace = TRUE)
    alt_bases <- character(num_rows) # Initialize an empty character vector
    
    # Ensure REF and ALT are different for each row
    for (i in 1:num_rows) {
      alt_bases[i] <- sample(bases[bases != ref_bases[i]], 1)
    }
    
    # Assign the new bases to the vcf object
    vcf@fix[, "REF"] <- ref_bases
    vcf@fix[, "ALT"] <- alt_bases
    
    return(vcf)
  }
  
  #Editing REF and ALT fields randomly temporarily if blank or AB
  ref_base <- temp_vcf@fix[1, "REF"]
  alt_base <- temp_vcf@fix[1, "ALT"]
  if (is.na(ref_base) || is.na(alt_base) || alt_base == "B") {
    temp_vcf <- assign_random_bases(temp_vcf)
  }
  
  # Adding filtered VCF as a temp object
  temp_vcf_path <- tempfile(fileext = ".vcf.gz")
  vcfR::write.vcf(temp_vcf, file = temp_vcf_path)
  rm(temp_vcf)
  
  # vcfR gzipped the VCF, but polyRAD cannot accept gzipped, only unzipped or bgzipped.
  # Unzipping the VCF (see if there is a workaround because this is wasteful for memory)
  temp_unzipped_path <- tempfile(fileext = ".vcf")
  R.utils::gunzip(temp_vcf_path, destname = temp_unzipped_path, overwrite = TRUE)
  rm(temp_vcf_path)
  
  # Load the VCF file as RADobject
  # Retaining all markers with at least 1 individual with reads (minimal filtering at this stage)
  # Need to determine the best contamRate to use; currently using the default
  polyRAD_obj <- VCF2RADdata(file = temp_unzipped_path,
                             min.ind.with.reads = min_ind_read,
                             min.ind.with.minor.allele = min_ind_maf,
                             taxaPloidy = ploidy,
                             contamRate = contamRate,
                             tol = tol,
                             phaseSNPs = FALSE
                             
  )
  
  # Empty vcf from memory
  rm(temp_unzipped_path)
  
  # Perform QC
  
  # Test overdispersion
  OD <- estimate_overdispersion(polyRAD_obj)
  
  # Test HindHe
  myhindhe <- HindHe(polyRAD_obj)
  
  # Perform dosage calling
  # "If you expect that your species has high linkage disequilibrium, the functions IterateHWE_LD and IteratePopStructLD behave like IterateHWE and IteratePopStruct, respectively, but also update priors based on genotypes at linked loci."
  if (model == "IterateHWE") {
    rad <- IterateHWE(polyRAD_obj, tol = tol, overdispersion = OD)
  } else if (model == "IteratePopStruct") {
    rad <- IteratePopStruct(polyRAD_obj, tol = tol, overdispersion = OD)
  } else if (model == "IterateHWE_LD") {
    rad <- IterateHWE_LD(polyRAD_obj, tol = tol, overdispersion = OD)
  } else if (model == "IteratePopStruct_LD") {
    rad <- IteratePopStruct_LD(polyRAD_obj, tol = tol, overdispersion = OD)
  } else {
    if (!is.null(p1)) {
      # First check that p1 is a valid parent
      if (!p1 %in% rownames(data.frame(polyRAD_obj$taxaPloidy))) {
        stop("Parent 1 is not a valid parent")
      }
      
      polyRAD_obj <- SetDonorParent(polyRAD_obj, p1)
    }
    if (!is.null(p2)) {
      # First check that p2 is a valid parent
      if (!p2 %in% rownames(data.frame(polyRAD_obj$taxaPloidy))) {
        stop("Parent 2 is not a valid parent")
      }
      
      polyRAD_obj <- SetRecurrentParent(polyRAD_obj, p2)
    }
    
    rad <- PipelineMapping2Parents(polyRAD_obj,
                                   overdispersion = OD,
                                   n.gen.backcrossing = n.bx,
                                   n.gen.intermating = n.inter,
                                   n.gen.selfing = n.self)
  }
  
  
  # Get discrete genotypes
  pg <- GetProbableGenotypes(rad, naIfZeroReads = TRUE)$genotypes
  
  # Free up memory
  rm(rad)
  rm(polyRAD_obj)
  
  return(list(Genos = t(pg), RADHindHe = myhindhe))
  
}


#' Convert polyRAD Genotypes to a VCF File
#'
#' Takes genotype dosage estimates from polyRAD and appends them to an existing
#' VCF, adding Hind/He summary statistics and custom FORMAT/INFO fields.
#' The original variant metadata (CHROM, POS, REF, ALT, etc.) are preserved
#' from a template VCF, and new polyRAD-derived fields are written for each
#' sample and locus.
#'
#' @param geno A marker-by-sample dosage matrix (or object coercible to a
#'   matrix) containing polyRAD genotype estimates. Row names must match the
#'   variant IDs in \code{vcf_path} (before allele suffixes), and columns
#'   correspond to samples.
#' @param model Character string indicating the polyRAD model used to fit the
#'   genotypes (e.g. \code{"Kriging"}, \code{"SampleDepth"}, etc.). This value
#'   is recorded in the VCF header for provenance.
#' @param vcf_path Path to the template VCF file used to provide variant
#'   metadata (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, and original
#'   FORMAT fields).
#' @param hindhe.obj Object containing Hind/He values per locus (typically a
#'   matrix or data frame from polyRAD). Columns must correspond to alleles
#'   for each marker; the function computes the mean Hind/He per marker and
#'   adds it as an \code{INFO} field (\code{HH}) in the output VCF.
#' @param ploidy Integer giving the ploidy level used for the polyRAD analysis.
#'   This is used to convert dosage values into genotype strings (e.g.
#'   \code{"0/1"}, \code{"1/1"}, \code{"0/0/1"}, etc.).
#' @param output.file Base name for the output VCF file (without extension).
#'   The function will append \code{".vcf"} to this name and write the
#'   resulting VCF to disk in the current working directory (unless a path
#'   is included).
#' @param session Optional Shiny \code{session} object. Currently included for
#'   compatibility with Shiny modules; not used internally in this function.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Flips and formats polyRAD dosages using \code{\link[BIGr]{flip_dosage}}.
#'   \item Harmonizes marker IDs between \code{geno}, \code{hindhe.obj}, and
#'     the template VCF by stripping allele suffixes.
#'   \item Reads the template VCF with \code{\link[vcfR]{read.vcfR}} and
#'     updates the header to include:
#'     \itemize{
#'       \item \code{##INFO=<ID=HH,...>} for Hind/He per locus.
#'       \item \code{##FORMAT=<ID=GT,...>} for polyRAD-estimated genotypes
#'         (in allele-count space). If the template VCF already has a
#'         \code{GT} field, it is kept and renamed to \code{GTO} (original
#'         genotype) in both the header and the FORMAT column.
#'       \item \code{##FORMAT=<ID=UD,...>} for the dosage count of reference
#'         alleles.
#'       \item Version tags for \code{polyRAD} and \code{BIGapp}.
#'     }
#'   \item Converts dosage values to genotype strings (\code{GT}-like,
#'     \code{"0/0"}, \code{"0/1"}, \code{"1/1"}, etc.) given the specified
#'     \code{ploidy}.
#'   \item Builds a new FORMAT field for each variant combining:
#'     \code{GT} (polyRAD), \code{UD}, and the original FORMAT fields from the
#'     template (with any original \code{GT} and \code{UD} renamed to
#'     \code{GTO} and \code{UDO}).
#'   \item Writes the updated header and variant table to a VCF file.
#' }
#'
#' Missing dosages are converted to a missing genotype string (e.g.
#' \code{"./."} for diploids) and represented as \code{"."} in the dosage
#' field. Variants are sorted by chromosome and numeric position before
#' writing.
#'
#' @return
#' No value is returned. The function is called for its side effect of writing
#' a VCF file to \code{paste0(output.file, ".vcf")}.
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom BIGr flip_dosage
#' @importFrom vcfR read.vcfR extract.gt
#'
#' @export 
polyRAD2vcf <- function(geno, model, vcf_path, hindhe.obj, ploidy, output.file, session) {
  # Making the VCF
  # Appending the polyRAD info to the original VCF file and exporting
  
  #Format dosages
  gt <- data.frame(flip_dosage(geno, ploidy = ploidy), check.names = FALSE)
  #row.names(gt) <- sub("_[A,T,C,G,B]$", "", row.names(gt)) #Removing the appended allele to the SNP IDs
  row.names(gt) <- sub("_[^_]+$", "", row.names(gt))
  
  #Format HindHe
  #colnames(hindhe.obj) <- sub("_[A,T,C,G,B]$", "", colnames(hindhe.obj)) #Removing the appended allele to the column
  hh <- data.frame(t(colMeans(hindhe.obj, na.rm = TRUE)), check.names = FALSE)
  hh <- hh[, row.names(gt)]
  
  #Getting information from VCF
  og_vcf <- read.vcfR(vcf_path, verbose = FALSE)
  
  # Updating Header
  # Extracting the header (meta)
  meta <- og_vcf@meta
  # Custom INFO lines to add
  custom_info_lines <- c(
    '##INFO=<ID=HH,Number=1,Type=Float,Description="Hind/He for the locus in the RADdata object">'
  )
  # The polyRAD genotype always occupies GT (and dosage UD). If the input VCF
  # already has these fields, they are kept but renamed with an "O" suffix
  # (original): GT -> GTO, UD -> UDO.
  og_format <- og_vcf@gt[, "FORMAT"]
  renamed_lines <- c(
    GT = '##FORMAT=<ID=GTO,Number=1,Type=String,Description="Original genotype from the input VCF (replaced in GT by the polyRAD genotype)">',
    UD = '##FORMAT=<ID=UDO,Number=1,Type=Integer,Description="Original UD field from the input VCF (replaced in UD by the polyRAD dosage)">'
  )
  for (id in names(renamed_lines)) {
    og_header <- grep(paste0("^##FORMAT=<ID=", id, ","), meta)
    if (length(og_header) > 0) {
      meta[og_header] <- renamed_lines[[id]]
    } else if (any(grepl(paste0("(^|:)", id, "(:|$)"), og_format))) {
      meta <- append(meta, renamed_lines[[id]], after = max(c(0, grep("^##FORMAT=", meta))))
    }
  }
  # Getting the first FORMAT line, to insert the polyRAD fields before it
  FORMAT_index <- min(grep("^##FORMAT=", meta))

  # Custom FORMAT lines to add
  custom_format_lines <- c(
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype estimated by polyRAD, where 1 is the count of alternate alleles">',
    '##FORMAT=<ID=UD,Number=1,Type=Integer,Description="Dosage count of reference alleles from polyRAD, where 0 = homozygous alternate">'
  )
  custom_version_lines <- c(
    paste0('##polyRAD','_',model,'=', packageVersion("polyRAD")),
    paste0('##BIGapp_polyRAD2vcf=', packageVersion("BIGapp"))
  )
  
  #polyrad_meta <- paste0('##PolyRADCommandLine.multidog=<ID=Multidog,Version="',
  #                     updog_version, '",CommandLine="> multidog(refmat = matrices$ref_matrix, sizemat = matrices$size_matrix, ploidy = ',ploidy,
  #                     ', model = ',model_selected,')">')
  #bigapp_meta <- paste0('##BIGappCommandLine.polyRAD2vcf=<ID=polyRAD2vcf,Version="',
  #                    packageVersion("BIGapp"), '",Data="',
  #                    Sys.time(),'", CommandLine="> updog2vcf(',deparse(substitute(multidog.object)),',',
  #                    output.file, ',',
  #                    updog_version,')">')
  
  # Inserting the new header lines into the VCF header
  meta <- append(meta, custom_format_lines, after = FORMAT_index-1) #FORMAT
  INFO_index <- max(grep("^##INFO=", meta))
  meta <- append(meta, custom_info_lines, after = INFO_index) #INFO
  FORMAT_index <- max(grep("^##FORMAT=", meta))
  meta <- append(meta, custom_version_lines, after = FORMAT_index)
  
  
  #Make a header separate from the dataframe
  vcf_header <- meta
  
  #Make the VCF df
  vcf_df <- data.frame(
    CHROM = data.frame(og_vcf@fix)$CHROM,
    POS = data.frame(og_vcf@fix)$POS,
    ID = data.frame(og_vcf@fix)$ID,
    REF = data.frame(og_vcf@fix)$REF,
    ALT = data.frame(og_vcf@fix)$ALT,
    QUAL = data.frame(og_vcf@fix)$QUAL,
    FILTER = data.frame(og_vcf@fix)$FILTER,
    INFO = data.frame(og_vcf@fix)$INFO,
    FORMAT = data.frame(og_vcf@gt)$FORMAT
  )
  
  # Replace NA values in the data frame with "."
  vcf_df[is.na(vcf_df)] <- "."
  
  #Get FORMAT info
  format_df <- og_vcf@gt[,-1]
  
  # Add the IDs as row names for subsetting
  row.names(vcf_df) <- vcf_df$ID
  row.names(format_df) <- vcf_df$ID
  vcf_df <- vcf_df[row.names(gt),]
  format_df <- format_df[row.names(gt),]
  
  #Add the INFO column for each SNP
  vcf_df$INFO <- paste0(vcf_df$INFO,";",
                        "HH=",t(hh))
  
  #Add the FORMAT label for each SNP (polyRAD GT first; original GT/UD renamed to GTO/UDO)
  og_format_renamed <- vcf_df$FORMAT
  for (id in c("GT", "UD")) {
    og_format_renamed <- gsub(paste0("(^|:)", id, "(?=:|$)"), paste0("\\1", id, "O"),
                              og_format_renamed, perl = TRUE)
  }
  vcf_df$FORMAT <- paste("GT", "UD", og_format_renamed, sep = ":")
  
  #Convert genotypes from dosage to gt
  precompute_genotype_strings <- function(ploidy) {
    genotype_strings <- character(ploidy + 1)
    # Generate the genotype string based on the dosage and ploidy
    # Updog uses the ref counts, which is not typical, so this corrects it
    for (dosage in 0:ploidy) {
      ref_count <- dosage
      alt_count <- ploidy - dosage
      genotype_strings[dosage + 1] <- paste(c(rep("0", ref_count), rep("1", alt_count)), collapse = "/")
    }
    return(genotype_strings)
  }
  
  # Apply the precomputed genotype strings to the matrix
  convert_dosage2gt <- function(dosage_matrix, ploidy) {
    dosage_matrix <- as.matrix(dosage_matrix)
    genotype_strings <- precompute_genotype_strings(ploidy)
    
    # Create missing GT
    missing_gt_string <- paste(rep(".", ploidy), collapse = "/")
    
    # Handle missing values separately
    genotype_matrix <- matrix(genotype_strings[dosage_matrix + 1], nrow = nrow(dosage_matrix), ncol = ncol(dosage_matrix))
    genotype_matrix[is.na(dosage_matrix)] <- missing_gt_string # Handle missing values
    
    # Retain row and column names
    rownames(genotype_matrix) <- rownames(dosage_matrix)
    colnames(genotype_matrix) <- colnames(dosage_matrix)
    
    return(genotype_matrix)
  }
  
  # Convert the dosage matrix to genotypes
  gt_df <- convert_dosage2gt(gt, ploidy)
  # Replace NA values in the data frame with "."
  gt[is.na(gt)] <- "."
  
  #Combine info from the matrices to form the VCF information for each sample
  # Combine the matrices into a single matrix with elements separated by ":"
  make_vcf_format <- function(matrix1, matrix2, matrix3) {
    # Check if input are matrices (optional, but good practice)
    if (!is.matrix(matrix1) || !is.matrix(matrix2) || !is.matrix(matrix3)) {
      stop("All inputs must be matrices.")
    }
    
    # Check if matrices have the same dimensions (optional, but recommended)
    if (!all(dim(matrix1) == dim(matrix2), dim(matrix1) == dim(matrix3))) {
      stop("Input matrices must have the same dimensions.")
    }
    
    # Check if matrices have row and column names (optional, based on problem description)
    if (is.null(rownames(matrix1)) || is.null(colnames(matrix1)) ||
        is.null(rownames(matrix2)) || is.null(colnames(matrix2)) ||
        is.null(rownames(matrix3)) || is.null(colnames(matrix3))) {
      stop("Input matrices must have row and column names.")
    }
    
    # 1. Convert matrices to vectors
    vector1 <- as.vector(matrix1)
    vector2 <- as.vector(matrix2)
    vector3 <- as.vector(matrix3)
    
    # 2. Combine the vectors element-wise with ":" as separator
    combined_vector <- paste(vector1, vector2, vector3, sep = ":")
    
    # 3. Reshape the combined vector back into a matrix with the original dimensions
    combined_matrix <- matrix(combined_vector, nrow = nrow(matrix1), ncol = ncol(matrix1))
    
    # 4. Convert the combined matrix to a dataframe
    combined_dataframe <- as.data.frame(combined_matrix)
    
    # 5. Set the row and column names for the dataframe
    rownames(combined_dataframe) <- rownames(matrix1)
    colnames(combined_dataframe) <- colnames(matrix1)
    
    return(combined_dataframe)
  }
  
  # Combine the matrices
  # First get each item that exists in the VCF FORMAT
  #format_fields <- strsplit(og_vcf@gt[1, "FORMAT"], ":")[[1]]
  geno_df <- make_vcf_format(gt_df,
                             as.matrix(gt),
                             format_df
  )
  
  #Combine the dataframes together
  vcf_df <- cbind(vcf_df,geno_df)
  
  # Add # to the CHROM column name
  colnames(vcf_df)[1] <- "#CHROM"
  
  # Sort
  vcf_df <- vcf_df[order(vcf_df[,1],as.numeric(as.character(vcf_df[,2]))),]
  
  # Write the header to the file
  # Make sure that .vcf is at the end of the file name
  output.file <- paste0(output.file, ".vcf")
  writeLines(vcf_header, con = output.file)
  
  # Append the dataframe to the file in tab-separated format
  suppressWarnings(
    write.table(vcf_df, file = output.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE)
  )
  
  
}
