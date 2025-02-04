#' Function to perform dosage calling with polyRAD
#'
#' @param vcf ToDo
#' @param ploidy ToDo
#' @param model ToDo
#' @param popStruct ToDo
#' @import polyRAD
#' @importFrom vcfR read.vcfR is.biallelic write.vcf
#' @importFrom R.utils gunzip
#'
polyRAD_dosage_call <- function(vcf, ploidy, model, p1 = NULL, p2 = NULL, backcross.gen = 0, intercross.gen = 0, selfing.gen = 0,  session) {
  
  # Variables
  vcf_path <- vcf
  ploidy <- ploidy
  model <- model
  n.bx <- backcross.gen
  n.inter <- intercross.gen
  n.self <- selfing.gen
  
  # Having some issues formatting the output when multiallelic SNPs is input, so excluding for now
  temp_vcf <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
  temp_vcf <- temp_vcf[is.biallelic(temp_vcf),]
  
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
  polyRAD_obj <- polyRAD::VCF2RADdata(file = temp_unzipped_path,
                                      min.ind.with.reads = 1,
                                      min.ind.with.minor.allele = 0,
                                      taxaPloidy = ploidy,
                                      contamRate = 0.001,
                                      phaseSNPs = FALSE
                                      
  )
  
  # Empty vcf from memory
  rm(temp_unzipped_path)
  
  # Perform QC
  
  # Test overdispersion
  # Need to see if we can adapt the to_test values based on the input VCF or RAD object values
  # Currently using a range of 2:20, but an error will be provbided if the optimal value is the max or min of this range
  # Will iterate through a wider range if this is the case until the OD is estimated
  to_test_min <- 0
  to_test_max <- 20
  OD <- polyRAD::TestOverdispersion(polyRAD_obj, to_test = to_test_min:to_test_max)$optimal
  while (is.na(OD)) {
    message("Optimal value is at the edge of the range, retesting with a wider range")
    to_test_min <- to_test_min + 15
    to_test_max <- to_test_max + 20
    OD <- TestOverdispersion(polyRAD_obj, to_test = to_test_min:to_test_max)$optimal
  }

  # Test HindHe
  myhindhe <- HindHe(polyRAD_obj)

  # Perform dosage calling
  # "If you expect that your species has high linkage disequilibrium, the functions IterateHWE_LD and IteratePopStructLD behave like IterateHWE and IteratePopStruct, respectively, but also update priors based on genotypes at linked loci."
  if (model == "IterateHWE") {
    rad <- IterateHWE(polyRAD_obj, tol = 1e-5, overdispersion = OD)
  } else if (model == "IteratePopStruct") {
    rad <- IteratePopStruct(polyRAD_obj, tol = 1e-5, overdispersion = OD)
  } else if (model == "IterateHWE_LD") {
    rad <- IterateHWE_LD(polyRAD_obj, tol = 1e-5, overdispersion = OD)
  } else if (model == "IteratePopStruct_LD") {
    rad <- IteratePopStruct_LD(polyRAD_obj, tol = 1e-5, overdispersion = OD)
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


#' Function to convert polyRAD output to VCF
#'
#' @param geno ToDo
#' @param model ToDo
#' @param vcf_path ToDo
#' @param hindhe.obj ToDo
#' @param ploidy ToDo
#' @param output.file ToDo
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom BIGr flip_dosage
#' @importFrom reshape2 melt dcast
#' @importFrom vcfR read.vcfR extract.gt
#'
polyRAD2vcf <- function(geno, model, vcf_path, hindhe.obj, ploidy, output.file, session) {
  # Making the VCF
  # Appending the polyRAD info to the original VCF file and exporting
  
  #Format dosages
  gt <- data.frame(flip_dosage(geno, ploidy = ploidy), check.names = FALSE)
  row.names(gt) <- sub("_[A,T,C,G,B]$", "", row.names(gt)) #Removing the appended allele to the SNP IDs
  
  #Format HindHe
  colnames(hindhe.obj) <- sub("_[A,T,C,G,B]$", "", colnames(hindhe.obj)) #Removing the appended allele to the column
  hh <- data.frame(t(colMeans(hindhe.obj, na.rm = TRUE)), check.names = FALSE)
  hh <- hh[, row.names(gt)]
  
  #Getting information from VCF
  og_vcf <- read.vcfR(vcf_path, verbose = FALSE)
  
  # Updating Header
  # Extracting the header (meta)
  meta <- og_vcf@meta
  # Getting the last lines for INFO and FORMAT fields
  FORMAT_index <- min(grep("^##FORMAT=", meta))
  # Custom INFO lines to add
  custom_info_lines <- c(
    '##INFO=<ID=HH,Number=1,Type=Float,Description="Hind/He for the locus in the RADdata object">'
  )
  # Custom FORMAT lines to add
  custom_format_lines <- c(
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype, where 1 is the count of alternate alleles">',
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
  
  #Add the FORMAT label for each SNP
  vcf_df$FORMAT <- paste("GT","UD",og_vcf@gt[1, "FORMAT"],sep=":")
  
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
    
    # Handle missing values separately
    genotype_matrix <- matrix(genotype_strings[dosage_matrix + 1], nrow = nrow(dosage_matrix), ncol = ncol(dosage_matrix))
    genotype_matrix[is.na(dosage_matrix)] <- "./." # Handle missing values
    
    # Retain row and column names
    rownames(genotype_matrix) <- rownames(dosage_matrix)
    colnames(genotype_matrix) <- colnames(dosage_matrix)
    
    return(genotype_matrix)
  }
  
  # Convert the dosage matrix to genotypes
  gt_df <- convert_dosage2gt(gt, ploidy)
  
  #Combine info from the matrices to form the VCF information for each sample
  # Combine the matrices into a single matrix with elements separated by ":"
  make_vcf_format <- function(..., separator = ":") {
    matrices <- list(...)
    n <- length(matrices)
    
    # Convert matrices to long form
    long_forms <- lapply(matrices, function(mat) {
      suppressMessages(reshape2::melt(mat, varnames = c("Row", "Col"), value.name = "Value"))
    })
    
    # Concatenate the elements
    combined_long <- long_forms[[1]]
    combined_long$Combined <- combined_long$Value
    
    for (i in 2:n) {
      combined_long$Combined <- paste(combined_long$Combined, long_forms[[i]]$Value, sep = separator)
    }
    
    # Convert back to wide form
    combined_wide <- suppressMessages(reshape2::dcast(combined_long, Row ~ Col, value.var = "Combined"))
    
    # Restore row and column names
    rownames(combined_wide) <- combined_wide$Row
    combined_wide$Row <- NULL
    colnames(combined_wide) <- colnames(matrices[[1]])
    
    return(as.matrix(combined_wide))
  }
  
  # Replace NA values in the data frame with "."
  gt[is.na(gt)] <- "."
  
  # Combine the matrices
  # First get each item that exists in the VCF FORMAT
  #format_fields <- strsplit(og_vcf@gt[1, "FORMAT"], ":")[[1]]
  geno_df <- make_vcf_format(gt_df,
                             gt,
                             as.matrix(format_df)
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
