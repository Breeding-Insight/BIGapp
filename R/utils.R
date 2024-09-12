# Modify the get_counts function to accept the MADC file path as an argument
get_counts <- function(madc_file, output_name) {
  # This function takes the MADC file as input and generates a Ref and Alt counts dataframe as output
  # Note: This assumes that the first 7 rows are not useful here like in the Strawberry DSt23-8501_MADC file

  # Read the madc file
  madc_df <- read.csv(madc_file, sep = ',', skip = 7, check.names = FALSE)

  # Retain only the Ref and Alt haplotypes
  filtered_df <- madc_df[!grepl("\\|AltMatch|\\|RefMatch", madc_df$AlleleID), ]

  #Remove extra text after Ref and Alt (_001 or _002)
  filtered_df$AlleleID <- sub("\\|Ref.*", "|Ref", filtered_df$AlleleID)
  filtered_df$AlleleID <- sub("\\|Alt.*", "|Alt", filtered_df$AlleleID)

  return(filtered_df)
}


#Get the alt, ref, and size matrix for use in Updog
#Add functionality here to stop the script if indentical() is False
get_matrices <- function(result_df) {
  #This function takes the dataframe of ref and alt counts for each sample, and converts them to ref, alt, and size(total count) matrices for Updog

  update_df <- result_df

  # Filter rows where 'AlleleID' ends with 'Ref'
  ref_df <- subset(update_df, grepl("Ref$", AlleleID))

  # Filter rows where 'AlleleID' ends with 'Alt'
  alt_df <- subset(update_df, grepl("Alt$", AlleleID))

  #remove alt or ref rows that do not have a counterpart in the other dataframe
  if (nrow(ref_df) > nrow(alt_df)) {
    ref_df <- ref_df[ref_df$CloneID %in% alt_df$CloneID,]
  } else if (nrow(ref_df) < nrow(alt_df)) {
    alt_df <- alt_df[alt_df$CloneID %in% ref_df$CloneID,]
  } else {
    alt_df <- alt_df[alt_df$CloneID %in% ref_df$CloneID,]
  }

  #Ensure that each has the same SNPs and that they are in the same order
  identical(alt_df$CloneID,ref_df$CloneID)

  ###Convert the ref and alt counts into matrices with the CloneID as the index
  #Set SNP names as index
  row.names(ref_df) <- ref_df$CloneID
  row.names(alt_df) <- alt_df$CloneID

  #Remove unwanted columns and convert to matrix
  #Probably best to just remove the column names that aren't wanted instead of the first 16 columns.
  ref_matrix <- as.matrix(ref_df[, -c(1:16)])
  alt_matrix <- as.matrix(alt_df[, -c(1:16)])

  #Make the size matrix by combining the two matrices
  size_matrix <- (ref_matrix + alt_matrix)

  #Count the number of cells with 0 count to estimate missing data
  # Count the number of cells with the value 0
  count_zeros <- sum(size_matrix == 0)

  # Print the result
  ratio_missing_data <- count_zeros / length(size_matrix)
  cat("Ratio of missing data =", ratio_missing_data, "\n")

  # Return the ref and alt matrices as a list
  matrices_list <- list(ref_matrix = ref_matrix, size_matrix = size_matrix)
  return(matrices_list)
}

#' Extract VCF info IDs
#' @param info_string to be documented
extract_info_ids <- function(info_string) {
  # Split the INFO string by ';'
  info_parts <- strsplit(info_string, ":")[[1]]
  # Extract the part before the '=' in each segment
  info_ids <- gsub("=.*", "", info_parts)
  return(info_ids)
}

convert_to_dosage <- function(gt) {
  # Split the genotype string
  alleles <- strsplit(gt, "[|/]")
  # Sum the alleles, treating NA values appropriately
  sapply(alleles, function(x) {
    if (any(is.na(x))) {
      return(NA)
    } else {
      return(sum(as.numeric(x), na.rm = TRUE))
    }
  })
}


#' Internal function
#'
#' @param genotypeMatrix fill description
#' @param maxK fill description
#' @param ploidy fill description
#'
#' @importFrom adegenet pop
#'
findK <- function(genotypeMatrix, maxK, ploidy) {
  # Convert the genotype matrix to a genlight object
  genlight_new <- new("genlight", t(genotypeMatrix),
                      ind.names = row.names(t(genotypeMatrix)),
                      loc.names = colnames(t(genotypeMatrix)),
                      ploidy = ploidy,
                      NA.char = NA)

  #Assign the populations as the sample names since there is no assumed populations
  pop(genlight_new) <- genlight_new@ind.names

  #Estimate number of clusters
  #Retain all pca for the find.clusters step. Retain as few as possible while maximizing variance captured for DAPC step.
  #Choose is the option to allow adegenet to select the best cluster number based on the BIC minimum
  #The default criterion is "diffNgroup, which is not necessarily the minimum BIC, but based on the sharp decrease of the BIC value.
  #Either way, this is a suggestion, and the number of clusters to use should be made with biology considerations.
  graphics.off() #Prevent plot from automatically displaying
  grp <- find.clusters(genlight_new, max.n.clust = maxK,
                       n.pca = nInd(genlight_new),
                       stat = "BIC",
                       criterion = "diffNgroup",
                       parallel = FALSE,
                       choose = FALSE)

  # Identify the best K based on lowest BIC
  bestK <- length(grp$size)

  # Create a BIC dataframe
  bicDF <- data.frame(K = 1:maxK, BIC = as.data.frame(grp$Kstat)$`grp$Kstat`)

  return(list(bestK = as.numeric(bestK), grp = grp, BIC = bicDF))

}

performDAPC <- function(genotypeMatrix, selected_K, ploidy) {

  #Convert matrix to genlight
  genlight_new <- new("genlight", t(genotypeMatrix),
                      ind.names = row.names(t(genotypeMatrix)),
                      loc.names = colnames(t(genotypeMatrix)),
                      ploidy = ploidy,
                      NA.char = NA)

  #Get groups based on specified cluster number (K)
  graphics.off() #Prevent plot from automatically displaying
  grp <- find.clusters(genlight_new, n.clust = selected_K,
                       n.pca = nInd(genlight_new),
                       stat = "BIC",
                       criterion = "diffNgroup",
                       parallel = FALSE,
                       choose = FALSE)

  # Find the optimal number of principal components
  #NOTE: The default n.da is K-1, but I have read previously to use #Samples - 1?
  dapc1 <- dapc(genlight_new, grp$grp,
                n.pca = nInd(genlight_new),
                n.da = nInd(genlight_new)-1,
                parallel = FALSE)

  a.score <- optim.a.score(dapc1, plot = FALSE)
  n.pca <- a.score$best

  # Perform DAPC with the best K
  finalDapc <- dapc(genlight_new, grp$grp, n.pca = n.pca, n.da = selected_K-1, parallel= FALSE)

  # Extract the membership probabilities
  Q <- as.data.frame(finalDapc$posterior)

  # Add cluster assignments to Q dataframe
  Q$Cluster_Assignment <- finalDapc$assign

  #a data.frame giving the contributions of original variables (alleles in the case of genetic data) to the principal components of DAPC.
  #dapc$var.contr

  # Return list containing BIC dataframe, Q dataframe w/ dapc assignments
  return(list(Q = Q, dapc = finalDapc))
}

#Heterozygosity function
calculate_heterozygosity <- function(genotype_matrix, ploidy = 2) {
  # Determine the heterozygous values based on ploidy
  heterozygous_values <- seq(1, ploidy - 1)

  # Create a logical matrix where TRUE represents heterozygous loci
  is_heterozygous <- sapply(genotype_matrix, function(x) x %in% heterozygous_values)

  # Count the number of heterozygous loci per sample, ignoring NAs
  heterozygosity_counts <- colSums(is_heterozygous, na.rm = TRUE)

  # Calculate the total number of non-NA loci per sample
  total_non_na_loci <- colSums(!is.na(genotype_matrix))

  # Compute the proportion of heterozygous loci
  heterozygosity_proportion <- heterozygosity_counts / total_non_na_loci

  # Create a dataframe with Sample ID and Observed Heterozygosity
  result_df <- data.frame(
    SampleID = colnames(genotype_matrix),
    ObservedHeterozygosity = heterozygosity_proportion,
    row.names = NULL,
    check.names = FALSE
  )

  return(result_df)
}

#' Updated MAF function
#'
#' @param df fill description
#' @param ploidy fill description
#'
#' @importFrom tibble rownames_to_column
#'
calculateMAF <- function(df, ploidy) {
  if (is.matrix(df)) {
    df <- as.data.frame(df)
  }

  #Convert the elements to numeric if they are characters
  df[] <- lapply(df, function(x) if(is.character(x)) as.numeric(as.character(x)) else x)

  allele_frequencies <- apply(df, 1, function(row) {
    non_na_count <- sum(!is.na(row))
    allele_sum <- sum(row, na.rm = TRUE)
    if (non_na_count > 0) {
      allele_sum / (ploidy * non_na_count)
    } else {
      NA
    }
  })

  maf <- ifelse(allele_frequencies <= 0.5, allele_frequencies, 1 - allele_frequencies)

  df$AF <- allele_frequencies
  df$MAF <- maf

  maf_df <- df[,c("AF", "MAF"), drop = FALSE]

  #Make the row names (SNP ID) the first column
  maf_df <- maf_df %>%
    rownames_to_column(var = "SNP_ID")

  return(maf_df)
}

# Function to calculate percentages for each genotype in each sample
calculate_percentages <- function(matrix_data, ploidy) {
  apply(matrix_data, 2, function(col) {
    counts <- table(col)
    prop <- prop.table(counts) * 100
    prop[as.character(0:ploidy)]  # Adjust the range based on the max value (consider entering the ploidy value explicitly for max_val)
  })
}

convert_genotype_counts <- function(df, ploidy, is_reference = TRUE) {
  if (is_reference) {
    # Convert from reference to alternate alleles
    return(abs(df - ploidy))
  } else {
    # Data already represents alternate alleles
    return(df)
  }
}

#' Internal function
#'
#' @param mat fill description
#'
#' @importFrom matrixcalc is.positive.definite
#'
posdefmat <- function(mat) {
  if (is.positive.definite(round(mat, 18))) {
    g = mat
  }
  else {
    g <-nearPD(mat)$mat
    warning("The matrix was adjusted for the nearest positive definite matrix")
  }
  return(g)
}

# Function to split INFO column and expand it into multiple columns
split_info_column <- function(info) {
  # Split the INFO column by semicolon
  info_split <- str_split(info, ";")[[1]]

  # Create a named list by splitting each element by equals sign
  info_list <- set_names(map(info_split, ~ str_split(.x, "=")[[1]][2]),
                         map(info_split, ~ str_split(.x, "=")[[1]][1]))

  return(info_list)
}

#' Read geno file
#'
#' @param file_path character indicanting path to file
#' @param requires which information is required from the VCF. Define the FORMAT or INFO letters. Example: c("GT", "DP", "PL")
#'
#' @importFrom vcfR read.vcfR
#' @importFrom shinyalert shinyalert
#'
read_geno_file <- function(file_path, requires = c("GT")){
  if (grepl("\\.csv$", file_path)) {
    geno <- read.csv(geno_path, header = TRUE, row.names = 1, check.names = FALSE)
    n_snps <- nrow(geno)
    return(list(geno, n_snps))

  } else if (grepl("\\.vcf$", file_path) || grepl("\\.gz$", file_path)) {

    #Convert VCF file if submitted
    vcf <- read.vcfR(file_path)

    all_requires <- vector()
    for(i in 1:length(requires))  all_requires[i] <- grepl(requires[i], vcf@fix[1,8]) | grepl(requires[i], vcf@gt[1,1])

    if(!all(all_requires)) {
      shinyalert(
        title = "Oops",
        text = paste("The VCF file does not contain required information:", requires[which(!all_requires)]),
        size = "xs",
        closeOnEsc = TRUE,
        closeOnClickOutside = FALSE,
        html = TRUE,
        type = "warning",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
        confirmButtonCol = "#004192",
        showCancelButton = FALSE,
        imageUrl = "",
        animation = TRUE,
      )
      return()
    }

    n_snps <- nrow(vcf@gt)

    #Extract GT
    geno <- extract.gt(vcf, element = "GT")
    geno <- apply(geno, 2, convert_to_dosage)
    class(geno) <- "numeric"

    return(list(geno, n_snps))
  } else {
    # If condition is met, show notification toast
    shinyalert(
      title = "Oops",
      text = "No valid genotype file detected",
      size = "xs",
      closeOnEsc = TRUE,
      closeOnClickOutside = FALSE,
      html = TRUE,
      type = "warning",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#004192",
      showCancelButton = FALSE,
      imageUrl = "",
      animation = TRUE,
    )

    return()
  }
}


