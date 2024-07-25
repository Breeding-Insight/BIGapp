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

  # Save the csv file for review and use in R
  #df_name <- paste0(output_name,'_MADC_alt_ref_counts.csv')

  #write.csv(filtered_df, file = df_name, row.names = FALSE)
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
