context("Dosage Calling")

# Imports
#library(vcfR)

test_that("Convert DArT files to VCF",{
  # Get the uploaded file paths
  dosage_file <- system.file("iris_DArT_Allele_Dose_Report.csv", package = "BIGapp")
  counts_file <- system.file("iris_DArT_Counts.csv", package = "BIGapp")
  ploidy <- 2

  # Use a temporary file path without appending .vcf
  temp_base <- tempfile()

  dosage_file_df <- read.csv(dosage_file)
  snp_number <- length(dosage_file_df$X.[-c(1:7)])

  # Convert to VCF using the BIGr package
  cat("Running BIGr::dosage2vcf...\n")
  dosage2vcf(
    dart.report = dosage_file,
    dart.counts = counts_file,
    output.file = temp_base,
    ploidy = as.numeric(ploidy)
  )

  # The output file should be temp_base.vcf
  output_name <- paste0(temp_base, ".vcf")

  # Check if the VCF file was created
  if (file.exists(output_name)) {
    cat("VCF file created successfully.\n")

    # Compress the VCF file using gzip
    gzip_file <- paste0(output_name, ".gz")
    gz <- gzfile(gzip_file, "w")
    writeLines(readLines(output_name), gz)
    close(gz)

    # Check if the gzip file was created
    if (file.exists(gzip_file)) {
      cat("Gzip file created successfully.\n")

      # Delete the temporary files
      unlink(gzip_file)
      unlink(output_name)

      cat("Temporary files deleted successfully.\n")
    } else {
      stop("Error: Failed to create the gzip file.")
    }
  } else {
    stop("Error: Failed to create the VCF file.")
  }
})
