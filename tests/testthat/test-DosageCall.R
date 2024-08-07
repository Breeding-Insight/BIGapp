context("Dosage Calling")

# Imports
library(vcfR)

test_that("Dosage Calling from MADC file",{

  madc_file <- system.file("iris_DArT_MADC.csv", package="BIGapp")
  output_name <- "output"
  ploidy <- 2
  cores <- 2
  model_select <- "norm" # TODO: test for other models

  # Status
  #Import genotype info if genotype matrix format

  # Call the get_counts function with the specified MADC file path and output file path
  #Status
  result_df <- get_counts(madc_file, output_name)

  #Call the get_matrices function
  matrices <- get_matrices(result_df)

  mout <- updog::multidog(refmat = matrices$ref_matrix,
                          sizemat = matrices$size_matrix,
                          ploidy = as.numeric(ploidy),
                          model = model_select,
                          nc = cores)

  expect_equal(sum(mout$snpdf$bias), 402.7979, tolerance = 0.01)
  expect_equal(sum(mout$inddf$postmean), 95229.13, tolerance = 0.01)

  # Convert updog to VCF
  updog2vcf(
    multidog.object = mout,
    output.file = output_name,
    updog_version = packageVersion("updog"),
    compress = TRUE
  )

  vcf_result <- read.vcfR(paste0(output_name,".vcf.gz"))

  DP <- sum(as.numeric(extract.gt(vcf_result, "DP")))

  expect_equal(DP, 23618990)

  MPP <- sum(as.numeric(extract.gt(vcf_result, "MPP")))

  expect_equal(MPP, 74519.94, tolerance = 0.01)
})


test_that("Dosage Calling from VCF file",{

  madc_file <- system.file("vcf_example_out.vcf.gz", package = "BIGapp")
  ploidy <- 2
  model_select <- "norm"
  cores <- 2
  output_name <- "out"

  #Initialize matrices list
  matrices <- list()

  #Import genotype information if in VCF format
  vcf <- read.vcfR(madc_file)

  #Get items in FORMAT column
  info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT

  info_ids <- extract_info_ids(info[1])
  chrom <- vcf@fix[,1]
  pos <- vcf@fix[,2]

  if (("DP" %in% info_ids) && (("RA" %in% info_ids) | ("AD" %in% info_ids))) {
    #Extract DP and RA and convert to matrices
    matrices$size_matrix <- extract.gt(vcf, element = "DP")
    if("RA" %in% info_ids){
      matrices$ref_matrix <- extract.gt(vcf, element = "RA")
    } else {
      ad_matrix <- extract.gt(vcf, element = "AD")
      matrices$ref_matrix <- matrix(sapply(strsplit(ad_matrix, ","), "[[", 1), nrow = nrow(matrices$size_matrix))
      colnames(matrices$ref_matrix) <- colnames(matrices$size_matrix)
    }

    class(matrices$size_matrix) <- "numeric"
    class(matrices$ref_matrix) <- "numeric"
    rownames(matrices$size_matrix) <- rownames(matrices$ref_matrix) <- paste0(chrom, "_", pos)

  }

  mout <- updog::multidog(refmat = matrices$ref_matrix,
                          sizemat = matrices$size_matrix,
                          ploidy = as.numeric(ploidy),
                          model = model_select,
                          nc = cores)

  expect_equal(sum(mout$snpdf$bias), 24.22855, tolerance = 0.01)
  expect_equal(sum(mout$inddf$postmean), 2755.687, tolerance = 0.01)

  # Convert updog to VCF
  updog2vcf(
    multidog.object = mout,
    output.file = output_name,
    updog_version = packageVersion("updog"),
    compress = TRUE
  )

  vcf_result <- read.vcfR(paste0(output_name,".vcf.gz"))

  DP <- sum(as.numeric(extract.gt(vcf_result, "DP")))
  expect_equal(DP, 18667)

  MPP <- sum(as.numeric(extract.gt(vcf_result, "MPP")))
  expect_equal(MPP, 2259.218, tolerance = 0.01)

})
