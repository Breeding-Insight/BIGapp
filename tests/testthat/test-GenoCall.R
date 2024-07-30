context("Dosage Calling")

test_that("Dosage Calling",{

  # Imports
  library(vcfR)

  # Convert MADC to VCF - Couldn't find a example file to test it (ask Alex)
  # TODO


  # Run updog
  madc_file <- system.file("vcf_example_out.vcf.gz", package = "BIGapp")
  #madc_file <- "inst/data/vcf_example_out.vcf.gz"
  ploidy <- 2
  model_select <- "norm"
  cores <- 2
  output_name <- "out.vcf"

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

  # Convert updog to VCF
  BIGr::updog2vcf(
    multidog.object = mout,
    ploidy = ploidy,
    output.file = output_name
  )

})
