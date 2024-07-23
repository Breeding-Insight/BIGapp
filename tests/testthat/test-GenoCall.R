context("Dosage Calling")

test_that("all functions",{
  
  # Convert MADC to VCF - Couldn't find a example file to test it (ask Alex)
  dosage_file
  counts_file
  temp_base
  ploidy
  
  BIGr::dosage2vcf(
    dart.report = dosage_file,
    dart.counts = counts_file,
    output.file = temp_base,
    ploidy = as.numeric(ploidy)
  )
  
  # Run updog
  #madc_file <- system.file("data/vcf_sample_out.vcf.gz", package = "BIGapp")
  madc_file <- "inst/data/vcf_example_out.vcf.gz"
  ploidy <- 2
  model_select <- "norm"
  cores <- 2
  
  #Initialize matrices list
  matrices <- list()
  
  #Import genotype information if in VCF format
  vcf <- read.vcfR(madc_file)
  
  #Get items in FORMAT column
  info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT
  
  info_ids <- extract_info_ids(info[1])
  
  if (("DP" %in% info_ids) && (("RA" %in% info_ids) | ("AD" %in% info_ids))) {
    #Extract DP and RA and convert to matrices
    matrices$size_matrix <- extract.gt(vcf, element = "DP")
    if("RA" %in% info_ids){
      matrices$ref_matrix <- extract.gt(vcf, element = "RA")
    } else {
      ad_matrix <- extract.gt(vcf, element = "AD")
      
    }
    class(matrices$size_matrix) <- "numeric"
    class(matrices$ref_matrix) <- "numeric"
    rm(vcf) #Remove VCF
    
    snp_number <- (nrow(matrices$size_matrix) / 2)
    
  }else{
    ##Add user warning about read depth and allele read depth not found
    warning("Error: DP and RA FORMAT flags not found in VCF file")
    return()
  }
  
  
  mout <- multidog(refmat = matrices$ref_matrix, 
                   sizemat = matrices$size_matrix, 
                   ploidy = as.numeric(ploidy),  
                   model = model_select,
                   nc = cores)
  
  # Convert updog to VCF
  
})