context("Filtering")

#library(vcfR)
#library(BIGr)
#library(testthat)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)

test_that("Filtering with updog metrics",{

  #Variables
  filter_ploidy <- 2
  filter_maf <- 0.05
  size_depth <- 10
  snp_miss <- 50
  sample_miss <- 50
  OD_filter <- 0.05
  Bias <- c(0.5, 2)
  Bias_min <- Bias[1]
  Bias_max <- Bias[2]
  Prop_mis <- 0.05
  maxpostprob_filter <- 0.5
  max_post <- maxpostprob_filter
  output_name <- "out"
  snp_miss <- snp_miss/100
  sample_miss <- sample_miss/100
  ploidy <- filter_ploidy
  maf_filter <- filter_maf

  input <- filtering_files <- list()
  input$updog_rdata$datapath <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")

  temp_file <- tempfile(fileext = ".vcf.gz")

  #Input file
  vcf <- read.vcfR(input$updog_rdata$datapath, verbose = FALSE)

  # Identify if have updog parameters
  format_fields <- unique(vcf@gt[,1])
  info_fields <- vcf@fix[1,8]

  updog_par <- grepl("MPP", format_fields) & grepl("PMC", info_fields) & grepl("BIAS", info_fields)

  #Starting SNPs
  starting_snps <- nrow(vcf)
  #export INFO dataframe
  filtering_files$raw_vcf_df <- data.frame(vcf@fix)

  #Filtering
  vcf <- filterVCF(vcf.file = vcf,
                   ploidy=ploidy,
                   output.file=NULL,
                   filter.OD = OD_filter,
                   filter.BIAS.min = Bias_min,
                   filter.BIAS.max = Bias_max,
                   filter.DP = as.numeric(size_depth),
                   filter.PMC = Prop_mis,
                   filter.SAMPLE.miss = as.numeric(sample_miss),
                   filter.SNP.miss = as.numeric(snp_miss),
                   filter.MAF = as.numeric(maf_filter),
                   filter.MPP = max_post)

  #Getting missing data information
  #Add support for genotype matrix filtering?
  gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  filtering_files$snp_miss_df <- rowMeans(is.na(gt_matrix)) #SNP missing values
  filtering_files$sample_miss_df <- as.numeric(colMeans(is.na(gt_matrix))) #Sample missing values

  expect_true(all(table(gt_matrix[,10]) == c(20,13,8)))

  rm(gt_matrix) #Remove gt matrix

  #Writing file
  write.vcf(vcf, file = temp_file)

  #Get final_snps
  final_snps <- nrow(vcf)
  expect_equal(final_snps, 43)

})

