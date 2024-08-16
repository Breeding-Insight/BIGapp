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


test_that("Filtering without updog metrics",{

  #Variables
  filter_ploidy <- 2
  filter_maf <- 0.05
  size_depth <- 10
  snp_miss <- 100
  sample_miss <- 100
  OD_filter <- NULL
  Bias <- NULL
  Bias_min <- NULL
  Bias_max <- NULL
  Prop_mis <- 0.05
  maxpostprob_filter <- NULL
  max_post <- maxpostprob_filter
  output_name <- "out"
  snp_miss <- snp_miss/100
  sample_miss <- sample_miss/100
  ploidy <- filter_ploidy
  maf_filter <- filter_maf
  input$hist_bins <- 50

  input <- filtering_files <- list()
  input$updog_rdata$datapath <- system.file("vcf_example_out.vcf.gz", package = "BIGapp")

  temp_file <- tempfile(fileext = ".vcf.gz")

  #Input file
  vcf <- read.vcfR(input$updog_rdata$datapath, verbose = FALSE)
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

  if(length(vcf@gt) == 0) stop("All markers were filtered. Loose the parameters to access results in this tab.")

  #Getting missing data information
  #Add support for genotype matrix filtering?
  gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  filtering_files$snp_miss_df <- rowMeans(is.na(gt_matrix)) #SNP missing values
  filtering_files$sample_miss_df <- as.numeric(colMeans(is.na(gt_matrix))) #Sample missing values

  rm(gt_matrix) #Remove gt matrix

  #Writing file
  write.vcf(vcf, file = temp_file)

  #Get final_snps
  final_snps <- nrow(vcf)

  #export INFO dataframe
  filtering_files$raw_vcf_df

  # Apply the function to each row and bind the results into a new dataframe
  new_df <- data.frame(filtering_files$raw_vcf_df) %>%
    mutate(INFO_list = map(INFO, split_info_column)) %>%
    unnest_wider(INFO_list)

  #Save df to reactive value
  filtering_output <- list()
  filtering_output$df <- new_df

  ##Make plots

  #Missing data

    #Histogram
    hist(as.numeric(filtering_files$snp_miss_df),
         main = "Ratio of Missing Data per SNP After Filtering",
         xlab = "Proportion of Missing Data per SNP",
         ylab = "Number of SNPs",
         col = "lightblue",
         border = "black",
         xlim = c(0,1),
         breaks = as.numeric(input$hist_bins))
    axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

    # Add vertical lines
    abline(v = mean(as.numeric(filtering_files$snp_miss_df)), col = "red", lty = 2)  # Mean line
    abline(v = median(as.numeric(filtering_files$snp_miss_df)), col = "green", lty = 2)  # Median line
    abline(v = quantile(as.numeric(filtering_files$snp_miss_df), 0.95), col = "blue", lty = 2)
    legend("topright", legend=c("mean", "median", "quantile"),
           col=c("red", "green","blue"), lty=1:2, cex=0.8)

    #Histogram
    hist(as.numeric(filtering_files$sample_miss_df),
         main = "Ratio of Missing Data per Sample After Filtering",
         xlab = "Proportion of Missing Data per Sample",
         ylab = "Number of Samples",
         col = "lightblue",
         border = "black",
         xlim = c(0,1),
         breaks = as.numeric(input$hist_bins))
    axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

    # Add vertical lines
    abline(v = mean(as.numeric(filtering_files$sample_miss_df)), col = "red", lty = 2)  # Mean line
    abline(v = median(as.numeric(filtering_files$sample_miss_df)), col = "green", lty = 2)  # Median line
    abline(v = quantile(as.numeric(filtering_files$sample_miss_df), 0.95), col = "blue", lty = 2)
    legend("topright", legend=c("mean", "median", "quantile"),
           col=c("red", "green","blue"), lty=1:2, cex=0.8)


  ##Read Depth (I would prefer that this show the mean depth for SNPs or Samples instead of all loci/sample cells)
  quantile(as.numeric(new_df$DP), 0.95)


})
