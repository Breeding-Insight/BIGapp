context("GSAcc")

test_that("test Predictive Ability iris",{

  # packages
  library(vcfR)
  library(BIGapp)
  library(rrBLUP)
  library(dplyr)
  library(tidyr)
  library(ggplot2)

  # Inputs
  input <- list()

  input$trait_file$datapath <- system.file("iris_passport_file.csv", package = "BIGapp")
  input$pred_file$datapath <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")

  input$pred_color_select <- "red"
  input$pred_ploidy <- 2
  input$pred_trait_info <- "Petal.Length"
  input$pred_cv <- 5
  input$pred_fixed_info <- NULL
  input$pred_fixed_cat <- NULL
  input$pred_cores <- 3

  input$ped_file <- NULL

  fixed_traits <- input$pred_fixed_info

  #Close popup window when user "saves options"
  advanced_options <- list()
  advanced_options$ped_file <- input$ped_file

  ####Genomic Prediction Accuracy
  #This tab involved 3 observeEvents
  #1) to get the traits listed in the phenotype file
  #2) to input and validate the input files
  #3) to perform the genomic prediction

  #2) Error check for prediction and save input files
  continue_prediction <- NULL
  pred_inputs <- list(
    pheno_input = NULL,
    geno_input = NULL,
    pred_snps = NULL,
    pred_genos = NULL,
    pred_geno_pheno = NULL
  )

  pred_outputs <- list(
    corr_output = NULL,
    box_plot = NULL,
    violin_plot = NULL,
    comb_output = NULL,
    avg_GEBVs = NULL,
    all_GEBVs = NULL,
    colors = NULL
  )

  #Variables
  pheno <- read.csv(input$trait_file$datapath, header = TRUE, check.names = FALSE)
  row.names(pheno) <- pheno[,1]

  #Getting genotype matrix
  #Geno.file conversion if needed
  geno_snps <- BIGapp:::read_geno_file(input$pred_file$datapath, requires = "GT", ploidy = input$pred_ploidy, check = FALSE)
  geno <- geno_snps[[1]]
  pred_inputs$pred_snps <- geno_snps[[2]] # n markers
  pred_inputs$pred_genos <- ncol(geno) # n samples

  ##### input checks
  #Check that the ploidy entered is correct
  if (input$pred_ploidy != max(geno, na.rm = TRUE)) stop(paste0("The maximum value in the genotype file (",max(geno, na.rm = TRUE),") does not equal the ploidy entered"))

  #Make sure the trait file and genotype file are in the same order
  # Column names for geno (assuming these are the individual IDs)
  colnames_geno <- colnames(geno)
  # Assuming the first column in Pheno contains the matching IDs
  ids_pheno <- pheno[, 1]
  # Find common identifiers
  common_ids <- intersect(colnames_geno, ids_pheno)
  #Get number of id
  pred_inputs$pred_geno_pheno <- length(common_ids)

  #Throw an error if there are less matching samples in the phenotype file than the genotype file
  if (length(common_ids) == 0) {
    stop("All samples were missing from the phenotype file")
  } else {
    if (length(common_ids) < length(colnames_geno))
      warning(paste0((length(colnames_geno)-length(common_ids))," samples were removed for not having trait information"))
    if (length(common_ids) < length(ids_pheno))
      warning(paste0((length(ids_pheno)-length(common_ids))," samples were removed for not having genotypic information"))
  }

  # Subset and reorder geno and pheno to ensure they only contain and are ordered by common IDs
  geno_adj <- geno[, common_ids]  # Assuming that the columns can be directly indexed by IDs
  pheno <- pheno[match(common_ids, ids_pheno), ] # If there is pheno but not geno, the sample is also discarted

  # Check pedigree

  ## Make ouput as checked inputs pred_inputs
  pred_inputs$pheno_input <- pheno
  pred_inputs$geno_input <- geno_adj


  ####### rrBLUP
  ##Need to add ability for the use of parallelism for the for cross-validation
  ##Example at this tutorial also: https://www.youtube.com/watch?v=ARWjdQU6ays

  input$pred_model <- "rrBLUP"

  # Convert genotype matrix according to ploidy and model used
  geno_formated <- BIGapp:::format_geno_matrix(pred_inputs$geno_input,input$pred_model, input$pred_matrix, input$pred_ploidy)

  results <- BIGapp:::run_predictive_model(geno = geno_formated,
                                  pheno = pred_inputs$pheno_input,
                                  selected_traits = input$pred_trait_info,
                                  predictive_model = input$pred_model,
                                  relationship_matrix_type = input$pred_matrix,
                                  pedigree = pred_inputs$ped_input,
                                  fixed_effects = input$pred_fixed_info,
                                  categorical_fixed_effects = input$pred_fixed_cat,
                                  ploidy = input$pred_ploidy,
                                  cores = input$pred_cores,
                                  cycles = input$pred_cv,
                                  folds = 5)

  #Save to reactive value
  pred_outputs_rrBLUP <- pred_outputs
  pred_outputs_rrBLUP$corr_output <- results$PredictionAccuracy
  pred_outputs_rrBLUP$all_GEBVs <- results$GEBVs

  # Convert trait columns to numeric
  results$GEBVs <- results$GEBVs %>%
    mutate(across(all_of(input$pred_trait_info), ~ as.numeric(.x)))

  # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
  average_gebvs_df <- results$GEBVs %>%
    group_by(Sample) %>%
    summarize(across(all_of(input$pred_trait_info), mean, na.rm = TRUE))

  pred_outputs_rrBLUP$avg_GEBVs <- average_gebvs_df

  columns <- setdiff(colnames(results$PredictionAccuracy), c("Iter","Fold"))
  average_accuracy_df <- results$PredictionAccuracy %>%
    group_by(Iter) %>%
    summarize(across(all_of(columns), mean, na.rm = TRUE))

  pred_outputs_rrBLUP$comb_output <- average_accuracy_df

  #########

  ######### GBLUP
  #Note: should wrap the GBLUP into a function too
  # Define variables
  #train_size <- floor(percentage / 100 * total_population)
  #Cross validation number for progress bar (not involved in the calculations, just shiny visuals)

  # Convert genotype matrix according to ploidy
  input$pred_model <- "GBLUP"
  input$pred_matrix <- "Gmatrix"
  advanced_options$pred_matrix <- input$pred_matrix

  # Convert genotype matrix according to ploidy and model used
  geno_formated <- BIGapp:::format_geno_matrix(pred_inputs$geno_input,input$pred_model, input$pred_matrix, input$pred_ploidy)

  # Main function
  results <- BIGapp:::run_predictive_model(geno = geno_formated,
                                  pheno = pred_inputs$pheno_input,
                                  selected_traits = input$pred_trait_info,
                                  predictive_model = input$pred_model,
                                  relationship_matrix_type = input$pred_matrix,
                                  pedigree = pred_inputs$ped_input,
                                  fixed_effects = input$pred_fixed_info,
                                  categorical_fixed_effects = input$pred_fixed_cat,
                                  ploidy = input$pred_ploidy,
                                  cores = input$pred_cores,
                                  cycles = input$pred_cv,
                                  folds = 5)

  #Save to reactive value
  pred_outputs_gBLUP <- pred_outputs
  pred_outputs_gBLUP$corr_output <- results$PredictionAccuracy
  pred_outputs_gBLUP$all_GEBVs <- results$GEBVs

  # Convert trait columns to numeric
  GEBVs <- results$GEBVs %>%
    mutate(across(all_of(input$pred_trait_info), ~ as.numeric(.x)))

  # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
  average_gebvs_df <- GEBVs %>%
    group_by(Sample) %>%
    summarize(across(all_of(input$pred_trait_info), mean, na.rm = TRUE))

  pred_outputs_gBLUP$avg_GEBVs <- average_gebvs_df

  columns <- setdiff(colnames(results$PredictionAccuracy), c("Iter","Fold"))
  average_accuracy_df <- results$PredictionAccuracy %>%
    group_by(Iter) %>%
    summarize(across(all_of(columns), mean, na.rm = TRUE))

  pred_outputs_gBLUP$comb_output <- average_accuracy_df

  # Compare rrBLUP and GBLUP
  expect_equal(pred_outputs_gBLUP$corr_output[,1], pred_outputs_rrBLUP$corr_output[,1], tolerance = 0.02)
  expect_equal(pred_outputs_gBLUP$comb_output[,2], pred_outputs_rrBLUP$comb_output[,2], tolerance = 0.02)
  expect_equal(sum(pred_outputs_gBLUP$avg_GEBVs[,2]), -0.594, tolerance = 0.01)
  expect_equal(sum(as.numeric(pred_outputs_gBLUP$all_GEBVs[,1])), -2.971, tolerance = 0.1)
})

# test_that("test Predictive Ability wheat (BGLR dataset)",{
#
#   # packages
#   library(vcfR)
#   library(BIGapp)
#   library(rrBLUP)
#   library(dplyr)
#   library(tidyr)
#   library(ggplot2)
#
#   library(BGLR)
#
#   # Inputs
#   input <- list()
#
#   input$pred_color_select <- "red"
#   input$pred_ploidy <- 2
#   input$pred_trait_info <- "Pheno2"
#   input$pred_cv <- 5
#   input$pred_fixed_info <- NULL
#   input$pred_fixed_cat <- NULL
#   input$pred_cores <- 3
#
#   input$ped_file <- NULL
#
#   fixed_traits <- input$pred_fixed_info
#
#   #Close popup window when user "saves options"
#   advanced_options <- list()
#   advanced_options$ped_file <- input$ped_file
#
#   ####Genomic Prediction Accuracy
#   #This tab involved 3 observeEvents
#   #1) to get the traits listed in the phenotype file
#   #2) to input and validate the input files
#   #3) to perform the genomic prediction
#
#   #2) Error check for prediction and save input files
#   continue_prediction <- NULL
#   pred_inputs <- list(
#     pheno_input = NULL,
#     geno_input = NULL,
#     pred_snps = NULL,
#     pred_genos = NULL,
#     pred_geno_pheno = NULL
#   )
#
#   pred_outputs <- list(
#     corr_output = NULL,
#     box_plot = NULL,
#     violin_plot = NULL,
#     comb_output = NULL,
#     avg_GEBVs = NULL,
#     all_GEBVs = NULL,
#     colors = NULL
#   )
#
#   data(wheat)
#   #Variables
#   pheno <- wheat.Y
#   #row.names(pheno) <- pheno[,1]
#   colnames(pheno) <- paste0("Pheno", 1:4)
#   pheno <- data.frame(Sample_ID = rownames(pheno), pheno) # First column as the sample names is required
#
#   #Getting genotype matrix
#   #Geno.file conversion if needed
#   geno_raw <- wheat.X
#
#   # Codification 0 homozygous ref, 1 homozygous alt
#   geno <- matrix(as.numeric(gsub(1,2,geno_raw)), nrow = nrow(geno_raw))
#   colnames(geno) <- colnames(geno_raw)
#   rownames(geno) <- rownames(pheno)
#   geno <- t(geno)
#
#   pred_inputs$pred_snps <- nrow(geno) # n markers
#   pred_inputs$pred_genos <- ncol(geno) # n samples
#
#   # Update colors based on input
#   pred_outputs$colors <- assign_colors(input$pred_color_select)
#
#   ##### input checks
#   #Check that the ploidy entered is correct
#   if (input$pred_ploidy != max(geno, na.rm = TRUE)) stop(paste0("The maximum value in the genotype file (",max(geno, na.rm = TRUE),") does not equal the ploidy entered"))
#
#   #Make sure the trait file and genotype file are in the same order
#   # Column names for geno (assuming these are the individual IDs)
#   colnames_geno <- colnames(geno)
#   # Assuming the first column in Pheno contains the matching IDs
#   ids_pheno <- rownames(pheno)
#   # Find common identifiers
#   common_ids <- intersect(colnames_geno, ids_pheno)
#   #Get number of id
#   pred_inputs$pred_geno_pheno <- length(common_ids)
#
#   #Throw an error if there are less matching samples in the phenotype file than the genotype file
#   if (length(common_ids) == 0) {
#     stop("All samples were missing from the phenotype file")
#   } else {
#     if (length(common_ids) < length(colnames_geno))
#       warning(paste0((length(colnames_geno)-length(common_ids))," samples were removed for not having trait information"))
#     if (length(common_ids) < length(ids_pheno))
#       warning(paste0((length(ids_pheno)-length(common_ids))," samples were removed for not having genotypic information"))
#   }
#
#   # Subset and reorder geno and pheno to ensure they only contain and are ordered by common IDs
#   geno_adj <- geno[, common_ids]  # Assuming that the columns can be directly indexed by IDs
#   pheno <- pheno[match(common_ids, ids_pheno), ] # If there is pheno but not geno, the sample is also discarted
#
#   # Check pedigree
#
#   ## Make ouput as checked inputs pred_inputs
#   pred_inputs$pheno_input <- pheno
#   pred_inputs$geno_input <- geno_adj
#
#   ####### rrBLUP
#   ##Need to add ability for the use of parallelism for the for cross-validation
#   ##Example at this tutorial also: https://www.youtube.com/watch?v=ARWjdQU6ays
#
#   input$pred_model <- "rrBLUP"
#
#   # Convert genotype matrix according to ploidy and model used
#   geno_formated <- format_geno_matrix(pred_inputs$geno_input,input$pred_model, input$pred_matrix, input$pred_ploidy)
#
#   results <- run_predictive_model(geno = geno_formated,
#                                   pheno = pred_inputs$pheno_input,
#                                   selected_traits = input$pred_trait_info,
#                                   predictive_model = input$pred_model,
#                                   relationship_matrix_type = input$pred_matrix,
#                                   pedigree = pred_inputs$ped_input,
#                                   fixed_effects = input$pred_fixed_info,
#                                   categorical_fixed_effects = input$pred_fixed_cat,
#                                   ploidy = input$pred_ploidy,
#                                   cores = input$pred_cores,
#                                   cycles = input$pred_cv,
#                                   folds = 5)
#
#   #Save to reactive value
#   pred_outputs_rrBLUP <- pred_outputs
#   pred_outputs_rrBLUP$corr_output <- results$PredictionAccuracy
#   pred_outputs_rrBLUP$all_GEBVs <- results$GEBVs
#
#   # Convert trait columns to numeric
#   results$GEBVs <- results$GEBVs %>%
#     mutate(across(all_of(input$pred_trait_info), ~ as.numeric(.x)))
#
#   # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
#   average_gebvs_df <- results$GEBVs %>%
#     group_by(Sample) %>%
#     summarize(across(all_of(input$pred_trait_info), mean, na.rm = TRUE))
#
#   pred_outputs_rrBLUP$avg_GEBVs <- average_gebvs_df
#
#   columns <- setdiff(colnames(results$PredictionAccuracy), c("Iter","Fold"))
#   average_accuracy_df <- results$PredictionAccuracy %>%
#     group_by(Iter) %>%
#     summarize(across(all_of(columns), mean, na.rm = TRUE))
#
#   pred_outputs_rrBLUP$comb_output <- average_accuracy_df
#
#   #########
#
#   ######### GBLUP
#   #Note: should wrap the GBLUP into a function too
#   # Define variables
#   #train_size <- floor(percentage / 100 * total_population)
#   #Cross validation number for progress bar (not involved in the calculations, just shiny visuals)
#
#   # Convert genotype matrix according to ploidy
#   input$pred_model <- "GBLUP"
#   input$pred_matrix <- "Gmatrix"
#   advanced_options$pred_matrix <- input$pred_matrix
#
#   # Convert genotype matrix according to ploidy and model used
#   geno_formated <- format_geno_matrix(pred_inputs$geno_input,input$pred_model, input$pred_matrix, input$pred_ploidy)
#
#   # Main function
#   results <- run_predictive_model(geno = geno_formated,
#                                   pheno = pred_inputs$pheno_input,
#                                   selected_traits = input$pred_trait_info,
#                                   predictive_model = input$pred_model,
#                                   relationship_matrix_type = input$pred_matrix,
#                                   pedigree = pred_inputs$ped_input,
#                                   fixed_effects = input$pred_fixed_info,
#                                   categorical_fixed_effects = input$pred_fixed_cat,
#                                   ploidy = input$pred_ploidy,
#                                   cores = input$pred_cores,
#                                   cycles = input$pred_cv,
#                                   folds = 5)
#
#   #Save to reactive value
#   pred_outputs_gBLUP <- pred_outputs
#   pred_outputs_gBLUP$corr_output <- results$PredictionAccuracy
#   pred_outputs_gBLUP$all_GEBVs <- results$GEBVs
#
#   # Convert trait columns to numeric
#   GEBVs <- results$GEBVs %>%
#     mutate(across(all_of(input$pred_trait_info), ~ as.numeric(.x)))
#
#   # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
#   average_gebvs_df <- GEBVs %>%
#     group_by(Sample) %>%
#     summarize(across(all_of(input$pred_trait_info), mean, na.rm = TRUE))
#
#   pred_outputs_gBLUP$avg_GEBVs <- average_gebvs_df
#
#   columns <- setdiff(colnames(results$PredictionAccuracy), c("Iter","Fold"))
#   average_accuracy_df <- results$PredictionAccuracy %>%
#     group_by(Iter) %>%
#     summarize(across(all_of(columns), mean, na.rm = TRUE))
#
#   pred_outputs_gBLUP$comb_output <- average_accuracy_df
#
#   # Checks
#   expect_equal(pred_outputs_gBLUP$corr_output[,1], pred_outputs_rrBLUP$corr_output[,1], tolerance = 0.01)
#   expect_equal(pred_outputs_gBLUP$comb_output[,2], pred_outputs_rrBLUP$comb_output[,2], tolerance = 0.01)
#   expect_equal(sum(pred_outputs_gBLUP$corr_output[,1]), 12.44776, tolerance = 0.01)
#   expect_equal(sum(pred_outputs_gBLUP$comb_output[,2]$Pheno2), 2.489551, tolerance = 0.01)
#   expect_equal(sum(pred_outputs_gBLUP$avg_GEBVs[,2]$Pheno2), -2.913457, tolerance = 0.01)
#   expect_equal(sum(as.numeric(pred_outputs_gBLUP$all_GEBVs[,1])), -14.56729, tolerance = 0.1)
#
#   #### A matrix
#   input$pred_matrix <- "Amatrix"
#   advanced_options$pred_matrix <- input$pred_matrix
#
#   # Convert genotype matrix according to ploidy and model used
#   geno_formated <- format_geno_matrix(pred_inputs$geno_input,input$pred_model, input$pred_matrix, input$pred_ploidy)
#
#   # Main function
#   results <- run_predictive_model(geno = geno_formated,
#                                   pheno = pred_inputs$pheno_input,
#                                   selected_traits = input$pred_trait_info,
#                                   predictive_model = input$pred_model,
#                                   relationship_matrix_type = input$pred_matrix,
#                                   pedigree = pred_inputs$ped_input,
#                                   fixed_effects = input$pred_fixed_info,
#                                   categorical_fixed_effects = input$pred_fixed_cat,
#                                   ploidy = input$pred_ploidy,
#                                   cores = input$pred_cores,
#                                   cycles = input$pred_cv,
#                                   folds = 5,
#                                   relationship_matrix = wheat.A)
#
#   #Save to reactive value
#   pred_outputs_gBLUPA <- pred_outputs
#   pred_outputs_gBLUPA$corr_output <- results$PredictionAccuracy
#   pred_outputs_gBLUPA$all_GEBVs <- results$GEBVs
#
#   # Convert trait columns to numeric
#   GEBVs <- results$GEBVs %>%
#     mutate(across(all_of(input$pred_trait_info), ~ as.numeric(.x)))
#
#   # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
#   average_gebvs_df <- GEBVs %>%
#     group_by(Sample) %>%
#     summarize(across(all_of(input$pred_trait_info), mean, na.rm = TRUE))
#
#   pred_outputs_gBLUPA$avg_GEBVs <- average_gebvs_df
#
#   columns <- setdiff(colnames(results$PredictionAccuracy), c("Iter","Fold"))
#   average_accuracy_df <- results$PredictionAccuracy %>%
#     group_by(Iter) %>%
#     summarize(across(all_of(columns), mean, na.rm = TRUE))
#
#   pred_outputs_gBLUPA$comb_output <- average_accuracy_df
#
#   # Checks
#   expect_equal(sum(pred_outputs_gBLUPA$corr_output[,1]), 10.335, tolerance = 0.01)
#   expect_equal(sum(pred_outputs_gBLUPA$comb_output[,2]$Pheno2), 2.06, tolerance = 0.01)
#   expect_equal(sum(pred_outputs_gBLUPA$avg_GEBVs[,2]$Pheno2), 318.62, tolerance = 0.01)
#   expect_equal(sum(as.numeric(pred_outputs_gBLUPA$all_GEBVs[,1])), 1593.09, tolerance = 0.1)
#
# })
#
#
# test_that("test Predictive Ability Josue",{
#
#   # packages
#   library(vcfR)
#   library(BIGapp)
#   library(rrBLUP)
#   library(dplyr)
#   library(tidyr)
#   library(ggplot2)
#
#   # Inputs
#   input <- list()
#
#   #input$trait_file$datapath <- "BIG_pheno2.csv"
#   input$trait_file$datapath <- "BIG_phenos.csv"
#   input$pred_file$datapath <- "BIG_genos.vcf"
#   #input$ped_file$datapath <- "sealice_ped.csv"
#   input$ped_file$datapath <- "BIG_ped.csv"
#
#   input$pred_color_select <- "red"
#   input$pred_ploidy <- 2
#   input$pred_trait_info <- "licedensity"
#   input$pred_cv <- 5
#   input$pred_fixed_info <- NULL
#   input$pred_fixed_cat <- NULL
#   input$pred_cores <- 3
#
#   fixed_traits <- input$pred_fixed_info
#
#   #Close popup window when user "saves options"
#   advanced_options <- list()
#   advanced_options$ped_file <- input$ped_file
#
#   ####Genomic Prediction Accuracy
#   #This tab involved 3 observeEvents
#   #1) to get the traits listed in the phenotype file
#   #2) to input and validate the input files
#   #3) to perform the genomic prediction
#
#   #2) Error check for prediction and save input files
#   continue_prediction <- NULL
#   pred_inputs <- list(
#     pheno_input = NULL,
#     geno_input = NULL,
#     pred_snps = NULL,
#     pred_genos = NULL,
#     pred_geno_pheno = NULL,
#     ped_input = NULL
#   )
#
#   pred_outputs <- list(
#     corr_output = NULL,
#     box_plot = NULL,
#     violin_plot = NULL,
#     comb_output = NULL,
#     avg_GEBVs = NULL,
#     all_GEBVs = NULL,
#     colors = NULL
#   )
#
#   #Variables
#   pheno <- read.csv(input$trait_file$datapath, header = TRUE, check.names = FALSE)
#   #pheno <- pheno[,-c(1:(which(colnames(pheno) == "Sample_ID") - 1))] # Sample_ID must be the first column
#   rownames(pheno) <- pheno[,1]
#
#   #Getting genotype matrix
#   #Geno.file conversion if needed
#   geno_snps <- read_geno_file(input$pred_file$datapath, requires = "GT")
#   geno <- geno_snps[[1]]
#   pred_inputs$pred_snps <- geno_snps[[2]] # n markers
#   pred_inputs$pred_genos <- ncol(geno) # n samples
#
#   ##### input checks
#   #Check that the ploidy entered is correct
#   if (input$pred_ploidy != max(geno, na.rm = TRUE)) stop(paste0("The maximum value in the genotype file (",max(geno, na.rm = TRUE),") does not equal the ploidy entered"))
#
#   #Make sure the trait file and genotype file are in the same order
#   # Column names for geno (assuming these are the individual IDs)
#   colnames_geno <- colnames(geno)
#   # Assuming the first column in Pheno contains the matching IDs
#   ids_pheno <- pheno[, 1]
#   # Find common identifiers
#   common_ids <- intersect(colnames_geno, ids_pheno)
#   #Get number of id
#   pred_inputs$pred_geno_pheno <- length(common_ids)
#
#   #Throw an error if there are less matching samples in the phenotype file than the genotype file
#   if (length(common_ids) == 0) {
#     stop("All samples were missing from the phenotype file")
#   } else {
#     if (length(common_ids) < length(colnames_geno))
#       warning(paste0((length(colnames_geno)-length(common_ids))," samples were removed for not having trait information"))
#     if (length(common_ids) < length(ids_pheno))
#       warning(paste0((length(ids_pheno)-length(common_ids))," samples were removed for not having genotypic information"))
#   }
#
#   # Subset and reorder geno and pheno to ensure they only contain and are ordered by common IDs
#   geno_adj <- geno[, common_ids]  # Assuming that the columns can be directly indexed by IDs
#   pheno <- pheno[match(common_ids, ids_pheno), ] # If there is pheno but not geno, the sample is also discarted
#
#   # Check pedigree
#   #Import pedigree file, where pedigree data name (3-column way format). Unknown value should be equal 0
#   if(!is.null(advanced_options$ped_file$datapath)){
#     ped <- read.csv(advanced_options$ped_file$datapath, check.names = FALSE, colClasses = "factor")
#     colnames(ped) <- c("Ind", "P1", "P2")
#     #Convert NAs to 0
#     ped[is.na(ped)] <- 0
#
#     common_ped <- intersect(ped$Ind, pheno[,1])
#     #Throw an error if there are less matching samples in the phenotype file than the pedigree file
#     if (length(common_ped) == 0) {
#       stop("All samples were missing from the phenotype file")
#     } else {
#       rm_unr <- remove_unrelated(ped, samples_with_trait_info = pheno[,1])
#       extended_ped <- rm_unr[[1]]
#       gen <- rm_unr[[2]]
#       cat(paste0("You have pedigree information until the ", gen,"th generation\n"))
#
#       if (length(common_ped) < length(pheno[,1])){
#         warning(paste0((length(pheno[,1])-length(common_ped))," samples were removed from the phenotype data for not having pedigree information"))
#         if(length(which(!pheno[,1] %in% extended_ped$Ind)) > 0) pheno <- pheno[-which(!pheno[,1] %in% extended_ped$Ind),]
#         if(length(which(!colnames(geno_adj) %in% extended_ped$Ind)) > 0) geno_adj <- geno_adj[,-which(!colnames(geno_adj) %in% extended_ped$Ind)]
#       }
#       if (length(ped$Ind) > length(extended_ped$Ind))
#         warning(paste0((length(ped$Ind)-length(extended_ped$Ind))," samples in the pedigree file were unrelated to the samples with phenotype information. They were removed from the analysis.")) # samples not removed
#
#       ped_temp <- tempfile()
#       ped_temp_file <- extended_ped
#       colnames(ped_temp_file) <- c("id", "sire", "dam")
#       write.table(ped_temp_file, file = ped_temp)
#       ped_check <- BIGr::check_ped(ped_temp)
#       if(dim(ped_check$repeated_ids)[1] != 0) stop("Check for repeated IDs in the pedigree file")
#       if(dim(ped_check$messy_parents)[1] != 0) stop(paste("We found inconsistencies in the pedigree file for the individuals:", paste0(ped_check$messy_parents$id, collapse = ", ")))
#     }
#
#   }
#
#   ## Make ouput as checked inputs pred_inputs
#   pred_inputs$pheno_input <- pheno
#   pred_inputs$geno_input <- geno_adj
#   pred_inputs$ped_input <- extended_ped
#
#   ####### rrBLUP
#   ##Need to add ability for the use of parallelism for the for cross-validation
#   ##Example at this tutorial also: https://www.youtube.com/watch?v=ARWjdQU6ays
#
#   input$pred_model <- "rrBLUP"
#
#   # Convert genotype matrix according to ploidy and model used
#   geno_formated <- format_geno_matrix(pred_inputs$geno_input,input$pred_model, input$pred_matrix, input$pred_ploidy)
#
#   results <- run_predictive_model(geno = geno_formated,
#                                   pheno = pred_inputs$pheno_input,
#                                   selected_traits = input$pred_trait_info,
#                                   predictive_model = input$pred_model,
#                                   relationship_matrix_type = input$pred_matrix,
#                                   pedigree = pred_inputs$ped_input,
#                                   fixed_effects = input$pred_fixed_info,
#                                   categorical_fixed_effects = input$pred_fixed_cat,
#                                   ploidy = input$pred_ploidy,
#                                   cores = input$pred_cores,
#                                   cycles = input$pred_cv,
#                                   folds = 5)
#
#   #Save to reactive value
#   pred_outputs_rrBLUP <- pred_outputs
#   pred_outputs_rrBLUP$corr_output <- results$PredictionAccuracy
#   pred_outputs_rrBLUP$all_GEBVs <- results$GEBVs
#
#   # Convert trait columns to numeric
#   results$GEBVs <- results$GEBVs %>%
#     mutate(across(all_of(input$pred_trait_info), ~ as.numeric(.x)))
#
#   # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
#   average_gebvs_df <- results$GEBVs %>%
#     group_by(Sample) %>%
#     summarize(across(all_of(input$pred_trait_info), mean, na.rm = TRUE))
#
#   pred_outputs_rrBLUP$avg_GEBVs <- average_gebvs_df
#
#   columns <- setdiff(colnames(results$PredictionAccuracy), c("Iter","Fold"))
#   average_accuracy_df <- results$PredictionAccuracy %>%
#     group_by(Iter) %>%
#     summarize(across(all_of(columns), mean, na.rm = TRUE))
#
#   pred_outputs_rrBLUP$comb_output <- average_accuracy_df
#
#   # Checks
#   expect_equal(mean(pred_outputs_rrBLUP$corr_output[,1]), 0.1848, tolerance = 0.01)
#   expect_equal(mean(pred_outputs_rrBLUP$comb_output[,2]$licedensity), 0.1848, tolerance = 0.01)
#   expect_equal(sum(pred_outputs_rrBLUP$avg_GEBVs[,2]$licedensity), 19.378, tolerance = 0.01)
#   expect_equal(sum(as.numeric(pred_outputs_rrBLUP$all_GEBVs[,1])), 96.89, tolerance = 0.1)
#   #########
#
#   ######### GBLUP
#   #Note: should wrap the GBLUP into a function too
#   # Define variables
#   #train_size <- floor(percentage / 100 * total_population)
#   #Cross validation number for progress bar (not involved in the calculations, just shiny visuals)
#
#   # Convert genotype matrix according to ploidy
#   input$pred_model <- "GBLUP"
#   input$pred_matrix <- "Gmatrix"
#   advanced_options$pred_matrix <- input$pred_matrix
#
#   # Convert genotype matrix according to ploidy and model used
#   geno_formated <- format_geno_matrix(pred_inputs$geno_input,input$pred_model, input$pred_matrix, input$pred_ploidy)
#
#   # Main function
#   results <- run_predictive_model(geno = geno_formated,
#                                   pheno = pred_inputs$pheno_input,
#                                   selected_traits = input$pred_trait_info,
#                                   predictive_model = input$pred_model,
#                                   relationship_matrix_type = input$pred_matrix,
#                                   pedigree = pred_inputs$ped_input,
#                                   fixed_effects = input$pred_fixed_info,
#                                   categorical_fixed_effects = input$pred_fixed_cat,
#                                   ploidy = input$pred_ploidy,
#                                   cores = input$pred_cores,
#                                   cycles = input$pred_cv,
#                                   folds = 5)
#
#   #Save to reactive value
#   pred_outputs_gBLUP <- pred_outputs
#   pred_outputs_gBLUP$corr_output <- results$PredictionAccuracy
#   pred_outputs_gBLUP$all_GEBVs <- results$GEBVs
#
#   # Convert trait columns to numeric
#   GEBVs <- results$GEBVs %>%
#     mutate(across(all_of(input$pred_trait_info), ~ as.numeric(.x)))
#
#   # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
#   average_gebvs_df <- GEBVs %>%
#     group_by(Sample) %>%
#     summarize(across(all_of(input$pred_trait_info), mean, na.rm = TRUE))
#
#   pred_outputs_gBLUP$avg_GEBVs <- average_gebvs_df
#
#   columns <- setdiff(colnames(results$PredictionAccuracy), c("Iter","Fold"))
#   average_accuracy_df <- results$PredictionAccuracy %>%
#     group_by(Iter) %>%
#     summarize(across(all_of(columns), mean, na.rm = TRUE))
#
#   pred_outputs_gBLUP$comb_output <- average_accuracy_df
#
#   # Compare rrBLUP and GBLUP
#   expect_equal(pred_outputs_gBLUP$corr_output[,1], pred_outputs_rrBLUP$corr_output[,1], tolerance = 0.01)
#   expect_equal(pred_outputs_gBLUP$comb_output[,2], pred_outputs_rrBLUP$comb_output[,2], tolerance = 0.01)
#   expect_equal(mean(pred_outputs_gBLUP$corr_output[,1]), 0.1848, tolerance = 0.01)
#   expect_equal(mean(pred_outputs_gBLUP$comb_output[,2]$Pheno3), 0.1848, tolerance = 0.01)
#   expect_equal(sum(pred_outputs_gBLUP$avg_GEBVs[,2]$Pheno3), -0.8889, tolerance = 0.01)
#   expect_equal(sum(as.numeric(pred_outputs_gBLUP$all_GEBVs[,1])), -4.444, tolerance = 0.1)
#
#   ## Using A matrix
#   # Convert genotype matrix according to ploidy
#   input$pred_model <- "GBLUP"
#   input$pred_matrix <- "Amatrix"
#   advanced_options$pred_matrix <- input$pred_matrix
#
#   # Convert genotype matrix according to ploidy and model used
#   geno_formated <- format_geno_matrix(pred_inputs$geno_input,input$pred_model, input$pred_matrix, input$pred_ploidy)
#
#   # Main function
#   results <- run_predictive_model(geno = geno_formated,
#                                   pheno = pred_inputs$pheno_input,
#                                   selected_traits = input$pred_trait_info,
#                                   predictive_model = input$pred_model,
#                                   relationship_matrix_type = input$pred_matrix,
#                                   pedigree = pred_inputs$ped_input,
#                                   fixed_effects = input$pred_fixed_info,
#                                   categorical_fixed_effects = input$pred_fixed_cat,
#                                   ploidy = input$pred_ploidy,
#                                   cores = input$pred_cores,
#                                   cycles = input$pred_cv,
#                                   folds = 5)
#
#   str(results)
#   #Save to reactive value
#   pred_outputs_gBLUPA <- pred_outputs
#   pred_outputs_gBLUPA$corr_output <- results$PredictionAccuracy
#   pred_outputs_gBLUPA$all_GEBVs <- results$GEBVs
#
#   # Convert trait columns to numeric
#   GEBVs <- results$GEBVs %>%
#     mutate(across(all_of(input$pred_trait_info), ~ as.numeric(.x)))
#
#   # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
#   average_gebvs_df <- GEBVs %>%
#     group_by(Sample) %>%
#     summarize(across(all_of(input$pred_trait_info), mean, na.rm = TRUE))
#
#   pred_outputs_gBLUPA$avg_GEBVs <- average_gebvs_df
#
#   columns <- setdiff(colnames(results$PredictionAccuracy), c("Iter","Fold"))
#   average_accuracy_df <- results$PredictionAccuracy %>%
#     group_by(Iter) %>%
#     summarize(across(all_of(columns), mean, na.rm = TRUE))
#
#   pred_outputs_gBLUPA$comb_output <- average_accuracy_df
#
#   # Checks
#   expect_equal(mean(pred_outputs_gBLUPA$corr_output[,1]), 0.153, tolerance = 0.01)
#   expect_equal(mean(pred_outputs_gBLUPA$comb_output[,2]$Pheno3), 0.153, tolerance = 0.01)
#   expect_equal(sum(pred_outputs_gBLUPA$avg_GEBVs[,2]$Pheno3), 2.165, tolerance = 0.01)
#   expect_equal(sum(as.numeric(pred_outputs_gBLUPA$all_GEBVs[,1])), 10.826, tolerance = 0.1)
#
#   ## Using H matrix
#   # Convert genotype matrix according to ploidy
#   input$pred_model <- "GBLUP"
#   input$pred_matrix <- "Hmatrix"
#   advanced_options$pred_matrix <- input$pred_matrix
#
#   # Convert genotype matrix according to ploidy and model used
#   geno_formated <- format_geno_matrix(pred_inputs$geno_input,input$pred_model,input$pred_matrix, input$pred_ploidy)
#
#   # Main function
#   results <- run_predictive_model(geno = geno_formated,
#                                   pheno = pred_inputs$pheno_input,
#                                   selected_traits = input$pred_trait_info,
#                                   predictive_model = input$pred_model,
#                                   relationship_matrix_type = input$pred_matrix,
#                                   pedigree = pred_inputs$ped_input,
#                                   fixed_effects = input$pred_fixed_info,
#                                   categorical_fixed_effects = input$pred_fixed_cat,
#                                   ploidy = input$pred_ploidy,
#                                   cores = input$pred_cores,
#                                   cycles = input$pred_cv,
#                                   folds = 5)
#
#   #Save to reactive value
#   pred_outputs_gBLUPH <- pred_outputs
#   pred_outputs_gBLUPH$corr_output <- results$PredictionAccuracy
#   pred_outputs_gBLUPH$all_GEBVs <- results$GEBVs
#
#   # Convert trait columns to numeric
#   GEBVs <- results$GEBVs %>%
#     mutate(across(all_of(input$pred_trait_info), ~ as.numeric(.x)))
#
#   # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
#   average_gebvs_df <- GEBVs %>%
#     group_by(Sample) %>%
#     summarize(across(all_of(input$pred_trait_info), mean, na.rm = TRUE))
#
#   pred_outputs_gBLUPH$avg_GEBVs <- average_gebvs_df
#
#   columns <- setdiff(colnames(results$PredictionAccuracy), c("Iter","Fold"))
#   average_accuracy_df <- results$PredictionAccuracy %>%
#     group_by(Iter) %>%
#     summarize(across(all_of(columns), mean, na.rm = TRUE))
#
#   pred_outputs_gBLUPH$comb_output <- average_accuracy_df
#
#   # Checks
#   expect_equal(mean(pred_outputs_gBLUPH$corr_output[,1]), 0.195, tolerance = 0.01)
#   expect_equal(mean(pred_outputs_gBLUPH$comb_output[,2]$Pheno3), 0.195, tolerance = 0.01)
#   expect_equal(sum(pred_outputs_gBLUPH$avg_GEBVs[,2]$Pheno3), -0.187, tolerance = 0.01)
#   expect_equal(sum(as.numeric(pred_outputs_gBLUPH$all_GEBVs[,1])), -0.934, tolerance = 0.1)
# })
