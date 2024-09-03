context("GSAcc")

test_that("test Predictive Ability",{

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
  input$pred_color_select <- "red"
  input$pred_ploidy <- 2
  input$pred_file$datapath <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")

  input$pred_trait_info <- "Petal.Length"
  input$pred_cv <- 5

  input$pred_fixed_info <- NULL
  input$pred_fixed_cat <- NULL
  input$pred_cores <- 3

  input$pred_model <- "rrBLUP"
  input$pred_matrix <- "Gmatrix"
  input$ped_file <- NULL

  #Close popup window when user "saves options"
  advanced_options <- list()
  advanced_options$pred_model <- input$pred_model
  advanced_options$pred_matrix <- input$pred_matrix
  advanced_options$ped_file <- input$ped_file

  ####Genomic Prediction Accuracy
  #This tab involved 3 observeEvents
  #1) to get the traits listed in the phenotype file
  #2) to input and validate the input files
  #3) to perform the genomic prediction

  #1) Get traits
  info_df <- read.csv(input$trait_file$datapath, header = TRUE, check.names = FALSE, nrow = 0)
  trait_var <- colnames(info_df)
  trait_var <- trait_var[2:length(trait_var)]

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

  # Update colors based on input
  pred_outputs$colors <- switch(input$pred_color_select,
                                "red" = "#F8766D",
                                "blue" = "#00BFC4",
                                "green" = "#00BA38",
                                input$pred_color_select)


  #Variables
  ploidy <- as.numeric(input$pred_ploidy)
  geno_path <- input$pred_file$datapath
  pheno <- read.csv(input$trait_file$datapath, header = TRUE, check.names = FALSE)
  row.names(pheno) <- pheno[,1]
  traits <- input$pred_trait_info
  CVs <- as.numeric(input$pred_cv)

  #Getting genotype matrix

  #Geno.file conversion if needed
  geno_snps <- read_geno_file(geno_path, requires = "GT")
  geno <- geno_snps[[1]]
  pred_inputs$pred_snps <- geno_snps[[2]]

  #Save number of samples in file
  pred_inputs$pred_genos <- ncol(geno)

  #Check that the ploidy entered is correct
  if (ploidy != max(geno, na.rm = TRUE)) stop(paste0("The maximum value in the genotype file (",max(geno, na.rm = TRUE),") does not equal the ploidy entered"))

  # Convert genotype matrix according to ploidy
  geno_adj_init <- 2 * (geno / as.numeric(ploidy)) - 1

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
    if (length(common_ids) < length(colnames_geno))  stop(paste0((length(colnames_geno)-length(common_ids))," samples were removed for not having trait information"))
  }

  # Subset and reorder geno and pheno to ensure they only contain and are ordered by common IDs
  geno_adj <- geno_adj_init[, common_ids]  # Assuming that the columns can be directly indexed by IDs
  pheno <- pheno[match(common_ids, ids_pheno), ]

  ##Save to reactive values
  #Gmatrix needs the original allele count values, so the user matrix selection determines the genotype matrix used
  pred_inputs$pheno_input <- pheno
  if (advanced_options$pred_matrix == "Gmatrix" || is.null(advanced_options$pred_matrix)) {
    pred_inputs$geno_input <- geno_adj
  } else if (advanced_options$pred_matrix == "Hmatrix") {
    pred_inputs$geno_input <- geno[, common_ids]
  } else {
    pred_inputs$geno_input <- geno_adj
  }

  #3) Analysis only proceeds once continue_prediction is converted to TRUE
  #Variables
  ploidy <- as.numeric(input$pred_ploidy)
  geno_adj <- pred_inputs$geno_input
  pheno <- pred_inputs$pheno_input
  traits <- input$pred_trait_info
  CVs <- as.numeric(input$pred_cv)
  fixed_traits <- input$pred_fixed_info
  fixed_cat <- input$pred_fixed_cat
  fixed_cov <- if (is.null(input$pred_fixed_info) || length(input$pred_fixed_info) == length(input$pred_fixed_cat)) {
    NULL
  } else {
    setdiff(input$pred_fixed_info, input$pred_fixed_cat)
  }
  cores <- input$pred_cores

  #Assign colors
  if (input$pred_color_select == "red"){
    pred_outputs$colors <- "#F8766D"
  } else if (input$pred_color_select == "blue") {
    pred_outputs$colors <- "#00BFC4"
  } else if (input$pred_color_select == "green") {
    pred_outputs$colors <- "#00BA38"
  } else{
    pred_outputs$colors <- input$pred_color_select
  }

  #Control whether rrBLUP or GBLUP run depending on user input

  ####### rrBLUP
  ##Need to add ability for the use of parallelism for the for cross-validation
  ##Example at this tutorial also: https://www.youtube.com/watch?v=ARWjdQU6ays

  pred_inputs$geno_input <- geno_adj

  # Example call to the function
  #This is slow when using 3k markers and 1.2k samples...will need to parallelize if using this script...
  results <- genomic_prediction(geno = geno_adj,
                                pheno = pheno,
                                traits = traits,
                                fixed_effects = fixed_traits,
                                iters = input$pred_cv,
                                cores = cores)

  #With fixed effects (need to inforporate the ability for fixed effects into the prediction?)
  #results <- genomic_prediction(geno_matrix, phenotype_df, c("height", "weight"), "~ age + sex")

  #Save to reactive value
  pred_outputs_rrBLUP <- pred_outputs
  pred_outputs_rrBLUP$corr_output <- results$PredictionAccuracy
  pred_outputs_rrBLUP$all_GEBVs <- results$GEBVs

  # Convert trait columns to numeric
  results$GEBVs <- results$GEBVs %>%
    mutate(across(all_of(traits), ~ as.numeric(.x)))

  # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
  average_gebvs_df <- results$GEBVs %>%
    group_by(Sample) %>%
    summarize(across(all_of(traits), mean, na.rm = TRUE))

  pred_outputs_rrBLUP$avg_GEBVs <- average_gebvs_df

  columns <- setdiff(colnames(results$PredictionAccuracy), c("Iter","Fold"))
  average_accuracy_df <- results$PredictionAccuracy %>%
    group_by(Iter) %>%
    summarize(across(all_of(columns), mean, na.rm = TRUE))

  pred_outputs_rrBLUP$comb_output <- average_accuracy_df

  df <- pred_outputs_rrBLUP$corr_output
  df <- df %>% dplyr::select(-Fold, -Iter)

  #Probably want to add the ability for the user to select which trait(s) to display here

  #Convert to long format for ggplot
  df_long <- pivot_longer(
    df,
    cols = colnames(df),  # Exclude the Cycle column from transformation
    names_to = "Trait",  # New column for trait names
    values_to = "Correlation"  # New column for correlation values
  )

  #This can be adapted if we start comparing more than one GP model
  #Also consider a violin plot to show each cor value
  #plot <- ggplot(df_long, aes(x = factor(Trait), y = Correlation, fill = "red"), fill = "red") +
  pred_outputs_rrBLUP$box_plot <- ggplot(df_long, aes(x = "rrBLUP", y = Correlation, fill = "red"), fill = "red") +
    #geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 0.7, outlier.size = 0.2) +
    geom_boxplot() +
    facet_wrap(~ Trait, nrow = 1) +  # Facet by trait, allowing different y-scales
    labs(title = "Predictive Ability by Trait",
         x = " ",
         y = "Pearson Correlation") +
    #theme_minimal() +                      # Using a minimal theme
    theme(legend.position = "none",
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2),
          strip.text.x = element_text(face = "bold"),
          axis.text.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank())

  pred_outputs_rrBLUP$violin_plot <- ggplot(df_long, aes(x = "rrBLUP", y = Correlation, fill = "red")) +
    geom_violin(trim = TRUE) +  # Add violin plot
    geom_point(position = position_jitter(width = 0.1), color = "black", size = 1.5) +  # Add jittered points
    facet_wrap(~ Trait, nrow = 1) +  # Facet by trait, allowing different y-scales
    labs(title = "Predictive Ability by Trait",
         x = " ",  # x-label is blank because it's not relevant per facet
         y = "Pearson Correlation") +
    theme(legend.position = "none",
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2),
          strip.text.x = element_text(face = "bold"),
          axis.text.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank())

  plots <- list(pred_outputs_rrBLUP$box_plot, pred_outputs_rrBLUP$violin_plot)

  #Output the genomic prediction correlation box plots
  plots[[1]]  + scale_fill_manual(values = pred_outputs_rrBLUP$colors)

  #Output the genomic prediction correlation box plots
  plots[[2]] + scale_fill_manual(values = pred_outputs_rrBLUP$colors)

  #Output the prediction tables
  pred_outputs_rrBLUP$corr_output
  pred_outputs_rrBLUP$comb_output
  pred_outputs_rrBLUP$all_GEBVs
  pred_outputs_rrBLUP$avg_GEBVs

  #########

  ######### GBLUP
  #Note: should wrap the GBLUP into a function too
  # Define variables
  cycles <- input$pred_cv
  folds <- 5
  total_population <- ncol(pred_inputs$geno_input)
  #train_size <- floor(percentage / 100 * total_population)
  cores <- as.numeric(cores)
  #Cross validation number for progress bar (not involved in the calculations, just shiny visuals)
  pb_value = 10

  if (advanced_options$pred_matrix == "Gmatrix") {
    #Convert normalized genotypes to relationship matrix
    #By default, it removes SNPs with more than 50% missing data and imputes using the mean
    Geno.mat <- A.mat(t(pred_inputs$geno_input))

  }else if (advanced_options$pred_matrix == "Amatrix") {

    #Import pedigree file, where pedigree data name (3-column way format). Unknown value should be equal 0
    ped <- read.csv(advanced_options$ped_file$datapath, header = TRUE, check.names = FALSE, colClasses = "factor")
    colnames(ped) <- c("Ind", "Sire", "Dam")
    #Convert NAs to 0
    ped[is.na(ped)] <- 0
    #Ensure Sire and Dam are also listed as individuals
    missing_parents <- unique(c(ped$Sire, ped$Dam))
    # Filter out parents already listed as individuals and non-zero values
    missing_parents <- missing_parents[!missing_parents %in% ped$Ind & missing_parents != 0]
    # Create new rows for missing parents and setting their parents to 0 (unknown)
    new_rows <- data.frame(Ind = missing_parents, Sire = 0, Dam = 0)
    # Combine the original dataframe with the new rows and remove duplicates
    ped_extended <- unique(rbind(ped, new_rows))

    #Converting to Amatrix
    #Using the default additive relationship options (Amatrix only works for even numbered ploidy)
    Geno.mat <- Amatrix(data = ped_extended, ploidy = ploidy)

    #Filter and order the ped file based on the phenotype file (make sure this is valid to subset after generating)
    pheno_ids <- as.character(rownames(pred_inputs$pheno_input))
    valid_ids <- intersect(pheno_ids, rownames(Geno.mat))
    pred_inputs$pheno_input <- pred_inputs$pheno_input[valid_ids, ]
    Geno.mat <- Geno.mat[valid_ids, valid_ids]

    #Update variable
    total_population <- ncol(Geno.mat)
    print("check15")
  }else if (advanced_options$pred_matrix == "Hmatrix") {
    print("check16")
    #Import pedigree file, where pedigree data name (3-column way format). Unknown value should be equal 0
    ped <- read.csv(advanced_options$ped_file$datapath, header = TRUE, check.names = FALSE, colClasses = "factor")
    colnames(ped) <- c("Ind", "Sire", "Dam")
    #Convert NAs to 0
    ped[is.na(ped)] <- 0
    #Ensure Sire and Dam are also listed as individuals
    missing_parents <- unique(c(ped$Sire, ped$Dam))
    # Filter out parents already listed as individuals and non-zero values
    missing_parents <- missing_parents[!missing_parents %in% ped$Ind & missing_parents != 0]
    # Create new rows for missing parents and setting their parents to 0 (unknown)
    new_rows <- data.frame(Ind = missing_parents, Sire = 0, Dam = 0)
    # Combine the original dataframe with the new rows and remove duplicates
    ped_extended <- unique(rbind(ped, new_rows))

    #Converting to Amatrix
    #Using the default additive relationship options (Amatrix only works for even numbered ploidy)
    Ped.mat <- Amatrix(data = ped_extended, ploidy = ploidy)

    #Filter and order the ped file based on the phenotype file (make sure this is valid to subset after generating)
    pheno_ids <- as.character(rownames(pred_inputs$pheno_input))
    valid_ids <- intersect(pheno_ids, rownames(Ped.mat))
    pred_inputs$pheno_input <- pred_inputs$pheno_input[valid_ids, ]
    Ped.mat <- Ped.mat[valid_ids, valid_ids]

    #Update variable
    total_population <- ncol(Ped.mat)

    #Using Gmatrix to get the Gmatrix instead of A.mat for consistency
    #Should I be using the raw dosage values or is it okay to use the scaled genotype data that is used for A.mat()?
    G.mat <- Gmatrix(t(pred_inputs$geno_input[ ,valid_ids]), method = "VanRaden", ploidy = as.numeric(ploidy), missingValue = "NA")
    G.mat <- round(G.mat,3) #to be easy to invert

    #Computing H matrix (Martini) - Using the name Geno.mat for consistency
    Geno.mat <- Hmatrix(A=Ped.mat, G=G.mat, method="Martini",
                        ploidy= ploidy,
                        maf=0.05)
    #Clean memory
    rm(G.mat)
    rm(Ped.mat)
    rm(ped_filtered)
  }
  # Establish accuracy results matrix
  results <- matrix(nrow = cycles*Folds, ncol = length(traits) + 2)
  colnames(results) <- c(paste0(traits), "Iter", "Fold")  # Set the column names to be the traits

  # Initialize a list to store GEBVs for all traits and cycles
  GEBVs <- list()

  #Establish heritability_scores_df () Maybe get h2 values
  # Establish results matrix
  heritability_scores <- matrix(nrow = cycles*Folds, ncol = length(traits) + 2)
  colnames(heritability_scores) <- c(paste0(traits,"_h2"), "Iter", "Fold")  # Set the column names to be the traits


  # For loop
  for (r in 1:cycles) {
    set.seed(r)
    fold_ids <- sample(rep(1:Folds, length.out = total_population))
    fold_df <- data.frame(Sample = row.names(Geno.mat), FoldID = fold_ids) #Randomly assign each sample to a fold
    fold_results <- matrix(nrow = Folds, ncol = length(traits))
    colnames(fold_results) <- traits

    #Initialize GEBV object for each cycle
    GEBVs_cycle <-list()

    #Status
    updateProgressBar(session = session, id = "pb_prediction", value = as.numeric(pb_value), title = paste0("Performing iteration:", r, "of", cycles))

    for (fold in 1:Folds) {

      #Status bar length
      pb_value = pb_value + (70 / as.numeric(cycles*Folds))

      #Subset training and testing samples
      train <- fold_df %>%
        dplyr::filter(FoldID != fold) %>%
        pull(Sample)
      test <- setdiff(row.names(Geno.mat),train)

      Fixed_train = NULL

      # Initialize a matrix to store GEBVs for this fold
      GEBVs_fold <- matrix(nrow = length(test), ncol = length(traits)+3)
      colnames(GEBVs_fold) <- c(traits,"Sample","Iter","Fold")
      rownames(GEBVs_fold) <- paste("Iter", r,"Fold",fold,"Ind", test, sep="_")

      #Evaluate each trait using the same train and testing samples for each
      for (trait_idx in 1:length(traits)) {
        #Mask phenotypes in testing group
        Pheno_test <- pred_inputs$pheno_input
        Pheno_test[test, traits[trait_idx]] <- NA
        #Kin.blup
        traitpred <- kin.blup(data = Pheno_test, geno = names(pred_inputs$pheno_input)[1], pheno = traits[trait_idx], fixed = fixed_cat, covariate = fixed_cov, K=Geno.mat)
        #Cor between test values and predicted breeding values
        results[(((r-1)*5)+fold), trait_idx] <- cor(pred_inputs$pheno_input[test, traits[trait_idx]], traitpred$g[test], use = "complete.obs")
        results[(((r-1)*5)+fold), (length(traits)+1)] <- r
        results[(((r-1)*5)+fold), (length(traits)+2)] <- fold

        # Extract GEBVs
        GEBVs_fold[, trait_idx] <- traitpred$g[test] #Confirm it is accuract to calculate the GEBVs for testing group from the trained model


        # Calculate heritability (these are wrong)
        Vu <- traitpred$Vg
        Ve <- traitpred$Ve
        heritability_scores[(((r-1)*5)+fold), trait_idx] <- Vu / (Vu + Ve)

      }
      #Add iter and fold information for each trait/result
      heritability_scores[(((r-1)*5)+fold), (length(traits)+1)] <- r
      heritability_scores[(((r-1)*5)+fold), (length(traits)+2)] <- fold

      #Add sample, iteration, and fold information to GEBVs_fold
      GEBVs_fold[,"Iter"] = r
      GEBVs_fold[,"Fold"] = fold
      GEBVs_fold[,"Sample"] <- test

      # Store GEBVs for this fold
      GEBVs_cycle[[fold]] <- GEBVs_fold

    }

    # Store GEBVs for this cycle
    GEBVs[[r]] <- do.call(rbind, GEBVs_cycle)

  }

  # Combine all GEBVs into a single DataFrame
  GEBVs_df <- as.data.frame(do.call(rbind, GEBVs))

  results <- as.data.frame(results)
  heritability_scores <- as.data.frame(heritability_scores)

  # Combine results and heritability_scores using cbind
  combined_results <- cbind(results, heritability_scores)

  #Save to reactive value
  pred_outputs_gBLUP <- pred_outputs
  pred_outputs_gBLUP$corr_output <- results
  pred_outputs_gBLUP$all_GEBVs <- results$GEBVs_df

  # Convert trait columns to numeric
  GEBVs <- GEBVs_df %>%
    mutate(across(all_of(traits), ~ as.numeric(.x)))

  # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
  average_gebvs_df <- GEBVs %>%
    group_by(Sample) %>%
    summarize(across(all_of(traits), mean, na.rm = TRUE))

  pred_outputs_gBLUP$avg_GEBVs <- average_gebvs_df

  columns <- setdiff(colnames(results), c("Iter","Fold"))
  average_accuracy_df <- results %>%
    group_by(Iter) %>%
    summarize(across(all_of(columns), mean, na.rm = TRUE))


  pred_outputs_gBLUP$comb_output <- average_accuracy_df

  df <- pred_outputs_gBLUP$corr_output
  df <- df %>% dplyr::select(-Fold, -Iter)

  #Probably want to add the ability for the user to select which trait(s) to display here

  #Convert to long format for ggplot
  df_long <- pivot_longer(
    df,
    cols = colnames(df),  # Exclude the Cycle column from transformation
    names_to = "Trait",  # New column for trait names
    values_to = "Correlation"  # New column for correlation values
  )

  #This can be adapted if we start comparing more than one GP model
  #Also consider a violin plot to show each cor value
  #plot <- ggplot(df_long, aes(x = factor(Trait), y = Correlation, fill = "red"), fill = "red") +
  plot <- ggplot(df_long, aes(x = "rrBLUP", y = Correlation, fill = "red"), fill = "red") +
    #geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 0.7, outlier.size = 0.2) +
    geom_boxplot() +
    facet_wrap(~ Trait, nrow = 1) +  # Facet by trait, allowing different y-scales
    labs(title = "Predictive Ability by Trait",
         x = " ",
         y = "Pearson Correlation") +
    #theme_minimal() +                      # Using a minimal theme
    theme(legend.position = "none",
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2),
          strip.text.x = element_text(face = "bold"),
          axis.text.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank())

  plot_violin <- ggplot(df_long, aes(x = "rrBLUP", y = Correlation, fill = "red")) +
    geom_violin(trim = TRUE) +  # Add violin plot
    geom_point(position = position_jitter(width = 0.1), color = "black", size = 1.5) +  # Add jittered points
    facet_wrap(~ Trait, nrow = 1) +  # Facet by trait, allowing different y-scales
    labs(title = "Predictive Ability by Trait",
         x = " ",  # x-label is blank because it's not relevant per facet
         y = "Pearson Correlation") +
    theme(legend.position = "none",
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2),
          strip.text.x = element_text(face = "bold"),
          axis.text.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank())


  #Output the genomic prediction correlation box plots
  plots()[[1]]  + scale_fill_manual(values = pred_outputs_gBLUP$colors)

  #Output the genomic prediction correlation box plots
  plots()[[2]] + scale_fill_manual(values = pred_outputs_gBLUP$colors)

  #Output the prediction tables
  pred_outputs_gBLUP$comb_output

  all_GEBVs()

  pred_outputs_gBLUP$comb_output

  comb_output()

  pred_outputs_gBLUP$avg_GEBVs

})




