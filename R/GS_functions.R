#' Function to perform genomic prediction
#'
#' @param geno ToDo
#' @param pheno ToDo
#' @param traits ToDo
#' @param fixed_effects ToDo
#' @param fixed_cat categorial fixed effects
#' @param folds ToDo
#' @param iters ToDo
#' @param cores ToDo
#' @param session ToDo
#'
rrBLUP_genomic_prediction <- function(geno, pheno, traits, fixed_effects = NULL, fixed_cat = NULL,folds = 5, iters = 5, cores = 1, session) {

  # Define variables
  cycles <- as.numeric(iters)
  folds <- as.numeric(folds)
  total_population <- ncol(geno)
  #train_size <- floor(percentage / 100 * total_population)
  fixed_traits <- fixed_effects
  cores <- as.numeric(cores)

  # Establish accuracy results matrix
  results <- matrix(nrow = cycles*folds, ncol = length(traits) + 2)
  colnames(results) <- c(paste0(traits), "Iter", "Fold")  # Set the column names to be the traits

  # Initialize a list to store GEBVs for all traits and cycles
  GEBVs <- list()

  #Establish heritability_scores_df () Maybe get h2 values
  # Establish results matrix
  heritability_scores <- matrix(nrow = cycles*folds, ncol = length(traits) + 2)
  colnames(heritability_scores) <- c(paste0(traits,"_h2"), "Iter", "Fold")  # Set the column names to be the traits

  #Cross validation number for progress bar (not involved in the calculations, just shiny visuals)
  pb_value = 10
  #Remove the fixed traits from the pheno file
  if (length(fixed_traits) == 0) {
    pheno <- pheno
  } else {
    #Subset fixed traits
    Fixed <- subset(pheno, select = fixed_traits)

    #pheno <- subset(pheno, select = -fixed_traits)
    convert_categorical_to_factor <- function(df, fixed_cat) {
      for (col in names(df)) {
        if (col %in% fixed_cat) {
          df[[col]] <- as.factor(df[[col]])
        }
      }
      return(df)
    }
    # Convert all columns to factor if they are not numeric or integer
    Fixed <- convert_categorical_to_factor(Fixed, fixed_cat)

    #Fixed <- as.data.frame(lapply(Fixed, as.factor)) #convert to factor
    row.names(Fixed) <- row.names(pheno)

    #Make the matrix
    formula_str <- paste("~", paste(fixed_traits, collapse = " + "))
    formula <- as.formula(formula_str)

    # Create the design matrix using the constructed formula
    Fixed <- model.matrix(formula, data = Fixed)
  }

  #Make kinship matrix of all individuals?
  #Kin_mat <- A.mat(t(geno), n.core = 1) ##Need to explore whether or not to use a Kinship matrix and if it makes a significant improvement to accuracy
  #If wanting to use Kkinship matrix, will then need to see how to implement it here

  #For now, I am just imputing the missing sites using mean, but EM is more accurate, but slower (can use multiple cores).
  impute = (A.mat(t(geno), max.missing=0.5,impute.method="mean",return.imputed=TRUE))
  geno <- impute$imputed

  # For loop
  for (r in 1:cycles) {
    set.seed(r)
    fold_ids <- sample(rep(1:folds, length.out = total_population))
    fold_df <- data.frame(Sample = row.names(geno), FoldID = fold_ids) #Randomly assign each sample to a fold
    fold_results <- matrix(nrow = folds, ncol = length(traits))
    colnames(fold_results) <- traits

    #Initialize GEBV object for each cycle
    GEBVs_cycle <-list()

    #Status
    if(!is.null(session)) updateProgressBar(session = session, id = "pb_prediction", value = as.numeric(pb_value), title = paste0("Performing iteration:", r, "of", cycles))

    for (fold in 1:folds) {

      #Status bar length
      pb_value = pb_value + (70 / as.numeric(cycles*folds))

      train <- fold_df %>%
        dplyr::filter(FoldID != fold) %>%
        pull(Sample)
      test <- setdiff(row.names(geno),train)

      #Subset datasets
      if (length(fixed_traits) == 0) {
        Fixed_train = NULL
      } else{
        Fixed_train <- data.frame(Fixed[train, ])
        Fixed_train <- as.matrix(Fixed_train)
        row.names(Fixed_train) <- train

        #Fixed (testing)
        Fixed_test<- data.frame(Fixed[test, ])
        Fixed_test <- as.matrix(Fixed_test)
        row.names(Fixed_test) <- test

      }

      pheno_train <- pheno[train, ] # Subset the phenotype df to only retain the relevant samples from the training set
      m_train <- geno[train, ]
      pheno_test <- pheno[test, ]
      #Fixed_test <- Fixed[test, ] #Where would the Fixed_test be used?
      m_valid <- geno[test, ]

      # Initialize a matrix to store GEBVs for this fold
      GEBVs_fold <- matrix(nrow = length(test), ncol = length(traits)+3)
      colnames(GEBVs_fold) <- c(traits,"Sample","Iter","Fold")
      rownames(GEBVs_fold) <- paste("Iter", r,"Fold",fold,"Ind", test, sep="_")

      #Evaluate each trait using the same train and testing samples for each
      for (trait_idx in 1:length(traits)) {
        trait <- pheno_train[, traits[trait_idx]] # Get the trait of interest
        trait_answer <- mixed.solve(y= trait, Z = m_train, K = NULL, X = Fixed_train, SE = FALSE, return.Hinv = FALSE)
        TRT <- trait_answer$u
        e <- as.matrix(TRT)
        pred_trait_test <- m_valid %*% e
        pred_trait <- pred_trait_test[, 1] + c(trait_answer$beta) # Make sure this still works when using multiple traits
        trait_test <- pheno_test[, traits[trait_idx]]
        results[(((r-1)*5)+fold), trait_idx] <- cor(pred_trait, trait_test, use = "complete")
        results[(((r-1)*5)+fold), (length(traits)+1)] <- r
        results[(((r-1)*5)+fold), (length(traits)+2)] <- fold

        # Extract GEBVs
        # Check if Fixed_train is not NULL and include beta if it is
        if (!is.null(Fixed_train) && !is.null(trait_answer$beta)) {
          # Calculate GEBVs including fixed effects
          GEBVs_fold[, trait_idx] <- m_valid %*% trait_answer$u + Fixed_test %*% trait_answer$beta
        } else {
          # Calculate GEBVs without fixed effects
          GEBVs_fold[, trait_idx] <- m_valid %*% trait_answer$u #Confirm it is accurate to calculate the GEBVs for testing group from the trained model
        }

        # Calculate heritability for the current trait
        Vu <- trait_answer$Vu
        Ve <- trait_answer$Ve
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

  return(list(GEBVs = GEBVs_df, PredictionAccuracy = results, CombinedResults = combined_results))
}


#' Compute relationship matrix
#'
#' @param type character defining which type
#' @param ploidy numeric indicating species ploidy
#' @param geno_input matrix with individuals in the row and markers in the columns
#' @param ped_file pedigree file
#' @param pheno ToDo
#'
#' @importFrom rrBLUP A.mat
#' @importFrom AGHmatrix Gmatrix Amatrix Hmatrix
#'
#'
get_relationship_mat <- function(geno_input, ped_file, type = c("Gmatrix", "Amatrix", "Hmatrix"), ploidy, pheno){

  if (type == "Gmatrix") {
    #Convert normalized genotypes to relationship matrix
    #By default, it removes SNPs with more than 50% missing data and imputes using the mean
    Geno.mat <- Gmatrix(t(geno_input),
                        method = "VanRaden",
                        ploidy = ploidy,
                        ploidy.correction=TRUE,
                        ratio = FALSE,
                        missingValue = "NA")
    return(Geno.mat)

  }else if (type == "Amatrix") {

    #Converting to Amatrix
    #Using the default additive relationship options (Amatrix only works for even numbered ploidy)
    Geno.mat <- Amatrix(data = ped_file, ploidy = ploidy)

    #Filter and order the ped file based on the phenotype file (make sure this is valid to subset after generating)
    pheno_ids <- as.character(pheno[,1])
    valid_ids <- intersect(pheno_ids, rownames(Geno.mat))
    Geno.mat <- Geno.mat[valid_ids, valid_ids]

    return(Geno.mat)

  }else if (type == "Hmatrix") {

    #Converting to Amatrix
    #Using the default additive relationship options (Amatrix only works for even numbered ploidy)
    Ped.mat <- Amatrix(data = ped_file, ploidy = ploidy)

    #Filter and order the ped file based on the phenotype file (make sure this is valid to subset after generating)
    pheno_ids <- as.character(pheno[,1])
    valid_ids <- intersect(pheno_ids, rownames(Ped.mat))
    Ped.mat <- Ped.mat[valid_ids, valid_ids]

    #Using Gmatrix to get the Gmatrix instead of A.mat for consistency
    #Should I be using the raw dosage values or is it okay to use the scaled genotype data that is used for A.mat()?
    G.mat <- Gmatrix(t(geno_input[ ,valid_ids]),
                     method = "VanRaden",
                     ploidy = as.numeric(ploidy),
                     missingValue = "NA",
                     ploidy.correction = TRUE)
    G.mat <- round(G.mat,3) #to be easy to invert

    #Computing H matrix (Martini) - Using the name Geno.mat for consistency
    Geno.mat <- Hmatrix(A=Ped.mat, G=G.mat,
                        method="Martini",
                        ploidy= ploidy,
                        maf=0.05)
    return(Geno.mat)
  }
}


#' Performes GBLUP
#'
#' @param pheno_dat ToDo
#' @param Geno.mat ToDo
#' @param cycles ToDo
#' @param folds ToDo
#' @param traits ToDo
#' @param cores ToDo
#' @param fixed_cov ToDo
#' @param fixed_cat ToDo
#' @param session ToDo
#'
GBLUP_genomic_prediction <- function(pheno_dat, Geno.mat, cycles, folds, traits, cores, fixed_cov = NULL, fixed_cat = NULL, session = NULL){

  # Establish accuracy results matrix
  results <- matrix(nrow = cycles*folds, ncol = length(traits) + 2)
  colnames(results) <- c(paste0(traits), "Iter", "Fold")  # Set the column names to be the traits
  pb_value <- 10

  # Initialize a list to store GEBVs for all traits and cycles
  GEBVs <- list()

  #Establish heritability_scores_df () Maybe get h2 values
  # Establish results matrix
  heritability_scores <- matrix(nrow = cycles*folds, ncol = length(traits) + 2)
  colnames(heritability_scores) <- c(paste0(traits,"_h2"), "Iter", "Fold")  # Set the column names to be the traits

  # For loop
  for (r in 1:cycles) {
    set.seed(r)
    fold_ids <- sample(rep(1:folds, length.out = nrow(Geno.mat)))
    fold_df <- data.frame(Sample = row.names(Geno.mat), FoldID = fold_ids) #Randomly assign each sample to a fold
    fold_results <- matrix(nrow = folds, ncol = length(traits))
    colnames(fold_results) <- traits

    #Initialize GEBV object for each cycle
    GEBVs_cycle <-list()

    if(!is.null(session)) updateProgressBar(session = session, id = "pb_prediction", value = as.numeric(pb_value), title = paste0("Performing iteration:", r, "of", cycles))

    for (fold in 1:folds) {

      #Status bar length
      pb_value = pb_value + (70 / as.numeric(cycles*folds))

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
        Pheno_test <- pheno_dat
        Pheno_test[test, traits[trait_idx]] <- NA
        #Kin.blup
        traitpred <- kin.blup(data = Pheno_test,
                              geno = colnames(pheno_dat)[1],
                              pheno = traits[trait_idx],
                              fixed = fixed_cat,
                              covariate = fixed_cov,
                              K=Geno.mat,
                              n.core = cores)

        #Cor between test values and predicted breeding values
        results[(((r-1)*5)+fold), trait_idx] <- cor(pheno_dat[test, traits[trait_idx]], traitpred$g[test], use = "complete.obs")
        results[(((r-1)*5)+fold), (length(traits)+1)] <- r
        results[(((r-1)*5)+fold), (length(traits)+2)] <- fold

        # Extract GEBVs
        GEBVs_fold[, trait_idx] <- traitpred$g[test]

        # Calculate heritability (*confirm this calculation* - either way will not report to user)
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

  return(list(GEBVs = GEBVs_df, PredictionAccuracy = results, CombinedResults = combined_results))
}

assign_colors <- function(color){
  if (color == "red"){
    our_color <- "#F8766D"
  } else if (color == "blue") {
    our_color <- "#00BFC4"
  } else if (color == "green") {
    our_color <- "#00BA38"
  } else{
    our_color <- color
  }
  return(our_color)
}

format_geno_matrix <- function(geno, model, pred_matrix = NULL, ploidy){

  if(is.null(pred_matrix)) pred_matrix <- "none_selected"
  if(model == "rrBLUP" & ploidy == 2 & pred_matrix == "Gmatrix") {
    #if(model == "rrBLUP") {
    geno_formated <- 2 * (geno / as.numeric(ploidy)) - 1 # codification -1 0 1
  } else {
    geno_formated <- geno # codification 0 1 2 3 ..
  }

  return(geno_formated)
}

run_predictive_model <- function(geno, pheno, selected_traits, predictive_model, relationship_matrix_type, pedigree,
                                 fixed_effects, categorical_fixed_effects, ploidy, cores, cycles, folds, relationship_matrix = NULL,
                                 session = NULL){

  if(predictive_model == "rrBLUP"){
    results <- rrBLUP_genomic_prediction(geno = geno,
                                         pheno = pheno,
                                         traits = selected_traits,
                                         fixed_effects = fixed_effects,
                                         iters = cycles,
                                         cores = cores,
                                         session = session)
    return(results)
  } else if(predictive_model == "GBLUP"){
    fixed_cov <- if (is.null(fixed_effects) || length(fixed_effects) == length(categorical_fixed_effects)) {
      NULL
    } else {
      setdiff(fixed_effects, categorical_fixed_effects)
    }

    if(is.null(relationship_matrix)){
      Geno.mat <- get_relationship_mat(geno_input = geno,
                                       type = relationship_matrix_type,
                                       ped_file = pedigree,
                                       ploidy = ploidy,
                                       pheno = pheno)
    } else Geno.mat <- relationship_matrix

    results <- GBLUP_genomic_prediction(pheno_dat = pheno,
                                        Geno.mat = Geno.mat,
                                        cycles = cycles,
                                        folds = folds, #?
                                        traits = selected_traits,
                                        cores = cores,
                                        fixed_cov = fixed_cov,
                                        fixed_cat = categorical_fixed_effects,
                                        session = session)
    return(results)
  }
}

# Remove individuals unrelated to the ones that have phenotype info
# Add line with zeros for older generation
remove_unrelated <- function(pedigree, samples_with_trait_info){
  common_ped <- intersect(pedigree$Ind, samples_with_trait_info)
  ped_test <- pedigree[which(pedigree$Ind %in% common_ped),]
  all_gene <- as.character(unique(unlist(ped_test)))
  if(length(which(all_gene == "0")) > 0) all_gene <- all_gene[-which(all_gene == "0")]
  gen <- 1
  while(length(all_gene) != length(ped_test$Ind)){
    gen <- gen + 1
    ped_test <- pedigree[which(pedigree$Ind %in% all_gene),] # add previous generation
    dim(ped_test)
    missing_older_gen <- which(!all_gene %in% pedigree$Ind)
    length(missing_older_gen)
    if(length(missing_older_gen) > 0)  {
      add_previous_gen <- data.frame(Ind = all_gene[missing_older_gen], P1 = 0, P2 = 0) # Add missing previous generation
      ped_test <- rbind(ped_test, add_previous_gen)
    }
    all_gene <- as.character(unique(unlist(ped_test)))
    if(length(which(all_gene == "0")) > 0) all_gene <- all_gene[-which(all_gene == "0")]
    length(all_gene)
  }
  return(list(ped_test, gen))
}
