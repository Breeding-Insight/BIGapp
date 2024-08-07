#' GS UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import utils
#'
#'
#' @importFrom bs4Dash valueBox
#' @importFrom shiny NS tagList
#' @importFrom shinyWidgets virtualSelectInput progressBar
#'
#'
mod_GS_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 3,
             box(title="Inputs", width = 12, collapsible = TRUE, collapsed = FALSE, status = "info", solidHeader = TRUE,
                 fileInput(ns("pred_known_file"), "Choose Training Genotypes File", accept = c(".csv",".vcf",".gz")),
                 fileInput(ns("pred_trait_file"), "Choose Passport File", accept = ".csv"),
                 fileInput(ns("pred_est_file"), "Choose Prediction Genotypes File", accept = c(".csv",".vcf",".gz")),
                 numericInput(ns("pred_est_ploidy"), "Species Ploidy", min = 1, value = NULL),
                 virtualSelectInput(
                   inputId = ns("pred_trait_info2"),
                   label = "Select Trait (eg, Color):",
                   choices = NULL,
                   showValueAsTags = TRUE,
                   search = TRUE,
                   multiple = TRUE
                 ),
                 virtualSelectInput(
                   inputId = ns("pred_fixed_info2"),
                   label = "Select Fixed Effects (optional) (not validated):",
                   choices = NULL,
                   showValueAsTags = TRUE,
                   search = TRUE,
                   multiple = TRUE
                 ),
                 actionButton(ns("prediction_est_start"), "Run Analysis"),
                 div(style="display:inline-block; float:right",dropdownButton(
                   tags$h3("GP Parameters"),
                   "GP uses the rrBLUP package: It can impute missing data, adapt to different ploidy, perform 5-fold cross validations with different number of iterations, run multiple traits, and accept multiple fixed effects.",
                   circle = FALSE,
                   status = "warning",
                   icon = icon("info"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see info!")
                 ))

             )
      ),

      column(width = 6,
             box(title = "Results", status = "info", solidHeader = FALSE, width = 12, height = 600,
                 bs4Dash::tabsetPanel(
                   tabPanel("Predicted Trait Table", DTOutput(ns("pred_trait_table")), style = "overflow-y: auto; height: 500px"),
                   tabPanel("GEBVs Table", DTOutput(ns("pred_gebvs_table2")),style = "overflow-y: auto; height: 500px")

                 )
             )
      ),

      column(width = 3,
             valueBoxOutput("shared_snps", width = NULL),
             box(title = "Status", width = 12, collapsible = TRUE, status = "info",
                 progressBar(id = "pb_gp", value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
             ),
             box(title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
                 div(style="display:inline-block; float:left",dropdownButton(
                   tags$h3("Save Files"),
                   fluidRow(
                     downloadButton(ns("download_pred_results_file"), "Save Files")),
                   circle = FALSE,
                   status = "danger",
                   icon = icon("floppy-disk"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see inputs!")
                 ))
             )

      )

    )

  )
}

#' GS Server Functions
#'
#' @importFrom vcfR read.vcfR extract.gt
#' @importFrom rrBLUP mixed.solve A.mat
#' @importFrom stats cor
#' @importFrom shinyalert shinyalert
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @importFrom DT renderDT
#' @noRd
mod_GS_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    ###Genomic Prediction
    #This tab involved 3 observeEvents
    #1) to get the traits listed in the phenotype file
    #2) to input and validate the input files
    #3) to perform the genomic prediction


    #1) Get traits
    observeEvent(input$pred_trait_file, {
      info_df2 <- read.csv(input$pred_trait_file$datapath, header = TRUE, check.names = FALSE, nrow = 0)
      trait_var2 <- colnames(info_df2)
      trait_var2 <- trait_var2[2:length(trait_var2)]
      #updateSelectInput(session, "pred_trait_info", choices = c("All", trait_var))
      updateVirtualSelect("pred_fixed_info2", choices = trait_var2, session = session)
      updateVirtualSelect("pred_trait_info2", choices = trait_var2, session = session)

      #output$passport_table <- renderDT({info_df}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 4)
      #)
    })

    #2) Error check for prediction and save input files
    continue_prediction2 <- reactiveVal(NULL)
    pred_inputs2 <- reactiveValues(
      pheno_input = NULL,
      train_geno_input = NULL,
      est_geno_input = NULL,
      shared_snps = NULL,
      pred_genos = NULL,
      pred_geno_pheno = NULL
    )

    pred_outputs2 <- reactiveValues(
      corr_output = NULL,
      box_plot = NULL,
      violin_plot = NULL,
      comb_output = NULL,
      avg_GEBVs = NULL,
      all_GEBVs = NULL,
      colors = NULL,
      trait_output = NULL
    )

    #Reactive boxes
    output$shared_snps <- renderValueBox({
      valueBox(
        value = pred_inputs2$shared_snps,
        subtitle = "Common SNPs in Genotype files",
        icon = icon("dna"),
        color = "info"
      )
    })

    observeEvent(input$prediction_est_start, {
      #req(pred_inputs$pheno_input, pred_inputs$geno_input)

      #Status
      updateProgressBar(session = session, id = "pb_gp", value = 5, title = "Checking input files")

      #Variables
      ploidy <- as.numeric(input$pred_est_ploidy)
      train_geno_path <- input$pred_known_file$datapath
      est_geno_path <- input$pred_est_file$datapath
      pheno2 <- read.csv(input$pred_trait_file$datapath, header = TRUE, check.names = FALSE)
      row.names(pheno2) <- pheno2[,1]
      traits <- input$pred_trait_info2
      #CVs <- as.numeric(input$pred_cv)
      #train_perc <- as.numeric(input$pred_folds)


      #Make sure at least one trait was input
      if (length(traits) == 0) {

        # If condition is met, show notification toast
        shinyalert(
          title = "Oops",
          text = "No traits were selected",
          size = "xs",
          closeOnEsc = TRUE,
          closeOnClickOutside = FALSE,
          html = TRUE,
          type = "info",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          imageUrl = "",
          animation = TRUE,
        )


        # Stop the observeEvent gracefully
        return()

      }


      #Getting genotype matrix

      #Geno file path
      file_path <- train_geno_path

      #Geno.file conversion if needed
      if (grepl("\\.csv$", file_path)) {
        train_geno <- read.csv(train_geno_path, header = TRUE, row.names = 1, check.names = FALSE)
        est_geno <- read.csv(est_geno_path, header = TRUE, row.names = 1, check.names = FALSE)

        #Save number of SNPs
        #pred_inputs$pred_snps <- nrow(geno)

      } else if (grepl("\\.vcf$", file_path) || grepl("\\.gz$", file_path)) {

        #Function to convert GT to dosage calls (add to BIGr)
        convert_to_dosage <- function(gt) {
          # Split the genotype string
          alleles <- strsplit(gt, "[|/]")
          # Sum the alleles, treating NA values appropriately
          sapply(alleles, function(x) {
            if (any(is.na(x))) {
              return(NA)
            } else {
              return(sum(as.numeric(x), na.rm = TRUE))
            }
          })
        }

        #Convert VCF file if submitted
        train_vcf <- vcfR::read.vcfR(train_geno_path)
        est_vcf <- vcfR::read.vcfR(est_geno_path)

        #Get number of SNPs
        #pred_inputs$pred_snps <- nrow(vcf)

        #Extract GT
        train_geno <- extract.gt(train_vcf, element = "GT")
        train_geno <- apply(train_geno, 2, convert_to_dosage)
        est_geno <- extract.gt(est_vcf, element = "GT")
        est_geno <- apply(est_geno, 2, convert_to_dosage)
        class(train_geno) <- "numeric"
        class(est_geno) <- "numeric"
        rm(train_vcf)
        rm(est_vcf)

      } else {

        # If condition is met, show notification toast
        shinyalert(
          title = "Oops",
          text = "No valid genotype file detected",
          size = "xs",
          closeOnEsc = TRUE,
          closeOnClickOutside = FALSE,
          html = TRUE,
          type = "info",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          imageUrl = "",
          animation = TRUE,
        )

        #Stop the analysis
        return()
      }

      #Save number of samples in file
      #pred_inputs$pred_genos <- ncol(geno)

      #Check that the ploidy entered is correct
      if (ploidy != max(train_geno, na.rm = TRUE)) {
        # If condition is met, show notification toast
        shinyalert(
          title = "Ploidy Mismatch",
          text = paste0("The maximum value in the genotype file (",max(train_geno, na.rm = TRUE),") does not equal the ploidy entered"),
          size = "xs",
          closeOnEsc = FALSE,
          closeOnClickOutside = FALSE,
          html = TRUE,
          type = "warning",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          #closeOnConfirm = TRUE,
          #closeOnCancel = TRUE,
          imageUrl = "",
          animation = TRUE
        )


        # Stop the observeEvent gracefully
        #return()
      }


      # Function to convert genotype matrix according to ploidy
      convert_genotype <- function(genotype_matrix, ploidy) {
        normalized_matrix <- 2 * (genotype_matrix / ploidy) - 1
        return(normalized_matrix)
      }

      #tranforming genotypes
      train_geno_adj_init <- convert_genotype(train_geno, as.numeric(ploidy))
      est_geno_adj_init <- convert_genotype(est_geno, as.numeric(ploidy))

      #Make sure the trait file and genotype file are in the same order
      # Column names for geno (assuming these are the individual IDs)
      colnames_geno <- colnames(train_geno)
      # Assuming the first column in Pheno contains the matching IDs
      ids_pheno <- pheno2[, 1]
      # Find common identifiers
      common_ids <- intersect(colnames_geno, ids_pheno)
      #Get number of id
      pred_inputs2$pred_geno_pheno <- length(common_ids)

      #Throw an error if there are less matching samples in the phenotype file than the genotype file
      if (length(common_ids) == 0) {

        # If condition is met, show notification toast
        shinyalert(
          title = "Oops",
          text = "All samples were missing from the phenotype file",
          size = "xs",
          closeOnEsc = TRUE,
          closeOnClickOutside = FALSE,
          html = TRUE,
          type = "info",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          imageUrl = "",
          animation = TRUE,
        )


        # Stop the observeEvent gracefully
        return()

      } else if (length(common_ids) < length(colnames_geno)) {
        # If condition is met, show notification toast
        shinyalert(
          title = "Data Mismatch",
          text = paste0((length(colnames_geno)-length(common_ids))," samples were removed for not having trait information"),
          size = "xs",
          closeOnEsc = FALSE,
          closeOnClickOutside = FALSE,
          html = TRUE,
          type = "warning",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          #closeOnConfirm = TRUE,
          #closeOnCancel = TRUE,
          imageUrl = "",
          animation = TRUE
        )


        # Stop the observeEvent gracefully
        #return()
      }




      #Final check before performing analyses
      shinyalert(
        title = "Ready?",
        text = "Inputs have been checked",
        size = "xs",
        closeOnEsc = FALSE,
        closeOnClickOutside = FALSE,
        html = TRUE,
        type = "info",
        showConfirmButton = TRUE,
        confirmButtonText = "Proceed",
        confirmButtonCol = "#004192",
        showCancelButton = TRUE,
        #closeOnConfirm = TRUE,
        #closeOnCancel = TRUE,
        imageUrl = "",
        animation = TRUE,
        callbackR = function(value) {
          if (isTRUE(value)) {
            # Proceed with adjusted data
            continue_prediction2(TRUE)
          } else {
            # Stop or change the process
            continue_prediction2(FALSE)
          }
        }
      )


      # Subset and reorder geno and pheno to ensure they only contain and are ordered by common IDs
      train_geno_adj <- train_geno_adj_init[, common_ids]  # Assuming that the columns can be directly indexed by IDs
      pheno2 <- pheno2[match(common_ids, ids_pheno), ]

      #Save to reactive values
      pred_inputs2$pheno_input <- pheno2
      #pred_inputs$geno_adj_input <- geno_adj

      #Match training and testing genotype file SNPs
      common_markers <- intersect(rownames(train_geno_adj), rownames(est_geno_adj_init))
      train_geno_adj <- train_geno_adj[common_markers, ]
      est_geno_adj_init <- est_geno_adj_init[common_markers, ]

      #Save to reactive values
      pred_inputs2$shared_snps <- length(common_markers)
      pred_inputs2$train_geno_input <- train_geno_adj
      pred_inputs2$est_geno_input <- est_geno_adj_init

    })

    #3) Analysis only proceeds once continue_prediction is converted to TRUE
    observe({

      req(continue_prediction2(),pred_inputs2$pheno_input, pred_inputs2$train_geno_input)

      # Stop analysis if cancel was selected
      if (isFALSE(continue_prediction2())) {
        return()
      }

      #Variables
      ploidy <- as.numeric(input$pred_est_ploidy)
      train_geno_adj <- pred_inputs2$train_geno_input
      est_geno_adj <- pred_inputs2$est_geno_input
      pheno <- pred_inputs2$pheno_input
      traits <- input$pred_trait_info2
      #CVs <- as.numeric(input$pred_cv)
      #train_perc <- as.numeric(input$pred_folds)
      fixed_traits <- input$pred_fixed_info2
      cores <- input$pred_cores

      ##Need to add ability for the use of parallelism for the for cross-validation
      ##Example at this tutorial also: https://www.youtube.com/watch?v=ARWjdQU6ays

      # Function to perform genomic prediction
      ##Make sure this is correct (I think I need to be generating a relationship matrix A.mat() to account for missing data, but I am not sure how that works with this)
      genomic_prediction2 <- function(train_geno,est_geno, Pheno, traits, fixed_effects = NULL, cores = 1) {

        # Define variables
        traits <- traits
        #cycles <- as.numeric(Iters)
        #Folds <- as.numeric(Fold)
        total_population <- ncol(train_geno)
        #train_size <- floor(percentage / 100 * total_population)
        fixed_traits <- fixed_effects
        cores <- as.numeric(cores)

        # Initialize a list to store GEBVs for all traits and cycles
        GEBVs <- list()

        #Cross validation number for progress bar (not involved in the calculations, just shiny visuals)
        pb_value = 10

        #Remove the fixed traits from the Pheno file
        if (length(fixed_traits) == 0) {
          Pheno <- Pheno
        } else {
          #Subset fixed traits
          Fixed <- subset(Pheno, select = fixed_traits)

          #Pheno <- subset(Pheno, select = -fixed_traits)
          convert_all_to_factor_if_not_numeric <- function(df) {
            for (col in names(df)) {
              if (!is.numeric(df[[col]]) && !is.integer(df[[col]])) {
                df[[col]] <- as.factor(df[[col]])
              }
            }
            return(df)
          }
          # Convert all columns to factor if they are not numeric or integer
          Fixed <- convert_all_to_factor_if_not_numeric(Fixed)

          #Fixed <- as.data.frame(lapply(Fixed, as.factor)) #convert to factor
          row.names(Fixed) <- row.names(Pheno)

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
        impute = (A.mat(t(train_geno), max.missing=1,impute.method="mean",return.imputed=TRUE))
        train_geno <- impute$imputed
        impute = (A.mat(t(est_geno), max.missing=1,impute.method="mean",return.imputed=TRUE))
        est_geno <- impute$imputed

        #Match training and testing genotype file SNPs
        common_markers <- intersect(colnames(train_geno), colnames(est_geno))
        train_geno <- train_geno[ ,common_markers]
        est_geno <- est_geno[ ,common_markers]

        #Calculate predicted traits and GEBVs
        #fold_ids <- sample(rep(1:Folds, length.out = total_population))
        #fold_df <- data.frame(Sample = row.names(geno), FoldID = fold_ids) #Randomly assign each sample to a fold
        #fold_results <- matrix(nrow = Folds, ncol = length(traits))
        #colnames(fold_results) <- traits

        #Status
        updateProgressBar(session = session, id = "pb_gp", value = 50, title = "Estimating Predicted Values")

        train <- row.names(train_geno)

        #Subset datasets
        #if (length(fixed_traits) == 0) {
        #  Fixed_train = NULL
        #} else{
        #  Fixed_train <- data.frame(Fixed[train, ])
        #  Fixed_train <- as.matrix(Fixed_train)
        #  row.names(Fixed_train) <- train
        #colnames(Fixed_train) <- colnames(Fixed)
        Fixed_train = NULL

        #Fixed (testing)
        #  Fixed_test<- data.frame(Fixed[test, ])
        #  Fixed_test <- as.matrix(Fixed_test)
        #  row.names(Fixed_test) <- test
        #colnames(Fixed_test) <- colnames(Fixed)

        Pheno_train <- Pheno[train, ] # Subset the phenotype df to only retain the relevant samples from the training set
        m_train <- train_geno
        #Pheno_test <- Pheno[test, ]
        #Fixed_test <- Fixed[test, ] #Where would the Fixed_test be used?
        m_valid <- est_geno

        print(dim(m_train))
        print(dim(m_valid))

        # Initialize a matrix to store GEBVs for this fold
        GEBVs_fold <- matrix(nrow = nrow(est_geno), ncol = length(traits))
        colnames(GEBVs_fold) <- c(traits)
        rownames(GEBVs_fold) <- row.names(est_geno)

        Pred_results <- matrix(nrow = nrow(est_geno), ncol = length(traits))
        colnames(Pred_results) <- c(traits)
        rownames(Pred_results) <- row.names(est_geno)

        #Evaluate each trait using the same train and testing samples for each
        for (trait_idx in 1:length(traits)) {
          trait <- Pheno_train[, traits[trait_idx]] # Get the trait of interest
          trait_answer <- mixed.solve(y= trait, Z = m_train, K = NULL, X = Fixed_train, SE = FALSE, return.Hinv = FALSE)
          TRT <- trait_answer$u
          e <- as.matrix(TRT)
          pred_trait_test <- m_valid %*% e
          pred_trait <- pred_trait_test[, 1] + c(trait_answer$beta) # Make sure this still works when using multiple traits
          Pred_results[, trait_idx] <- pred_trait #save to dataframe

          # Extract GEBVs
          # Check if Fixed_train is not NULL and include beta if it is
          if (!is.null(Fixed_train) && !is.null(trait_answer$beta)) {
            # Calculate GEBVs including fixed effects
            #GEBVs_fold[, trait_idx] <- m_train %*% trait_answer$u + Fixed_train %*% matrix(trait_answer$beta, nrow = length(trait_answer$beta), ncol = 1)
            #GEBVs_fold[, trait_idx] <- m_valid %*% trait_answer$u + Fixed_test %*% matrix(trait_answer$beta, nrow = length(trait_answer$beta), ncol = 1)
            GEBVs_fold[, trait_idx] <- m_valid %*% trait_answer$u + Fixed_test %*% trait_answer$beta
          } else {
            # Calculate GEBVs without fixed effects
            #GEBVs_fold[, trait_idx] <- m_train %*% trait_answer$u
            GEBVs_fold[, trait_idx] <- m_valid %*% trait_answer$u #Confirm it is accuract to calculate the GEBVs for testing group from the trained model
          }

          # Calculate heritability for the current trait
          #Vu <- trait_answer$Vu
          #Ve <- trait_answer$Ve
          #heritability_scores[(((r-1)*5)+fold), trait_idx] <- Vu / (Vu + Ve)

        }
        #Add iter and fold information for each trait/result
        #heritability_scores[(((r-1)*5)+fold), (length(traits)+1)] <- r
        #heritability_scores[(((r-1)*5)+fold), (length(traits)+2)] <- fold

        #Add sample, iteration, and fold information to GEBVs_fold
        #GEBVs_fold[,"Iter"] = r
        #GEBVs_fold[,"Fold"] = fold
        #GEBVs_fold[,"Sample"] <- row.names(est_geno)

        # Store GEBVs for this fold
        GEBVs_df <- data.frame(GEBVs_fold)

        Pred_results <- data.frame(Pred_results)


        # Store GEBVs for this cycle
        #GEBVs[[r]] <- do.call(rbind, GEBVs_cycle)


        # Combine all GEBVs into a single DataFrame
        #GEBVs_df <- as.data.frame(do.call(rbind, GEBVs))

        #results <- as.data.frame(results)
        #heritability_scores <- as.data.frame(heritability_scores)

        # Combine results and heritability_scores using cbind
        #combined_results <- cbind(results, heritability_scores)

        return(list(GEBVs = GEBVs_df, Predictions = Pred_results))
      }

      # Example call to the function
      #This is slow when using 3k markers and 1.2k samples...will need to parallelize if using this script...
      results <- genomic_prediction2(train_geno_adj, est_geno_adj, pheno, traits = traits, fixed_effects = fixed_traits, cores = cores)

      #With fixed effects (need to inforporate the ability for fixed effects into the prediction?)
      #results <- genomic_prediction(geno_matrix, phenotype_df, c("height", "weight"), "~ age + sex")

      #Save to reactive value
      pred_outputs2$trait_output <- results$Predictions
      pred_outputs2$all_GEBVs <- results$GEBVs
      #TESTING!!!
      #write.csv(results$GEBVs, "GEBVs_test.csv")

      # Convert trait columns to numeric
      results$GEBVs <- results$GEBVs %>%
        mutate(across(all_of(traits), ~ as.numeric(.x)))


      #Get average accuracy and h2 for each iter accross the 5 folds

      #columns <- setdiff(colnames(results$CombinedResults), c("Iter","Fold"))
      #average_accuracy_df <- results$CombinedResults %>%
      #  group_by(Iter) %>%
      #  summarize(across(all_of(columns), mean, na.rm = TRUE))


      #Status
      updateProgressBar(session = session, id = "pb_gp", value = 90, title = "Generating Results")

      ##Figures and Tables

      #Status
      updateProgressBar(session = session, id = "pb_gp", value = 100, title = "Finished!")

      #End the event
      continue_prediction2(NULL)
    })

    #Output the prediction tables
    observe({
      #GEBVs from all iterations/folds
      req(pred_outputs2$all_GEBVs)

      output$pred_gebvs_table2 <- renderDT({pred_outputs2$all_GEBVs}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

    })

    observe({
      #GEBVs from all iterations/folds
      req(pred_outputs2$trait_output)

      output$pred_trait_table <- renderDT({pred_outputs2$trait_output}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

    })
  })
}

## To be copied in the UI
# mod_GS_ui("GS_1")

## To be copied in the server
# mod_GS_server("GS_1")
