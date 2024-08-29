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
#' @importFrom shinyWidgets virtualSelectInput
#'
mod_GSAcc_ui <- function(id){
  ns <- NS(id)
  tagList(
    # Add GWAS content here
    fluidRow(
      column(width = 3,
             box(title="Inputs", width = 12, collapsible = TRUE, collapsed = FALSE, status = "info", solidHeader = TRUE,
                 fileInput(ns("pred_file"), "Choose VCF File", accept = c(".csv",".vcf",".gz")),
                 fileInput(ns("trait_file"), "Choose Passport File", accept = ".csv"),
                 numericInput(ns("pred_ploidy"), "Species Ploidy", min = 1, value = NULL),
                 numericInput(ns("pred_cv"), "Iterations", min = 1, max=20, value = 5),
                 virtualSelectInput(
                   inputId = ns("pred_trait_info"),
                   label = "Select Trait (eg, Color):",
                   choices = NULL,
                   showValueAsTags = TRUE,
                   search = TRUE,
                   multiple = TRUE
                 ),
                 virtualSelectInput(
                   inputId = ns("pred_fixed_info"),
                   label = "Select Fixed Effects (optional) (not validated):",
                   choices = NULL,
                   showValueAsTags = TRUE,
                   search = TRUE,
                   multiple = TRUE
                 ),
                 actionButton(ns("prediction_start"), "Run Analysis"),
                 div(style="display:inline-block; float:right",dropdownButton(
                   tags$h3("GP Parameters"),
                   "You can download examples of the expected input input files here: \n",
                   downloadButton(ns('download_vcf'), "Download VCF Example File"),
                   downloadButton(ns('download_pheno'), "Download Passport Example File"),
                   #"GP uses the rrBLUP package: It can impute missing data, adapt to different ploidy, perform 5-fold cross validations with different number of iterations, run multiple traits, and accept multiple fixed effects.",
                   circle = FALSE,
                   status = "warning",
                   icon = icon("info"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see info!")
                 ))

             )
      ),

      column(width = 6,
             box(
               title = "Plots", status = "info", solidHeader = FALSE, width = 12, height = 600,
               bs4Dash::tabsetPanel(
                 tabPanel("Violin Plot", plotOutput(ns("pred_violin_plot"), height = "500px")),
                 tabPanel("Box Plot", plotOutput(ns("pred_box_plot"), height = "500px")),
                 tabPanel("Accuracy Table", DTOutput(ns("pred_acc_table")), style = "overflow-y: auto; height: 500px"),
                 tabPanel("GEBVs Table", DTOutput(ns("pred_gebvs_table")),style = "overflow-y: auto; height: 500px")

               )

             )

      ),

      column(width = 3,
             valueBoxOutput(ns("pred_snps"), width = NULL),
             valueBoxOutput(ns("pred_geno"), width = NULL),
             box(title = "Status", width = 12, collapsible = TRUE, status = "info",
                 progressBar(id = ns("pb_prediction"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
             ),
             box(title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
                 selectInput(ns("pred_color_select"), label = "Color Selection", choices = c("red","orange","yellow","green","blue","violet", "grey", "white")),
                 div(style="display:inline-block; float:left",dropdownButton(
                   tags$h3("Save Image"),
                   selectInput(inputId = ns('pred_figures'), label = 'Figure', choices = c("Violin Plot",
                                                                                           "Box Plot")),
                   selectInput(inputId = ns('pred_image_type'), label = 'File Type', choices = c("jpeg","tiff","png"), selected = "jpeg"),
                   sliderInput(inputId = ns('pred_image_res'), label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                   sliderInput(inputId = ns('pred_image_width'), label = 'Width', value = 10, min = 1, max = 20, step=0.5),
                   sliderInput(inputId = ns('pred_image_height'), label = 'Height', value = 6, min = 1, max = 20, step = 0.5),
                   fluidRow(
                     downloadButton(ns("download_pred_figure"), "Save Image"),
                     downloadButton(ns("download_pred_file"), "Save Files")),
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
#' @noRd
mod_GSAcc_server <- function(input, output, session, parent_session){

    ns <- session$ns
    ####Genomic Prediction Accuracy
    #This tab involved 3 observeEvents
    #1) to get the traits listed in the phenotype file
    #2) to input and validate the input files
    #3) to perform the genomic prediction


    #1) Get traits
    observeEvent(input$trait_file, {
      info_df <- read.csv(input$trait_file$datapath, header = TRUE, check.names = FALSE, nrow = 0)
      trait_var <- colnames(info_df)
      trait_var <- trait_var[2:length(trait_var)]
      updateVirtualSelect("pred_fixed_info", choices = trait_var, session = session)
      updateVirtualSelect("pred_trait_info", choices = trait_var, session = session)

    })

    #2) Error check for prediction and save input files
    continue_prediction <- reactiveVal(NULL)
    pred_inputs <- reactiveValues(
      pheno_input = NULL,
      geno_input = NULL,
      pred_snps = NULL,
      pred_genos = NULL,
      pred_geno_pheno = NULL
    )

    pred_outputs <- reactiveValues(
      corr_output = NULL,
      box_plot = NULL,
      violin_plot = NULL,
      comb_output = NULL,
      avg_GEBVs = NULL,
      all_GEBVs = NULL,
      colors = NULL
    )

    #Reactive boxes
    output$pred_snps <- renderValueBox({
      valueBox(
        value = pred_inputs$pred_snps,
        subtitle = "SNPs in Genotype File",
        icon = icon("dna"),
        color = "info"
      )
    })

    output$pred_geno <- renderValueBox({
      valueBox(
        value = pred_inputs$pred_geno_pheno,
        subtitle = "Samples with Phenotype Information",
        icon = icon("location-dot"),
        color = "info"
      )
    })

    observe({
      # Update colors based on input
      pred_outputs$colors <- switch(input$pred_color_select,
                                    "red" = "#F8766D",
                                    "blue" = "#00BFC4",
                                    "green" = "#00BA38",
                                    input$pred_color_select)
    })

    observeEvent(input$prediction_start, {

      toggleClass(id = "pred_ploidy", class = "borderred", condition = (is.na(input$pred_ploidy) | is.null(input$pred_ploidy)))

      if (is.null(input$pred_file$datapath) | is.null(input$trait_file$datapath)) {
        shinyalert(
          title = "Missing input!",
          text = "Upload VCF and phenotype files",
          size = "s",
          closeOnEsc = TRUE,
          closeOnClickOutside = FALSE,
          html = TRUE,
          type = "error",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          animation = TRUE
        )
      }
      req(input$pred_file$datapath,  input$pred_ploidy, input$trait_file$datapath)

      #Status
      updateProgressBar(session = session, id = "pb_prediction", value = 5, title = "Checking input files")

      #Variables
      ploidy <- as.numeric(input$pred_ploidy)
      geno_path <- input$pred_file$datapath
      pheno <- read.csv(input$trait_file$datapath, header = TRUE, check.names = FALSE)
      row.names(pheno) <- pheno[,1]
      traits <- input$pred_trait_info
      CVs <- as.numeric(input$pred_cv)

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
      file_path <- geno_path

      #Geno.file conversion if needed
      if (grepl("\\.csv$", file_path)) {
        geno <- read.csv(geno_path, header = TRUE, row.names = 1, check.names = FALSE)

        #Save number of SNPs
        pred_inputs$pred_snps <- nrow(geno)

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
        vcf <- vcfR::read.vcfR(file_path)

        #Get number of SNPs
        pred_inputs$pred_snps <- nrow(vcf)

        #Extract GT
        geno <- extract.gt(vcf, element = "GT")
        geno <- apply(geno, 2, convert_to_dosage)
        class(geno) <- "numeric"
        rm(vcf)

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
      pred_inputs$pred_genos <- ncol(geno)

      #Check that the ploidy entered is correct
      if (ploidy != max(geno, na.rm = TRUE)) {
        # If condition is met, show notification toast
        shinyalert(
          title = "Ploidy Mismatch",
          text = paste0("The maximum value in the genotype file (",max(geno, na.rm = TRUE),") does not equal the ploidy entered"),
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
      geno_adj_init <- convert_genotype(geno, as.numeric(ploidy))

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
            continue_prediction(TRUE)
          } else {
            # Stop or change the process
            continue_prediction(FALSE)
          }
        }
      )


      # Subset and reorder geno and pheno to ensure they only contain and are ordered by common IDs
      geno_adj <- geno_adj_init[, common_ids]  # Assuming that the columns can be directly indexed by IDs
      pheno <- pheno[match(common_ids, ids_pheno), ]

      #Save to reactive values
      pred_inputs$pheno_input <- pheno
      #pred_inputs$geno_adj_input <- geno_adj
      pred_inputs$geno_input <- geno_adj

    })

    #3) Analysis only proceeds once continue_prediction is converted to TRUE
    observe({

      req(continue_prediction(),pred_inputs$pheno_input, pred_inputs$geno_input)

      # Stop analysis if cancel was selected
      if (isFALSE(continue_prediction())) {
        return()
      }

      #Variables
      ploidy <- as.numeric(input$pred_ploidy)
      geno_adj <- pred_inputs$geno_input
      pheno <- pred_inputs$pheno_input
      traits <- input$pred_trait_info
      CVs <- as.numeric(input$pred_cv)
      fixed_traits <- input$pred_fixed_info
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

      ##Need to add ability for the use of parallelism for the for cross-validation
      ##Example at this tutorial also: https://www.youtube.com/watch?v=ARWjdQU6ays

      # Function to perform genomic prediction
      ##Make sure this is correct (I think I need to be generating a relationship matrix A.mat() to account for missing data, but I am not sure how that works with this)
      genomic_prediction <- function(geno, Pheno, traits, fixed_effects = NULL, Fold = 5, Iters = 5, cores = 1) {

        # Define variables
        traits <- traits
        cycles <- as.numeric(Iters)
        Folds <- as.numeric(Fold)
        total_population <- ncol(geno)
        #train_size <- floor(percentage / 100 * total_population)
        fixed_traits <- fixed_effects
        cores <- as.numeric(cores)

        # Establish accuracy results matrix
        results <- matrix(nrow = cycles*Folds, ncol = length(traits) + 2)
        colnames(results) <- c(paste0(traits), "Iter", "Fold")  # Set the column names to be the traits

        # Initialize a list to store GEBVs for all traits and cycles
        GEBVs <- list()

        #Establish heritability_scores_df () Maybe get h2 values
        # Establish results matrix
        heritability_scores <- matrix(nrow = cycles*Folds, ncol = length(traits) + 2)
        colnames(heritability_scores) <- c(paste0(traits,"_h2"), "Iter", "Fold")  # Set the column names to be the traits

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
        impute = (A.mat(t(geno), max.missing=0.5,impute.method="mean",return.imputed=TRUE))
        geno <- impute$imputed

        # For loop
        for (r in 1:cycles) {
          set.seed(r)
          fold_ids <- sample(rep(1:Folds, length.out = total_population))
          fold_df <- data.frame(Sample = row.names(geno), FoldID = fold_ids) #Randomly assign each sample to a fold
          fold_results <- matrix(nrow = Folds, ncol = length(traits))
          colnames(fold_results) <- traits

          #Initialize GEBV object for each cycle
          GEBVs_cycle <-list()

          #Status
          updateProgressBar(session = session, id = "pb_prediction", value = as.numeric(pb_value), title = paste0("Performing iteration:", r, "of", cycles))

          for (fold in 1:Folds) {

            #Status bar length
            pb_value = pb_value + (70 / as.numeric(cycles*Folds))

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

            Pheno_train <- Pheno[train, ] # Subset the phenotype df to only retain the relevant samples from the training set
            m_train <- geno[train, ]
            Pheno_test <- Pheno[test, ]
            #Fixed_test <- Fixed[test, ] #Where would the Fixed_test be used?
            m_valid <- geno[test, ]

            # Initialize a matrix to store GEBVs for this fold
            GEBVs_fold <- matrix(nrow = length(test), ncol = length(traits)+3)
            colnames(GEBVs_fold) <- c(traits,"Sample","Iter","Fold")
            rownames(GEBVs_fold) <- paste("Iter", r,"Fold",fold,"Ind", test, sep="_")

            #Evaluate each trait using the same train and testing samples for each
            for (trait_idx in 1:length(traits)) {
              trait <- Pheno_train[, traits[trait_idx]] # Get the trait of interest
              trait_answer <- mixed.solve(y= trait, Z = m_train, K = NULL, X = Fixed_train, SE = FALSE, return.Hinv = FALSE)
              TRT <- trait_answer$u
              e <- as.matrix(TRT)
              pred_trait_test <- m_valid %*% e
              pred_trait <- pred_trait_test[, 1] + c(trait_answer$beta) # Make sure this still works when using multiple traits
              trait_test <- Pheno_test[, traits[trait_idx]]
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

      # Example call to the function
      #This is slow when using 3k markers and 1.2k samples...will need to parallelize if using this script...
      results <- genomic_prediction(geno_adj, pheno, traits = traits, fixed_effects = fixed_traits, Iters = input$pred_cv, cores = cores)

      #With fixed effects (need to inforporate the ability for fixed effects into the prediction?)
      #results <- genomic_prediction(geno_matrix, phenotype_df, c("height", "weight"), "~ age + sex")

      #Save to reactive value
      pred_outputs$corr_output <- results$PredictionAccuracy
      pred_outputs$all_GEBVs <- results$GEBVs

      # Convert trait columns to numeric
      results$GEBVs <- results$GEBVs %>%
        mutate(across(all_of(traits), ~ as.numeric(.x)))

      # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
      average_gebvs_df <- results$GEBVs %>%
        group_by(Sample) %>%
        summarize(across(all_of(traits), mean, na.rm = TRUE))

      pred_outputs$avg_GEBVs <- average_gebvs_df

      #Get average accuracy and h2 for each iter accross the 5 folds

      #columns <- setdiff(colnames(results$CombinedResults), c("Iter","Fold"))
      #average_accuracy_df <- results$CombinedResults %>%
      #  group_by(Iter) %>%
      #  summarize(across(all_of(columns), mean, na.rm = TRUE))

      columns <- setdiff(colnames(results$PredictionAccuracy), c("Iter","Fold"))
      average_accuracy_df <- results$PredictionAccuracy %>%
        group_by(Iter) %>%
        summarize(across(all_of(columns), mean, na.rm = TRUE))


      pred_outputs$comb_output <- average_accuracy_df

      #Status
      updateProgressBar(session = session, id = "pb_prediction", value = 90, title = "Generating Results")

      ##Figures and Tables

      #Status
      updateProgressBar(session = session, id = "pb_prediction", value = 100, title = "Finished!")

      #End the event
      continue_prediction(NULL)
    })

    plots <- reactive({
      validate(
        need(!is.null(pred_outputs$corr_output), "Upload the input files, set the parameters and click 'run analysis' to access results in this session.")
      )

      df <- pred_outputs$corr_output
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
              strip.text.x = element_text(face = "bold"))

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
              strip.text.x = element_text(face = "bold"))

      list(plot, plot_violin)
    })

    #Output the genomic prediction correlation box plots
    output$pred_box_plot <- renderPlot({
      plots()[[1]]  + scale_fill_manual(values = pred_outputs$colors)
    })

    #Output the genomic prediction correlation box plots
    output$pred_violin_plot <- renderPlot({
      plots()[[2]] + scale_fill_manual(values = pred_outputs$colors)
    })

    #Output the prediction tables

    all_GEBVs <- reactive({
      validate(
        need(!is.null(pred_outputs$all_GEBVs), "Upload the input files, set the parameters and click 'run analysis' to access results in this session.")
      )
      pred_outputs$comb_output
    })

    output$pred_all_table <- renderDT({all_GEBVs()}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

    comb_output <- reactive({
      validate(
        need(!is.null(pred_outputs$comb_output), "Upload the input files, set the parameters and click 'run analysis' to access results in this session.")
      )
      pred_outputs$comb_output
    })

    output$pred_acc_table <- renderDT({comb_output()}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

    avg_GEBVs <- reactive({
      validate(
        need(!is.null(pred_outputs$avg_GEBVs), "Upload the input files, set the parameters and click 'run analysis' to access results in this session.")
      )
      pred_outputs$avg_GEBVs
    })

    output$pred_gebvs_table <- renderDT({avg_GEBVs()}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

    #Download files for GP
    output$download_pred_file <- downloadHandler(
      filename = function() {
        paste0("GS-results-", Sys.Date(), ".zip")
      },
      content = function(file) {
        # Temporary files list
        temp_dir <- tempdir()
        temp_files <- c()

        if (!is.null(pred_outputs$avg_GEBVs)) {
          # Create a temporary file for assignments
          gebv_file <- file.path(temp_dir, paste0("GEBVs-", Sys.Date(), ".csv"))
          write.csv(pred_outputs$avg_GEBVs, gebv_file, row.names = FALSE)
          temp_files <- c(temp_files, gebv_file)
        }

        if (!is.null(pred_outputs$comb_output)) {
          # Create a temporary file for BIC data frame
          acc_file <- file.path(temp_dir, paste0("GS-accuracy-statistics-", Sys.Date(), ".csv"))
          write.csv(pred_outputs$comb_output, acc_file, row.names = FALSE)
          temp_files <- c(temp_files, acc_file)
        }

        # Zip files only if there's something to zip
        if (length(temp_files) > 0) {
          zip(file, files = temp_files, extras = "-j") # Using -j to junk paths
        }

        # Optionally clean up
        file.remove(temp_files)
      }
    )

    #Download GP Figures
    output$download_pred_figure <- downloadHandler(

      filename = function() {
        if (input$pred_image_type == "jpeg") {
          paste("GS-", Sys.Date(), ".jpg", sep="")
        } else if (input$pred_image_type == "png") {
          paste("GS-", Sys.Date(), ".png", sep="")
        } else {
          paste("GS-", Sys.Date(), ".tiff", sep="")
        }
      },
      content = function(file) {
        #req(all_plots$pca_2d, all_plots$pca3d, all_plots$scree, input$pca_image_type, input$pca_image_res, input$pca_image_width, input$pca_image_height) #Get the plots
        req(input$pred_figures)

        if (input$pred_image_type == "jpeg") {
          jpeg(file, width = as.numeric(input$pred_image_width), height = as.numeric(input$pred_image_height), res= as.numeric(input$pred_image_res), units = "in")
        } else if (input$pred_image_type == "png") {
          png(file, width = as.numeric(input$pred_image_width), height = as.numeric(input$pred_image_height), res= as.numeric(input$pred_image_res), units = "in")
        } else {
          tiff(file, width = as.numeric(input$pred_image_width), height = as.numeric(input$pred_image_height), res= as.numeric(input$pred_image_res), units = "in")
        }

        # Conditional plotting based on input selection
        if (input$pred_figures == "Violin Plot") {
          req(pred_outputs$violin_plot)

          print(pred_outputs$violin_plot + scale_fill_manual(values = pred_outputs$colors))

        } else if (input$pred_figures == "Box Plot") {
          req(pred_outputs$box_plot)
          #Plot
          print(pred_outputs$box_plot  + scale_fill_manual(values = pred_outputs$colors))

        }

        dev.off()
      }

    )

    output$download_vcf <- downloadHandler(
      filename = function() {
        paste0("BIGapp_VCF_Example_file.vcf.gz")
      },
      content = function(file) {
        ex <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")
        file.copy(ex, file)
      })

    output$download_pheno <- downloadHandler(
      filename = function() {
        paste0("BIGapp_passport_Example_file.csv")
      },
      content = function(file) {
        ex <- system.file("iris_passport_file.csv", package = "BIGapp")
        file.copy(ex, file)
      })
}

## To be copied in the UI
# mod_GSAcc_ui("GSAcc_1")

## To be copied in the server
# mod_GSAcc_server("GSAcc_1")