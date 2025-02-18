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
#' @import shinydisconnect
#'
#'
mod_GS_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      disconnectMessage(
        text = "An unexpected error occurred, please reload the application and check the input file(s).",
        refresh = "Reload now",
        background = "white",
        colour = "grey",
        overlayColour = "grey",
        overlayOpacity = 0.3,
        refreshColour = "purple"
      ),
      column(width = 3,
             box(title="Inputs", width = 12, collapsible = TRUE, collapsed = FALSE, status = "info", solidHeader = TRUE,
                 "* Required",
                 fileInput(ns("pred_known_file"), "Choose VCF File*", accept = c(".csv",".vcf",".gz")),
                 fileInput(ns("pred_trait_file"), "Choose Trait File*", accept = ".csv"),
                 numericInput(ns("pred_est_ploidy"), "Species Ploidy*", min = 1, value = NULL),
                 virtualSelectInput(
                   inputId = ns("pred_trait_info2"),
                   label = "Select Trait*",
                   choices = NULL,
                   showValueAsTags = TRUE,
                   search = TRUE,
                   multiple = TRUE
                 ),
                 virtualSelectInput(
                   inputId = ns("pred_fixed_info2"),
                   label = span("Select Fixed Effects", bs4Badge("beta", position = "right", color = "success")),
                   choices = NULL,
                   showValueAsTags = TRUE,
                   search = TRUE,
                   multiple = TRUE
                 ),
                 conditionalPanel(
                   condition = "input.pred_fixed_info2.length > 0", ns = ns,
                   div(
                     "(unselected will be considered covariates)",
                     virtualSelectInput(
                       inputId = ns("pred_fixed_cat2"),
                       label = "Select Categorical Fixed Effects",
                       choices = NULL,
                       showValueAsTags = TRUE,
                       search = TRUE,
                       multiple = TRUE
                     )
                   )
                 ),
                 actionButton(ns("prediction_est_start"), "Run Analysis"),
                 div(style="display:inline-block; float:right",dropdownButton(
                   HTML("<b>Input files</b>"),
                   p(downloadButton(ns('download_vcf'),""), "VCF Example File"),
                   p(downloadButton(ns('download_pheno'),""), "Trait Example File"),
                   p(downloadButton(ns('download_vcfp'), ""), "Download Prediction VCF Example File"),hr(),
                   p(HTML("<b>Parameters description:</b>"), actionButton(ns("goPar"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>Results description:</b>"), actionButton(ns("goRes"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>How to cite:</b>"), actionButton(ns("goCite"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   actionButton(ns("pred_summary"), "Summary"),
                   circle = FALSE,
                   status = "warning",
                   icon = icon("info"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see info!")
                 )),
                 tags$hr(style="border-color: #d3d3d3; margin-top: 20px; margin-bottom: 20px;"),  # Lighter grey line
                 div(style="text-align: left; margin-top: 10px;",
                     actionButton(ns("advanced_options_pred"),
                                  label = HTML(paste(icon("cog", style = "color: #007bff;"), "Advanced Options (beta)")),
                                  style = "background-color: transparent; border: none; color: #007bff; font-size: smaller; text-decoration: underline; padding: 0;"
                     )
                 )

             )
      ),

      column(width = 6,
             box(title = "Results", status = "info", solidHeader = FALSE, width = 12, height = 600, maximizable = T,
                 bs4Dash::tabsetPanel(
                   tabPanel("Predicted Pheno Table", DTOutput(ns("pred_trait_table")), style = "overflow-y: auto; height: 500px"),
                   tabPanel("EBVs Table", DTOutput(ns("pred_gebvs_table2")),style = "overflow-y: auto; height: 500px")

                 )
             )
      ),

      column(width = 3,
             valueBoxOutput("shared_snps", width = NULL),
             box(title = "Status", width = 12, collapsible = TRUE, status = "info",
                 progressBar(id = ns("pb_gp"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
             ),
             box(title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
                 div(style="display:inline-block; float:left",dropdownButton(
                   tags$h3("Save Files"),
                   fluidRow(
                     downloadButton(ns("download_pred_results_file"), "Save Files")),
                   circle = FALSE,
                   status = "danger",
                   icon = icon("floppy-disk"), width = "300px", label = "Save",
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
mod_GS_server <- function(input, output, session, parent_session){

  ns <- session$ns

  # Help links
  observeEvent(input$goPar, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "Genomic_Prediction_tabset",
                      selected = "Genomic_Prediction_par")
    # expand specific box
    updateBox(id = "Genomic_Prediction_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goRes, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "Genomic_Prediction_tabset",
                      selected = "Genomic_Prediction_results")
    # expand specific box
    updateBox(id = "Genomic_Prediction_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goCite, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "Genomic_Prediction_tabset",
                      selected = "Genomic_Prediction_cite")
    # expand specific box
    updateBox(id = "Genomic_Prediction_box", action = "toggle", session = parent_session)
  })

  #Default model choices
  advanced_options_pred <- reactiveValues(
    pred_model = "GBLUP",
    pred_matrix = "Gmatrix",
    pred_est_file = NULL,
    ped_file = NULL
  )

  pred_outputs <- reactiveValues(corr_output = NULL,
                                 comb_output = NULL,
                                 all_GEBVs = NULL,
                                 avg_GEBVs = NULL)

  #List the ped file name if previously uploaded
  output$uploaded_file_name <- renderText({
    if (!is.null(advanced_options_pred$ped_file)) {
      paste("Previously uploaded file:", advanced_options_pred$ped_file$name)
    } else {
      ""  # Return an empty string if no file has been uploaded
    }
  })

  output$uploaded_file_name_pred <- renderText({
    if (!is.null(advanced_options_pred$pred_est_file)) {
      paste("Previously uploaded file:", advanced_options_pred$pred_est_file$name)
    } else {
      ""  # Return an empty string if no file has been uploaded
    }
  })

  #UI popup window for input
  observeEvent(input$advanced_options_pred, {
    showModal(modalDialog(
      title = "Advanced Options (beta)",
      selectInput(
        inputId = ns('pred_model'),
        label = 'Method Choice',
        choices = c("GBLUP"),
        selected = advanced_options_pred$pred_model  # Initialize with stored value
      ),
      conditionalPanel(
        condition = "input.pred_model == 'GBLUP'", ns = ns,
        div(
          selectInput(
            inputId = ns('pred_matrix'),
            label = 'Relationship Matrix Choice',
            #choices = c("Gmatrix", "Amatrix", "Hmatrix"),
            choices = c("Gmatrix"),
            selected = advanced_options_pred$pred_matrix  # Initialize with stored value
          )
        )
      ),
      conditionalPanel(
        condition = "input.pred_matrix != 'Amatrix'", ns = ns,
        div(
          fileInput(ns("pred_est_file"), "Choose Prediction VCF", accept = c(".vcf",".gz")),
          conditionalPanel(
            condition = "output.uploaded_file_name_pred !== ''", # Show only if there's content
            textOutput(ns("uploaded_file_name_pred"))  # Display the uploaded file name
          )
        )
      ),
      conditionalPanel(
        condition = "input.pred_matrix != 'Gmatrix'", ns = ns,
        div(
          fileInput(ns("ped_file"), "Choose Pedigree File", accept = ".csv"),
          conditionalPanel(
            condition = "output.uploaded_file_name !== ''", # Show only if there's content
            textOutput(ns("uploaded_file_name"))  # Display the uploaded file name
          )
        )
      ),
      footer = tagList(
        modalButton("Close"),
        actionButton(ns("save_advanced_options_pred"), "Save")
      )
    ))
  })



  #Close popup window when user "saves options"
  observeEvent(input$save_advanced_options_pred, {
    advanced_options_pred$pred_model <- input$pred_model
    advanced_options_pred$pred_matrix <- input$pred_matrix
    advanced_options_pred$pred_est_file <- input$pred_est_file
    advanced_options_pred$ped_file <- input$ped_file
    # Save other inputs as needed

    removeModal()  # Close the modal after saving
  })

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
    updateVirtualSelect("pred_fixed_info2", choices = trait_var2, session = session)
    updateVirtualSelect("pred_trait_info2", choices = trait_var2, session = session)
  })

  observeEvent(input$pred_fixed_info2, {
    updateVirtualSelect("pred_fixed_cat2", choices = input$pred_fixed_info2, session = session)
  })

  #2) Error check for prediction and save input files
  continue_prediction2 <- reactiveVal(NULL)
  pred_inputs2 <- reactiveValues(
    pheno_input = NULL,
    train_geno_input = NULL,
    est_geno_input = NULL,
    shared_snps = NULL,
    pred_genos = NULL,
    pred_geno_pheno = NULL,
    matrix = NULL
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

    toggleClass(id = "pred_est_ploidy", class = "borderred", condition = (is.na(input$pred_est_ploidy) | is.null(input$pred_est_ploidy)))

    if (is.null(input$pred_known_file$datapath) | is.null(input$pred_trait_file$datapath)) {
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
    req(input$pred_known_file$datapath, input$pred_trait_file$datapath, input$pred_est_ploidy)

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
        title = "Missing input!",
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
      train_vcf <- vcfR::read.vcfR(train_geno_path, verbose = FALSE)
      if (is.null(est_geno_path) || is.na(est_geno_path)){
        est_vcf <- NULL
      } else {
        est_vcf <- vcfR::read.vcfR(est_geno_path, verbose = FALSE)
      }

      #Get number of SNPs
      #pred_inputs$pred_snps <- nrow(vcf)

      #Extract GT
      if (is.null(est_vcf)) {
        est_geno <- NULL
      } else {
        est_geno <- extract.gt(est_vcf, element = "GT")
        est_geno <- apply(est_geno, 2, convert_to_dosage)
        class(est_geno) <- "numeric"
      }
      train_geno <- extract.gt(train_vcf, element = "GT")
      train_geno <- apply(train_geno, 2, convert_to_dosage)
      class(train_geno) <- "numeric"
      rm(train_vcf)
      rm(est_vcf)

    } else {

      # If condition is met, show notification toast
      shinyalert(
        title = "File Error",
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
    if (is.null(est_geno)) {
      est_geno_adj_init <- NULL
    } else {
      est_geno_adj_init <- convert_genotype(est_geno, as.numeric(ploidy))
    }

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
        title = "Data Mismatch",
        text = "All genotyped samples were missing from the phenotype file",
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
        title = "Data Check",
        text = paste0((length(common_ids))," samples in VCF File have trait information"),
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

    #Status
    updateProgressBar(session = session, id = "pb_gp", value = 40, title = "Generating Matrices")

    #Create relationship matrix depending on the input VCF files
    if (is.null(advanced_options_pred$pred_est_file)) {
      #Subset phenotype file by common IDs
      pheno2 <- pheno2[common_ids, ]

      #Matrix
      Kin_mat <- A.mat(t(train_geno_adj_init), max.missing=NULL,impute.method="mean",return.imputed=FALSE)

    } else{
      #Subset phenotype file and training file by common IDs
      pheno2 <- pheno2[common_ids, ]

      #Match training and testing genotype file SNPs
      common_markers <- intersect(rownames(train_geno_adj_init), rownames(est_geno_adj_init))
      #first remove samples from training genotype matrix that are not in the phenotype file
      train_geno_adj <- train_geno_adj_init[common_markers, common_ids]
      #Merge the training and prediction genotype matrices
      est_geno_adj_init <- est_geno_adj_init[common_markers, ]
      train_pred_geno <- cbind(train_geno_adj, est_geno_adj_init)

      #Matrix
      Kin_mat <- A.mat(t(train_pred_geno), max.missing=NULL,impute.method="mean",return.imputed=FALSE)

    }

    #Save to reactive values
    #pred_inputs2$shared_snps <- length(common_markers)
    pred_inputs2$matrix <- Kin_mat
    pred_inputs2$pheno_input <- pheno2
  })

  #3) Analysis only proceeds once continue_prediction is converted to TRUE
  observe({

    req(continue_prediction2(),pred_inputs2$pheno_input, pred_inputs2$matrix)

    # Stop analysis if cancel was selected
    if (isFALSE(continue_prediction2())) {
      return()
    }

    #Variables
    ploidy <- as.numeric(input$pred_est_ploidy)
    gmat <- pred_inputs2$matrix
    pheno2 <- pred_inputs2$pheno_input
    traits <- input$pred_trait_info2
    # Assign fixed_cat based on input$pred_fixed_cat2
    if (!is.null(input$pred_fixed_cat2) && length(input$pred_fixed_cat2) > 0) {
      fixed_cat <- input$pred_fixed_cat2
    } else {
      fixed_cat <- NULL
    }
    # Assign fixed_cov based on conditions
    if (!is.null(fixed_cat) && length(fixed_cat) == length(input$pred_fixed_info2)) {
      fixed_cov <- NULL
    } else if (length(input$pred_fixed_info2) > 0 && is.null(fixed_cat)) {
      fixed_cov <- input$pred_fixed_info2
    } else {
      fixed_cov <- setdiff(input$pred_fixed_info2, input$pred_fixed_cat2)
    }
    #fixed_cov <- setdiff(input$pred_fixed_info2, input$pred_fixed_cat2)
    cores <- 1
    total_population <- ncol(gmat)
    #train_size <- floor(percentage / 100 * total_population)

    #Status
    updateProgressBar(session = session, id = "pb_gp", value = 90, title = "Generating Results")

    #initialize dataframe
    results_GEBVs <- matrix(nrow = ncol(gmat), ncol = length(traits) + 1)
    results_TRAITs <- matrix(nrow = ncol(gmat), ncol = length(traits) + 1)
    colnames(results_TRAITs) <- c("Sample",paste0(traits))  # Set the column names to be the traits
    colnames(results_GEBVs) <- c("Sample",paste0(traits))  # Set the column names to be the traits

    #GBLUP for each trait
    for (trait_idx in 1:length(traits)) {
      traitpred <- kin.blup(data = pheno2,
                            geno = colnames(pheno2)[1],
                            pheno = traits[trait_idx],
                            fixed = fixed_cat,
                            covariate = fixed_cov,
                            K=gmat,
                            n.core = cores)

      results_GEBVs[, (trait_idx+1)] <- traitpred$g
      results_TRAITs[, (trait_idx+1)] <- traitpred$pred
    }
    #Organize dataframes
    results_GEBVs[,1] <- row.names(data.frame(traitpred$g))
    results_TRAITs[,1] <- row.names(data.frame(traitpred$pred))

    #Label the samples that already had phenotype information
    results_GEBVs <- data.frame(results_GEBVs)
    results_TRAITs <- data.frame(results_TRAITs)
    exists_in_df <- results_GEBVs[[1]] %in% pheno2[[1]]
    results_GEBVs <- cbind(results_GEBVs[1], "w/Pheno" = exists_in_df, results_GEBVs[-1])
    results_TRAITs <- cbind(results_TRAITs[1], "w/Pheno" = exists_in_df, results_TRAITs[-1])

    #Status
    updateProgressBar(session = session, id = "pb_gp", value = 100, title = "Finished!")

    ##Make output tables depending on 1 or 2 VCF/pedigree files used.
    #GEBVs
    if (!is.null(advanced_options_pred$pred_est_file)) {
      # Subset rows where 'w/Pheno' is FALSE and drop the 'w/Pheno' column
      pred_outputs2$all_GEBVs <- results_GEBVs[results_GEBVs$`w/Pheno` == FALSE, !names(results_GEBVs) %in% "w/Pheno"]
    } else{
      pred_outputs2$all_GEBVs <- results_GEBVs
    }

    #BLUPs of genetic values
    if (!is.null(advanced_options_pred$pred_est_file)) {
      # Subset rows where 'w/Pheno' is FALSE and drop the 'w/Pheno' column
      pred_outputs2$trait_output <- results_TRAITs[results_TRAITs$`w/Pheno` == FALSE, !names(results_TRAITs) %in% "w/Pheno"]
    } else{
      pred_outputs2$trait_output <- results_TRAITs
    }

    #End the event
    continue_prediction2(NULL)
  })

  #Output the prediction tables
  all_GEBVs <- reactive({
    validate(
      need(!is.null(pred_outputs2$all_GEBVs), "Upload the input files, set the parameters and click 'run analysis' to access results in this session.")
    )
    pred_outputs2$all_GEBVs
  })

  #GEBVs from all iterations/folds
  output$pred_gebvs_table2 <- renderDT({all_GEBVs()}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

  trait_output <- reactive({
    validate(
      need(!is.null(pred_outputs2$trait_output), "Upload the input files, set the parameters and click 'run analysis' to access results in this session.")
    )
    pred_outputs2$trait_output
  })

  #GEBVs from all iterations/folds
  output$pred_trait_table <- renderDT({trait_output()}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

  output$download_vcft <- downloadHandler(
    filename = function() {
      paste0("BIGapp_Training_VCF_Example_file.vcf.gz")
    },
    content = function(file) {
      ex <- system.file("test-dose.vcf.gz", package = "BIGapp")
      file.copy(ex, file)
    })

  output$download_vcfp <- downloadHandler(
    filename = function() {
      paste0("BIGapp_Predict_VCF_Example_file.vcf")
    },
    content = function(file) {
      ex <- system.file("test-dose-use-for-prediction.vcf", package = "BIGapp")
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

  #Download files for GP
  output$download_pred_results_file <- downloadHandler(
    filename = function() {
      paste0("Prediction-results-", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Temporary files list
      temp_dir <- tempdir()
      temp_files <- c()

      # Create a temporary file for data frames
      ebv_file <- file.path(temp_dir, paste0("GS-EBV-predictions-", Sys.Date(), ".csv"))
      write.csv(pred_outputs2$all_GEBVs, ebv_file, row.names = FALSE)
      temp_files <- c(temp_files, ebv_file)

      trait_file <- file.path(temp_dir, paste0("GS-Phenotype-predictions-", Sys.Date(), ".csv"))
      write.csv(pred_outputs2$trait_output, trait_file, row.names = FALSE)
      temp_files <- c(temp_files, trait_file)

      # Zip files only if there's something to zip
      if (length(temp_files) > 0) {
        zip(file, files = temp_files, extras = "-j") # Using -j to junk paths
      }

      # Optionally clean up
      file.remove(temp_files)
    }
  )

  ##Summary Info
  pred_summary_info <- function() {
    # Handle possible NULL values for inputs
    dosage_file_name <- if (!is.null(input$pred_known_file$name)) input$pred_known_file$name else "No file selected"
    est_file_name <- if (!is.null(input$pred_est_file$name)) input$pred_est_file$name else "No file selected"
    passport_file_name <- if (!is.null(input$pred_trait_file$name)) input$pred_trait_file$name else "No file selected"
    selected_ploidy <- if (!is.null(input$pred_est_ploidy)) as.character(input$pred_est_ploidy) else "Not selected"

    # Print the summary information
    cat(
      "BIGapp Selection Summary\n",
      "\n",
      paste0("Date: ", Sys.Date()), "\n",
      paste("R Version:", R.Version()$version.string), "\n",
      "\n",
      "### Input Files ###\n",
      "\n",
      paste("Input Genotype File 1:", dosage_file_name), "\n",
      paste("Input Genotype File 2:", est_file_name), "\n",
      paste("Input Passport File:", passport_file_name), "\n",
      "\n",
      "### User Selected Parameters ###\n",
      "\n",
      paste("Selected Ploidy:", selected_ploidy), "\n",
      paste("Selected Trait(s):", input$pred_trait_info2), "\n",
      paste("Selected Fixed Effects:", input$pred_fixed_info2), "\n",
      #paste("Selected Model:", input$pred_fixed_info2), "\n",
      #paste("Selected Matrix:", input$pred_fixed_info2), "\n",
      "\n",
      "### R Packages Used ###\n",
      "\n",
      paste("BIGapp:", packageVersion("BIGapp")), "\n",
      paste("AGHmatrix:", packageVersion("AGHmatrix")), "\n",
      paste("ggplot2:", packageVersion("ggplot2")), "\n",
      paste("rrBLUP:", packageVersion("rrBLUP")), "\n",
      paste("vcfR:", packageVersion("vcfR")), "\n",
      paste("dplyr:", packageVersion("dplyr")), "\n",
      paste("tidyr:", packageVersion("tidyr")), "\n",
      sep = ""
    )
  }

  # Popup for analysis summary
  observeEvent(input$pred_summary, {
    showModal(modalDialog(
      title = "Summary Information",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close"),
        downloadButton("download_pred_info", "Download")
      ),
      pre(
        paste(capture.output(pred_summary_info()), collapse = "\n")
      )
    ))
  })


  # Download Summary Info
  output$download_pred_info <- downloadHandler(
    filename = function() {
      paste("pred_summary_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      # Write the summary info to a file
      writeLines(paste(capture.output(pred_summary_info()), collapse = "\n"), file)
    }
  )
}

## To be copied in the UI
# mod_GS_ui("GS_1")

## To be copied in the server
# mod_GS_server("GS_1")
