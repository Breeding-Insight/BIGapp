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
#' @import shinydisconnect
#'
mod_GSAcc_ui <- function(id){
  ns <- NS(id)
  tagList(
    # Add GWAS content here
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
                 fileInput(ns("pred_file"), "Choose VCF File*", accept = c(".csv",".vcf",".gz")),
                 fileInput(ns("trait_file"), "Choose Trait File*", accept = ".csv"),
                 numericInput(ns("pred_ploidy"), "Species Ploidy*", min = 1, value = NULL),
                 numericInput(ns("pred_cv"), "Iterations*", min = 1, max=20, value = 5),
                 virtualSelectInput(
                   inputId = ns("pred_trait_info"),
                   label = "Select Trait(s)*:",
                   choices = NULL,
                   showValueAsTags = TRUE,
                   search = TRUE,
                   multiple = TRUE
                 ),
                 virtualSelectInput(
                   inputId = ns("pred_fixed_info"),
                   label = span("Select Fixed Effects", bs4Badge("beta", position = "right", color = "success")),
                   choices = NULL,
                   showValueAsTags = TRUE,
                   search = TRUE,
                   multiple = TRUE
                 ),
                 conditionalPanel(
                   condition = "input.pred_fixed_info.length > 0", ns = ns,
                   "(unselected will be considered covariates)",
                   div(
                     virtualSelectInput(
                       inputId = ns("pred_fixed_cat"),
                       label = "Select Categorical Fixed Effects",
                       choices = NULL,
                       showValueAsTags = TRUE,
                       search = TRUE,
                       multiple = TRUE
                     )
                   )
                 ),
                 actionButton(ns("prediction_start"), "Run Analysis"),
                 div(style="display:inline-block; float:right", dropdownButton(
                   HTML("<b>Input files</b>"),
                   p(downloadButton(ns('download_vcf'),""), "VCF Example File"),
                   p(downloadButton(ns('download_pheno'),""), "Trait Example File"), hr(),
                   p(HTML("<b>Parameters description:</b>"), actionButton(ns("goPar"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>Results description:</b>"), actionButton(ns("goRes"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>How to cite:</b>"), actionButton(ns("goCite"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   actionButton(ns("predAcc_summary"), "Summary"),
                   circle = FALSE,
                   status = "warning",
                   icon = icon("info"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see info!")
                 )),
                 tags$hr(style="border-color: #d3d3d3; margin-top: 20px; margin-bottom: 20px;"),  # Lighter grey line
                 div(style="text-align: left; margin-top: 10px;",
                     actionButton(ns("advanced_options"),
                                  label = HTML(paste(icon("cog", style = "color: #007bff;"), "Advanced Options (beta)")),
                                  style = "background-color: transparent; border: none; color: #007bff; font-size: smaller; text-decoration: underline; padding: 0;"
                     )
                 )
             )
      ),

      column(width = 6,
             box(
               title = "Plots", status = "info", solidHeader = FALSE, width = 12, height = 600, maximizable = T,
               bs4Dash::tabsetPanel(
                 tabPanel("Violin Plot", plotOutput(ns("pred_violin_plot"), height = "500px")),
                 tabPanel("Box Plot", plotOutput(ns("pred_box_plot"), height = "500px")),
                 tabPanel("P.A. Table", DTOutput(ns("pred_acc_table")), style = "overflow-y: auto; height: 500px")
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
#' @importFrom rrBLUP mixed.solve A.mat kin.blup
#' @importFrom stats cor
#' @importFrom shinyalert shinyalert
#' @importFrom AGHmatrix Gmatrix Amatrix Hmatrix
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @noRd
mod_GSAcc_server <- function(input, output, session, parent_session){

  ns <- session$ns

  # Help links
  observeEvent(input$goPar, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "Predictive_Ability_tabset",
                      selected = "Predictive_Ability_par")
    # expand specific box
    updateBox(id = "Predictive_Ability_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goRes, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "Predictive_Ability_tabset",
                      selected = "Predictive_Ability_results")
    # expand specific box
    updateBox(id = "Predictive_Ability_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goCite, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "Predictive_Ability_tabset",
                      selected = "Predictive_Ability_cite")
    # expand specific box
    updateBox(id = "Predictive_Ability_box", action = "toggle", session = parent_session)
  })

  #Default model choices
  advanced_options <- reactiveValues(
    pred_model = "GBLUP",
    pred_matrix = "Gmatrix",
    ped_file = NULL
  )

  pred_outputs <- reactiveValues(corr_output = NULL,
                                 comb_output = NULL,
                                 all_GEBVs = NULL,
                                 avg_GEBVs = NULL)

  #List the ped file name if previously uploaded
  output$uploaded_file_name <- renderText({
    if (!is.null(advanced_options$ped_file)) {
      paste("Previously uploaded file:", advanced_options$ped_file$name)
    } else {
      ""  # Return an empty string if no file has been uploaded
    }
  })

  #UI popup window for input
  observeEvent(input$advanced_options, {
    showModal(modalDialog(
      title = "Advanced Options (beta)",
      selectInput(
        inputId = ns('pred_model'),
        label = 'Method Choice',
        choices = c("GBLUP"),
        selected = advanced_options$pred_model  # Initialize with stored value
      ),
      conditionalPanel(
        condition = "input.pred_model == 'GBLUP'", ns = ns,
        div(
          selectInput(
            inputId = ns('pred_matrix'),
            label = 'Relationship Matrix Choice',
            choices = c("Gmatrix", "Amatrix", "Hmatrix"),
            selected = advanced_options$pred_matrix  # Initialize with stored value
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
        actionButton(ns("save_advanced_options"), "Save")
      )
    ))
  })



  #Close popup window when user "saves options"
  observeEvent(input$save_advanced_options, {
    advanced_options$pred_model <- input$pred_model
    advanced_options$pred_matrix <- input$pred_matrix
    advanced_options$ped_file <- input$ped_file
    # Save other inputs as needed

    removeModal()  # Close the modal after saving
  })



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

  colors <- reactiveValues(colors = NULL)
  values_boxes <- reactiveValues(pred_snps = 0, pred_geno_pheno = 0)

  #Reactive boxes
  output$pred_snps <- renderValueBox({
    valueBox(
      value = values_boxes$pred_snps,
      subtitle = "SNPs in Genotype File",
      icon = icon("dna"),
      color = "info"
    )
  })

  output$pred_geno <- renderValueBox({
    valueBox(
      value = values_boxes$pred_geno_pheno,
      subtitle = "Samples with Phenotype Information",
      icon = icon("location-dot"),
      color = "info"
    )
  })

  observe({
    # Update colors based on input
    colors$colors <- assign_colors(input$pred_color_select)
  })

  observeEvent(input$pred_fixed_info, {
    updateVirtualSelect("pred_fixed_cat", choices = input$pred_fixed_info, session = session)
  })

  #2) Error check for prediction and save input files
  pred_inputs <- eventReactive(input$prediction_start,{

    toggleClass(id = "pred_ploidy", class = "borderred", condition = (is.na(input$pred_ploidy) | is.null(input$pred_ploidy)))

    if(is.null(advanced_options$pred_matrix)) advanced_options$pred_matrix <- "none_selected"
    if (((is.null(input$pred_file$datapath) &  advanced_options$pred_matrix != "Amatrix") |
         (is.null(advanced_options$ped_file$datapath) &  advanced_options$pred_matrix == "Amatrix")) |
        is.null(input$trait_file$datapath)) {
      shinyalert(
        title = "Missing input!",
        text = "Upload VCF or a pedigree file and the phenotype file",
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
      return()
    }

    #Status
    updateProgressBar(session = session, id = "pb_prediction", value = 5, title = "Checking input files")

    #Variables
    pheno <- read.csv(input$trait_file$datapath, header = TRUE, check.names = FALSE)
    row.names(pheno) <- pheno[,1]
    # Assuming the first column in Pheno contains the matching IDs
    ids_pheno <- pheno[, 1]

    #Make sure at least one trait was input
    if (length(input$pred_trait_info) == 0) {

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

    pred_inputs <- list(
      pheno_input = NULL,
      geno_input = NULL,
      ped_input = NULL
    )

    #Getting genotype matrix
    #Geno.file conversion if needed
    if(!is.null(input$pred_file$datapath)){
      geno_snps <- read_geno_file(input$pred_file$datapath, requires = "GT")
      geno <- geno_snps[[1]]
      values_boxes$pred_snps <- geno_snps[[2]]

      #Check that the ploidy entered is correct
      if (input$pred_ploidy != max(geno, na.rm = TRUE)) {
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

      #Make sure the trait file and genotype file are in the same order
      # Column names for geno (assuming these are the individual IDs)
      colnames_geno <- colnames(geno)

      # Find common identifiers
      common_ids <- intersect(colnames_geno, ids_pheno)
      #Get number of id
      values_boxes$pred_geno_pheno <- length(common_ids)

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
      } else {
        if (length(common_ids) < length(colnames_geno))
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
        if (length(common_ids) < length(ids_pheno))
          shinyalert(
            title = "Data Mismatch",
            text = paste0((length(ids_pheno)-length(common_ids))," samples were removed for not having genotypic information"),
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
      }

      # Subset and reorder geno and pheno to ensure they only contain and are ordered by common IDs
      geno_adj <- geno[, common_ids]  # Assuming that the columns can be directly indexed by IDs
      pheno <- pheno[match(common_ids, ids_pheno), ] # If there is pheno but not geno, the sample is also discarded
    } else geno_adj <- NULL

    # Check pedigree
    #Import pedigree file, where pedigree data name (3-column way format). Unknown value should be equal 0
    if(!is.null(advanced_options$ped_file$datapath) & (advanced_options$pred_matrix == "Amatrix" | advanced_options$pred_matrix == "Hmatrix")){
      ped <- read.csv(advanced_options$ped_file$datapath, check.names = FALSE, colClasses = "factor")
      colnames(ped) <- c("Ind", "P1", "P2")
      #Convert NAs to 0
      ped[is.na(ped)] <- 0

      common_ped <- intersect(ped$Ind, pheno[,1])
      #Throw an error if there are less matching samples in the phenotype file than the pedigree file
      if (length(common_ped) == 0) {
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
      } else {
        rm_unr <- remove_unrelated(ped, samples_with_trait_info = pheno[,1])
        extended_ped <- rm_unr[[1]]
        gen <- rm_unr[[2]]
        cat(paste0("You have pedigree information until the ", gen,"th generation\n"))

        if (length(common_ped) < length(ids_pheno)){
          shinyalert(
            title = "Data Mismatch",
            text = paste0((length(pheno[,1])-length(common_ped))," samples were removed from the phenotype data for not having pedigree information"),
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
          if(length(which(!pheno[,1] %in% extended_ped$Ind)) > 0) pheno <- pheno[-which(!pheno[,1] %in% extended_ped$Ind),]
          if(!is.null(geno_adj) & length(which(!colnames(geno_adj) %in% extended_ped$Ind)) > 0) geno_adj <- geno_adj[,-which(!colnames(geno_adj) %in% extended_ped$Ind)]
        }
        if (length(ped$Ind) > length(extended_ped$Ind))
          shinyalert(
            title = "Data Mismatch",
            text = paste0((length(ped$Ind)-length(extended_ped$Ind))," samples in the pedigree file were unrelated to the samples with phenotype information. They were removed from the analysis."),
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

        ped_temp <- tempfile()
        ped_temp_file <- extended_ped
        colnames(ped_temp_file) <- c("id", "sire", "dam")
        write.table(ped_temp_file, file = ped_temp)
        ped_check <- BIGr::check_ped(ped_temp)
        if(dim(ped_check$repeated_ids)[1] != 0){
          shinyalert(
            title = "Oops",
            text = "Check for repeated IDs in the pedigree file",
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
        if(dim(ped_check$messy_parents)[1] != 0){
          shinyalert(
            title = "Oops",
            text = paste("We found inconsistencies in the pedigree file for the individuals:", paste0(ped_check$messy_parents$id, collapse = ", ")),
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
        }
      }
    } else extended_ped <- NULL

    #Status
    updateProgressBar(session = session, id = "pb_prediction", value = 8, title = "Inputs checked!")

    ## Make ouput as checked inputs pred_inputs
    pred_inputs$pheno_input <- pheno
    pred_inputs$geno_input <- geno_adj
    pred_inputs$ped_input <- extended_ped
    pred_inputs
  })

  observeEvent(pred_inputs(),{

    # Convert genotype matrix according to ploidy and model used
    if(!is.null(pred_inputs()$geno_input)){
      geno_formated <- format_geno_matrix(pred_inputs()$geno_input,advanced_options$pred_model, advanced_options$pred_matrix, input$pred_ploidy)
    } else geno_formated <- NULL

    #Status
    updateProgressBar(session = session, id = "pb_prediction", value = 10, title = paste("Genotype matrix formatted for", advanced_options$pred_model, advanced_options$pred_matrix))

    results <- run_predictive_model(geno = geno_formated,
                                    pheno = pred_inputs()$pheno_input,
                                    selected_traits = input$pred_trait_info,
                                    predictive_model = advanced_options$pred_model,
                                    relationship_matrix_type = advanced_options$pred_matrix,
                                    pedigree = pred_inputs()$ped_input,
                                    fixed_effects = input$pred_fixed_info,
                                    categorical_fixed_effects = input$pred_fixed_cat,
                                    ploidy = input$pred_ploidy,
                                    cores = input$pred_cores,
                                    cycles = input$pred_cv,
                                    folds = 5, session = session)

    updateProgressBar(session = session, id = "pb_prediction", value = 90, title = "Cross validation concluded")

    #Save to reactive value
    pred_outputs$corr_output <- results$PredictionAccuracy
    pred_outputs$all_GEBVs <- results$GEBVs

    # Convert trait columns to numeric
    results$GEBVs <- results$GEBVs %>%
      mutate(across(all_of(input$pred_trait_info), ~ as.numeric(.x)))

    # Calculate the average value for each column in the traits list for each SampleID, ignoring Iter and Fold
    average_gebvs_df <- results$GEBVs %>%
      group_by(Sample) %>%
      summarize(across(all_of(input$pred_trait_info), mean, na.rm = TRUE))

    pred_outputs$avg_GEBVs <- average_gebvs_df

    columns <- setdiff(colnames(results$PredictionAccuracy), c("Iter","Fold"))
    average_accuracy_df <- results$PredictionAccuracy %>%
      group_by(Iter) %>%
      summarize(across(all_of(columns), mean, na.rm = TRUE))

    pred_outputs$comb_output <- average_accuracy_df
    #Status
    updateProgressBar(session = session, id = "pb_prediction", value = 100, title = "Finished!")

    #pred_outputs
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

    plots <- list(box_plot = NULL, violin_plot = NULL)
    #This can be adapted if we start comparing more than one GP model
    #Also consider a violin plot to show each cor value
    #plot <- ggplot(df_long, aes(x = factor(Trait), y = Correlation, fill = "red"), fill = "red") +
    plots$box_plot <- ggplot(df_long, aes(x = "rrBLUP", y = Correlation, fill = "red"), fill = "red") +
      #geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 0.7, outlier.size = 0.2) +
      geom_boxplot() +
      facet_wrap(~ Trait, nrow = 1) +  # Facet by trait, allowing different y-scales
      labs(title = "Predictive Ability by Trait",
           x = " ",
           y = "Predictive Ability") +
      #theme_minimal() +                      # Using a minimal theme
      theme(legend.position = "none",
            strip.text = element_text(size = 12),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2),
            strip.text.x = element_text(face = "bold"),
            axis.text.x.bottom = element_blank(),
            axis.ticks.x.bottom = element_blank())

    plots$violin_plot <- ggplot(df_long, aes(x = "rrBLUP", y = Correlation, fill = "red")) +
      geom_violin(trim = TRUE) +  # Add violin plot
      geom_point(position = position_jitter(width = 0.1), color = "black", size = 1.5) +  # Add jittered points
      facet_wrap(~ Trait, nrow = 1) +  # Facet by trait, allowing different y-scales
      labs(title = "Predictive Ability by Trait",
           x = " ",  # x-label is blank because it's not relevant per facet
           y = "Predictive Ability") +
      theme(legend.position = "none",
            strip.text = element_text(size = 12),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2),
            strip.text.x = element_text(face = "bold"),
            axis.text.x.bottom = element_blank(),
            axis.ticks.x.bottom = element_blank())
    plots
  })

  #Output the genomic prediction correlation box plots
  output$pred_box_plot <- renderPlot({
    plots()$box_plot  + scale_fill_manual(values = colors$colors)
  })

  #Output the genomic prediction correlation box plots
  output$pred_violin_plot <- renderPlot({
    plots()$violin_plot + scale_fill_manual(values = colors$colors)
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

  #Download files for GP
  output$download_pred_file <- downloadHandler(
    filename = function() {
      paste0("GS-results-", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Temporary files list
      temp_dir <- tempdir()
      temp_files <- c()

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
        req(plots()$violin_plot)

        print(plots()$violin_plot + scale_fill_manual(values = colors$colors))

      } else if (input$pred_figures == "Box Plot") {
        req(plots()$box_plot)
        #Plot
        print(plots()$box_plot  + scale_fill_manual(values = colors$colors))

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

  ##Summary Info
  predAcc_summary_info <- function() {
    # Handle possible NULL values for inputs
    dosage_file_name <- if (!is.null(input$pred_file$name)) input$pred_file$name else "No file selected"
    passport_file_name <- if (!is.null(input$trait_file$name)) input$trait_file$name else "No file selected"
    selected_ploidy <- if (!is.null(input$pred_ploidy)) as.character(input$pred_ploidy) else "Not selected"

    # Print the summary information
    cat(
      "BIGapp Selection Model CV Summary\n",
      "\n",
      paste0("Date: ", Sys.Date()), "\n",
      paste("R Version:", R.Version()$version.string), "\n",
      "\n",
      "### Input Files ###\n",
      "\n",
      paste("Input Genotype File:", dosage_file_name), "\n",
      paste("Input Passport File:", passport_file_name), "\n",
      "\n",
      "### User Selected Parameters ###\n",
      "\n",
      paste("Selected Ploidy:", selected_ploidy), "\n",
      paste("Selected Trait(s):", input$pred_trait_info), "\n",
      paste("Selected Fixed Effects:", input$pred_fixed_info), "\n",
      paste("Selected Model:", advanced_options$pred_model), "\n",
      paste("Selected Matrix:", advanced_options$pred_matrix), "\n",
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
  observeEvent(input$predAcc_summary, {
    showModal(modalDialog(
      title = "Summary Information",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close"),
        downloadButton("download_predAcc_info", "Download")
      ),
      pre(
        paste(capture.output(predAcc_summary_info()), collapse = "\n")
      )
    ))
  })


  # Download Summary Info
  output$download_predAcc_info <- downloadHandler(
    filename = function() {
      paste("predAcc_summary_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      # Write the summary info to a file
      writeLines(paste(capture.output(predAcc_summary_info()), collapse = "\n"), file)
    }
  )
}

## To be copied in the UI
# mod_GSAcc_ui("GSAcc_1")

## To be copied in the server
# mod_GSAcc_server("GSAcc_1")
