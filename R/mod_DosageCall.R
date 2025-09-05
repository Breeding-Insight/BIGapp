#' DosageCall UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom future availableCores
#' @importFrom bs4Dash renderValueBox
#' @import shinydisconnect
#'
#'
mod_DosageCall_ui <- function(id){
  ns <- NS(id)
  tagList(
    useShinyjs(),
    fluidPage(
      disconnectMessage(
        text = "An unexpected error occurred, please reload the application and check the input file(s).",
        refresh = "Reload now",
        background = "white",
        colour = "grey",
        overlayColour = "grey",
        overlayOpacity = 0.3,
        refreshColour = "purple"
      ),
      fluidRow(
        box(
          title = "Inputs", status = "info", solidHeader = TRUE, collapsible = FALSE, collapsed = FALSE,
          "* Required",
          selectInput(ns("Rpackage"), "Dosage Calling Method", choices = c("Updog", "polyRAD"), selected = "Updog"),
          fileInput(ns("madc_file"), "Choose VCF File*", accept = c(".csv",".vcf",".gz")),
          fileInput(ns("madc_passport"), "Choose Trait File", accept = c(".csv")),
          conditionalPanel(
            condition = "output.passportTablePopulated",
            ns = ns,
            tags$div(
              style = "padding-left: 20px;",  # Add padding/indentation
              virtualSelectInput(
                inputId = ns("cat_madc"),
                label = "Select Category Subset:",
                choices = NULL,
                showValueAsTags = TRUE,
                search = TRUE,
                multiple = FALSE
              ),
              virtualSelectInput(
                inputId = ns("item_madc"),
                label = "Select Subset Values:",
                choices = NULL,
                showValueAsTags = TRUE,
                search = TRUE,
                multiple = TRUE
              )
            )
          ),
          conditionalPanel(
            condition = "input.Rpackage == 'Updog'",
            ns = ns,
            textInput(ns("output_name"), "Output File Name*"),
            numericInput(ns("ploidy"), "Species Ploidy*", min = 1, value = NULL),
            selectInput(ns("updog_model"), "Updog Model", choices = c("norm","hw","bb","s1","s1pp","f1","f1pp","flex","uniform"), selected = "norm"),
            conditionalPanel(
              condition = "input.updog_model == 'f1' | input.updog_model == 'f1pp'",
              ns = ns,
              tags$div(
                style = "padding-left: 20px;",  # Add padding/indentation
                textInput(
                  inputId = ns("parent1"),
                  label = "Enter parent1 ID:",
                  value = NULL
                ),
                textInput(
                  inputId = ns("parent2"),
                  label = "Enter parent2 ID:",
                  value = NULL
                )
              )
            ),
            conditionalPanel(
              condition = "input.updog_model == 's1' | input.updog_model == 's1pp'",
              ns = ns,
              tags$div(
                style = "padding-left: 20px;",  # Add padding/indentation
                textInput(
                  inputId = ns("parent"),
                  label = "Enter parent ID:",
                  value = NULL
                )
              )
            ),
            sliderInput(ns("cores"), "Number of CPU Cores*", min = 1, max = (availableCores() - 1), value = 1, step = 1),
          ),
          conditionalPanel(
            condition = "input.Rpackage == 'polyRAD'",
            ns = ns,
            textInput(ns("output_name"), "Output File Name*"),
            numericInput(ns("ploidy"), "Species Ploidy*", min = 1, value = NULL),
            #selectInput(ns("polyRAD_model"), "polyRAD Model", choices = c("IterateHWE","Iterate_LD","IteratePopStruct","IteratePopStruct_LD","PipelineMapping2Parents"), selected = "IterateHWE"),
            selectInput(ns("polyRAD_model"), "polyRAD Model", choices = c("IterateHWE","IteratePopStruct","PipelineMapping2Parents"), selected = "IterateHWE"),
            conditionalPanel(
              condition = "input.polyRAD_model == 'PipelineMapping2Parents'",
              ns = ns,
              tags$div(
                style = "padding-left: 20px;",  # Add padding/indentation
                textInput(
                  inputId = ns("parent1"),
                  label = "Donor Parent ID:",
                  value = NULL
                ),
                textInput(
                  inputId = ns("parent2"),
                  label = "Recurrent Parent ID:",
                  value = NULL
                ),
                numericInput(
                  inputId = ns("bx.gen"),
                  label = "No. of Backcross Generations:",
                  value = 0
                ),
                numericInput(
                  inputId = ns("inter.gen"),
                  label = "No. of Intermating Generations:",
                  value = 0
                ),
                numericInput(
                  inputId = ns("self.gen"),
                  label = "No. of Selfing Geneerations:",
                  value = 0
                )
              )
            )
          ),
          actionButton(ns("run_analysis"), "Run Analysis"),
          uiOutput(ns('mybutton')),

          div(style="display:inline-block; float:right",dropdownButton(
            HTML("<b>Input files</b>"),
            p(downloadButton(ns('download_vcf'),""), "VCF Example File"),hr(),
            p(HTML("<b>Parameters description:</b>"), actionButton(ns("goPar"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
            p(HTML("<b>Results description:</b>"), actionButton(ns("goRes"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
            p(HTML("<b>How to cite:</b>"), actionButton(ns("goCite"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
            p(HTML("<b>Updog tutorial:</b>"), actionButton(ns("goUpdog"), icon("arrow-up-right-from-square", verify_fa = FALSE), onclick ="window.open('https://dcgerard.github.io/updog/', '_blank')" )), hr(),
            actionButton(ns("dosage_summary"), "Summary"),
            circle = FALSE,
            status = "warning",
            icon = icon("info"), width = "500px",
            tooltip = tooltipOptions(title = "Click to see info!")
          )),
          br(),
          tags$hr(style="border-color: #d3d3d3; margin-top: 20px; margin-bottom: 20px;"),  # Lighter grey line
          div(style="text-align: left; margin-top: 10px;",
              actionButton(ns("advanced_options"),
                           label = HTML(paste(icon("cog", style = "color: #007bff;"), "Advanced Options")),
                           style = "background-color: transparent; border: none; color: #007bff; font-size: smaller; text-decoration: underline; padding: 0;"
              )
          )
        ),
        column(width=4,
               valueBoxOutput(ns("MADCsnps"), width=12),
               box(title = "Status", width = 12, collapsible = TRUE, status = "info",
                   progressBar(id = ns("pb_madc"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
               )
        )
      )
    )
  )
}

#' DosageCall Server Functions
#'
#' @import vcfR
#' @import updog
#' @importFrom BIGr updog2vcf
#' @import dplyr
#'
#' @noRd
mod_DosageCall_server <- function(input, output, session, parent_session){

  ns <- session$ns

  # Help links
  observeEvent(input$goPar, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "Updog_Dosage_Calling_tabset",
                      selected = "Updog_Dosage_Calling_par")
    # expand specific box
    updateBox(id = "Updog_Dosage_Calling_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goRes, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "Updog_Dosage_Calling_tabset",
                      selected = "Updog_Dosage_Calling_results")
    # expand specific box
    updateBox(id = "Updog_Dosage_Calling_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goCite, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "Updog_Dosage_Calling_tabset",
                      selected = "Updog_Dosage_Calling_cite")
    # expand specific box
    updateBox(id = "Updog_Dosage_Calling_box", action = "toggle", session = parent_session)
  })

  # Update dropdown menu choices based on uploaded passport file
  passport_table <- reactive({
    validate(
      need(!is.null(input$madc_passport), "Upload Trait File to access results in this section."),
    )
    info_df <- read.csv(input$madc_passport$datapath, header = TRUE, check.names = FALSE)
    info_df[,1] <- as.character(info_df[,1]) #Makes sure that the sample names are characters instead of numeric

    updateVirtualSelect("cat_madc", choices = colnames(info_df), session = session)
    info_df
  })

  # Server logic to check if passport_table() has data
  output$passportTablePopulated <- reactive({
    !is.null(passport_table()) && nrow(passport_table()) > 0  # Check if the table has rows
  })
  outputOptions(output, "passportTablePopulated", suspendWhenHidden = FALSE)

  #MADC specific category selection
  observeEvent(input$cat_madc, {

    # Get selected column name
    selected_col <- input$cat_madc

    # Extract unique values from the selected column
    unique_values <- unique(passport_table()[[selected_col]])

    #Add category selection
    updateVirtualSelect("item_madc", choices = unique_values, session = session)

  })

  snp_number <- reactiveVal(0)

  #SNP counts value box
  output$MADCsnps <- renderValueBox({
    valueBox(snp_number(), "Markers in uploaded file", icon = icon("dna"), color = "info")
  })

  #Default model choices
  advanced_options <- reactiveValues(
    contamRate = 0.001,
    min_ind_read = 1,
    min_ind_maf = 0,
    tol = 1e-05
  )

  #UI popup window for input
  observeEvent(input$advanced_options, {
    showModal(modalDialog(
      title = "Advanced Options",
      conditionalPanel(
        condition = "input.Rpackage == 'polyRAD'", ns = ns,
        div(
          numericInput(
            inputId = ns('contamRate'),
            label = 'contamRate',
            min = 0,
            max = 1,
            value = advanced_options$contamRate
          ),
          numericInput(
            inputId = ns('min_ind_read'),
            label = 'min.ind.with.reads',
            min = 1,
            value = advanced_options$min_ind_read
          ),
          numericInput(
            inputId = ns('min_ind_maf'),
            label = 'min.ind.with.minor.allele',
            min = 0,
            value = advanced_options$min_ind_maf
          ),
          numericInput(
            inputId = ns('tol'),
            label = 'tol',
            min = 0,
            max = 1,
            value = advanced_options$tol
          )
        )
      ),
      conditionalPanel(
        condition = "input.Rpackage == 'Updog'", ns = ns,
        "Currently, there are no Advanced Options for Updog"
      ),
      footer = tagList(
        modalButton("Close"),
        actionButton(ns("save_advanced_options"), "Save")
      )
    ))
  })



  #Close popup window when user "saves options"
  observeEvent(input$save_advanced_options, {
    advanced_options$contamRate <- input$contamRate
    advanced_options$min_ind_read <- input$min_ind_read
    advanced_options$min_ind_maf <- input$min_ind_maf
    advanced_options$tol <- input$tol
    # Save other inputs as needed

    removeModal()  # Close the modal after saving
  })

  disable("download_updog_vcf")
  
  #Default model choices
  polyrad_items <- reactiveValues(
    vcf_path = NULL,
  )

  ##This is for performing Dosage Calling

  updog_out <- eventReactive(input$run_analysis,{
    # Updog
    if (input$Rpackage == "Updog") {
      # Missing input with red border and alerts
      toggleClass(id = "ploidy", class = "borderred", condition = (is.na(input$ploidy) | is.null(input$ploidy)))
      toggleClass(id = "output_name", class = "borderred", condition = (is.na(input$output_name) | is.null(input$output_name) | input$output_name == ""))
      toggleClass(id = "parent", class = "borderred", condition = ((input$updog_model == "s1" | input$updog_model == "s1pp") & (is.null(input$parent) | input$parent == "")))
      toggleClass(id = "parent1", class = "borderred", condition = ((input$updog_model == "f1" | input$updog_model == "f1pp") & (is.null(input$parent1) | input$parent1 == "")))
      toggleClass(id = "parent2", class = "borderred", condition = ((input$updog_model == "f1" | input$updog_model == "f1pp") & (is.null(input$parent2) | input$parent2 == "")))

      if (is.null(input$madc_file$datapath)) {
        shinyalert(
          title = "Missing input!",
          text = "Upload VCF File",
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
      req(input$madc_file$datapath, input$output_name, input$ploidy)

      # Get inputs
      madc_file <- input$madc_file$datapath
      output_name <- input$output_name
      ploidy <- input$ploidy
      cores <- input$cores

      # Status
      updateProgressBar(session = session, id = "pb_madc", value = 0, title = "Formatting Input Files")
      #Import genotype info if genotype matrix format
      if (grepl("\\.csv$", madc_file)) {
        # Call the get_counts function with the specified MADC file path and output file path
        #Status
        result_df <- get_counts(madc_file, output_name)

        #Call the get_matrices function
        matrices <- get_matrices(result_df)

        #Number of SNPs
        snp_number <- (nrow(result_df) / 2)

        #SNP counts value box
        output$MADCsnps <- renderValueBox({
          valueBox(snp_number, "Markers in MADC File", icon = icon("dna"), color = "info")
        })

      } else {

        #Initialize matrices list
        matrices <- list()

        #Import genotype information if in VCF format
        #### VCF sanity check
        checks <- vcf_sanity_check(madc_file, depth_support_fields = c("AD","RA"), max_markers = 10000)
        
        error_if_false <- c(
          "VCF_header", "VCF_columns", "unique_FORMAT", 
          "samples", "chrom_info", "pos_info", "allele_counts", "VCF_compressed"
        )
        
        error_if_true <- c(
          "duplicated_samples", "duplicated_markers"
        )
        warning_if_false <- c("ref_alt","max_markers")
        
        checks_result <- vcf_sanity_messages(checks, 
                                             error_if_false, 
                                             error_if_true, 
                                             warning_if_false = warning_if_false, 
                                             warning_if_true = NULL,
                                             input_ploidy = NULL)
        
        if(checks_result) return() # Stop the analysis if checks fail
        #########
        
        vcf <- read.vcfR(madc_file, verbose = FALSE)

        #Get items in FORMAT column
        info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT
        chrom <- vcf@fix[,1]
        pos <- vcf@fix[,2]

        info_ids <- extract_info_ids(info[1])

        if (("DP" %in% info_ids) && (("RA" %in% info_ids) | ("AD" %in% info_ids))) {
          #Extract DP and RA and convert to matrices
          matrices$size_matrix <- extract.gt(vcf, element = "DP")
          if("RA" %in% info_ids){
            matrices$ref_matrix <- extract.gt(vcf, element = "RA")
          } else {
            ad_matrix <- extract.gt(vcf, element = "AD")
            matrices$ref_matrix <- matrix(sapply(strsplit(ad_matrix, ","), "[[", 1), nrow = nrow(matrices$size_matrix))
            colnames(matrices$ref_matrix) <- colnames(matrices$size_matrix)
            rownames(matrices$ref_matrix) <- rownames(matrices$size_matrix)
          }

          class(matrices$size_matrix) <- "numeric"
          class(matrices$ref_matrix) <- "numeric"
          rownames(matrices$size_matrix) <- rownames(matrices$ref_matrix) <- paste0(chrom, "_", pos)

          rm(vcf) #Remove VCF

          snp_number <- (nrow(matrices$size_matrix))

          #SNP counts value box
          output$MADCsnps <- renderValueBox({
            valueBox(snp_number, "Markers in VCF File", icon = icon("dna"), color = "info")
          })

        }else{
          ##Add user warning about read depth and allele read depth not found
          stop(safeError("Error: DP and RA/AD FORMAT flags not found in VCF file"))
        }
      }

      #Subset samples from the matrices if the user selected items in the passport file
      if (!is.null(input$item_madc) && length(input$item_madc) > 0){

        #First getting the samples that are both in the passport and the MADC/VCF file
        #**Assuming the first column of the passport table is the sample IDs
        shared_samples <- intersect(passport_table()[[1]], colnames(matrices$ref_matrix))

        # Filter the passport dataframe
        filtered_shared_samples <- passport_table() %>%
          filter(passport_table()[[1]] %in% shared_samples,
                 passport_table()[[input$cat_madc]] %in% input$item_madc) %>%
          pull(1)

        #Give warning if no samples were subset
        if (length(filtered_shared_samples) < 1) {
          shinyalert(
            title = "Data Warning!",
            text = "No samples remain after subsetting options",
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

        #Subset the matrices
        matrices$ref_matrix <- matrices$ref_matrix[, filtered_shared_samples]
        matrices$size_matrix <- matrices$size_matrix[, filtered_shared_samples]

      }

      # Select parents
      if(input$updog_model == "s1" | input$updog_model == "s1pp"){
        parents <- c(input$parent, NULL)
      } else if(input$updog_model == "f1" | input$updog_model == "f1pp"){
        parents <- c(input$parent1, input$parent2)
      } else {
        parents <- c(NULL, NULL)
      }

      if (!all(parents %in% colnames(matrices$size_matrix))) {
        shinyalert(
          title = "Data Warning!",
          text = "Parent(s) not found. Check the genotype input file parent(s) ID, make sure they match with the input parent(s) ID.",
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

      if (nrow(matrices$ref_matrix) == 0 || nrow(matrices$size_matrix) == 0) {
        shinyalert(
          title = "Data Warning!",
          text = "All markers are missing read count information for reference and alternate alleles",
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

      #Run Updog
      #I am also taking the ploidy from the max value in the
      updateProgressBar(session = session, id = "pb_madc", value = 40, title = "Dosage Calling in Progress")
      mout <- multidog(refmat = matrices$ref_matrix,
                       sizemat = matrices$size_matrix,
                       ploidy = as.numeric(ploidy),
                       p1_id = parents[1],
                       p2_id = parents[2],
                       model = input$updog_model,
                       nc = cores)
      #Status
      updateProgressBar(session = session, id = "pb_madc", value = 100, title = "Finished")
      mout
    } else {
      # PolyRAD
      
      #Variables
      polyrad_items$vcf_path <- input$madc_file$datapath
      
      #Status
      updateProgressBar(session = session, id = "pb_madc", value = 15, title = "Formatting Input File")
      
      #Subset samples from the VCF if the user selected items in the passport file
      if (!is.null(input$item_madc) && length(input$item_madc) > 0){
        
        #Status
        updateProgressBar(session = session, id = "pb_madc", value = 25, title = "Subsetting Input File")
        
        #First getting the samples that are both in the passport and the MADC/VCF file
        #**Assuming the first column of the passport table is the sample IDs
        #*#populate preview_data
        preview_vcf <- read.vcfR(input$madc_file$datapath, verbose = FALSE, nrows = 1)
        
        #Get names
        sample_list <- names(data.frame(preview_vcf@gt, check.names=FALSE)[,-1])
        
        shared_samples <- intersect(passport_table()[[1]], sample_list)
        
        # Filter the passport dataframe (samples to kee[])
        filtered_shared_samples <- passport_table() %>%
          filter(passport_table()[[1]] %in% shared_samples,
                 passport_table()[[input$cat_madc]] %in% input$item_madc) %>%
          pull(1)
        
        # Find samples to remove (samples in the vcf that are not in the filtered passport table)
        filtered_shared_samples_remove <- setdiff(sample_list, filtered_shared_samples)
        
        #Give warning if no samples were subset
        if (length(filtered_shared_samples) < 1) {
          shinyalert(
            title = "Data Warning!",
            text = "No samples remain after subsetting options",
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
        
        #Subset the VCF
        temp_vcf <- read.vcfR(polyrad_items$vcf_path, verbose = FALSE)
        temp_vcf <- subset_vcf(temp_vcf, remove.sample.list = filtered_shared_samples_remove)[[1]]
        #Write to temp location
        temp_file_subset <- tempfile(fileext = ".vcf.gz")
        write.vcf(temp_vcf, file = temp_file_subset)
        
        #Update path
        polyrad_items$vcf_path <- temp_file_subset
        
        #reduce memory 
        rm(temp_vcf)
      }
      
      updateProgressBar(session = session, id = "pb_madc", value = 35, title = "Performing Dosage Calling")
      polyrad_results <- polyRAD_dosage_call(vcf = polyrad_items$vcf_path,
                                           ploidy = input$ploidy,
                                           model = input$polyRAD_model,
                                           p1 = input$parent1,
                                           p2 = input$parent2,
                                           backcross.gen = input$bx.gen,
                                           intercross.gen = input$inter.gen,
                                           selfing.gen = input$self.gen,
                                           contamRate = advanced_options$contamRate,
                                           min_ind_read = advanced_options$min_ind_read,
                                           min_ind_maf = advanced_options$min_ind_maf,
                                           tol=advanced_options$tol)
      updateProgressBar(session = session, id = "pb_madc", value = 100, title = "Finished")
      polyrad_results$vcf_path <- polyrad_items$vcf_path
      polyrad_results
    }
  })

  # Only make available the download button when analysis is finished
  output$mybutton <- renderUI({
    if(isTruthy(updog_out()))
      downloadButton(ns("download_updog_vcf"), "Download VCF file", class = "butt")
  })

  output$download_updog_vcf <- downloadHandler(
    filename = function() {
      output_name <- gsub("\\.vcf$", "", input$output_name)
      paste0(output_name, ".vcf.gz")
    },
    content = function(file) {

      #Save Updog output as VCF file
      temp <- tempfile()
      if (input$Rpackage == "Updog") {
        updog2vcf(
          multidog.object = updog_out(),
          output.file = temp,
          updog_version = packageVersion("updog"),
          compress = FALSE
        )

      } else {
        polyRAD2vcf(updog_out()$Genos,
                    model = input$polyRAD_model,
                    vcf_path = polyrad_items$vcf_path,
                    hindhe.obj = updog_out()$RADHindHe,
                    ploidy = input$ploidy,
                    output.file = temp
        )
      }
      bgzip_compress(paste0(temp, ".vcf"), file)

      # Delete the temporary file
      unlink(temp)
    })

  output$download_vcf <- downloadHandler(
    filename = function() {
      paste0("BIGapp_VCF_Example_file.vcf.gz")
    },
    content = function(file) {
      ex <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")
      file.copy(ex, file)
    })

  output$download_madc <- downloadHandler(
    filename = function() {
      paste0("BIGapp_MADC_Example_file.csv")
    },
    content = function(file) {
      ex <- system.file("iris_DArT_MADC.csv", package = "BIGapp")
      file.copy(ex, file)
    })

  ##Summary Info
  dosage_summary_info <- function() {
    #Handle possible NULL values for inputs
    genotype_file_name <- if (!is.null(input$madc_file$name)) input$madc_file$name else "No file selected"
    report_file_name <- if (!is.null(input$madc_passport$name)) input$madc_passport$name else "No file selected"
    selected_ploidy <- if (!is.null(input$ploidy)) as.character(input$ploidy) else "Not selected"

    #Print the summary information
    cat(
      "BIGapp Dosage Calling Summary\n",
      "\n",
      paste0("Date: ", Sys.Date()), "\n",
      paste("R Version:", R.Version()$version.string), "\n",
      "\n",
      "### Input Files ###\n",
      "\n",
      paste("Input Genotype File:", genotype_file_name), "\n",
      paste("Input Passport File:", report_file_name), "\n",
      "\n",
      "### User Selected Parameters ###\n",
      "\n",
      paste("Selected Ploidy:", selected_ploidy), "\n",
      paste("Selected Updog Model:", input$updog_model), "\n",
      "\n",
      "### R Packages Used ###\n",
      "\n",
      paste("BIGapp:", packageVersion("BIGapp")), "\n",
      paste("BIGr:", packageVersion("BIGr")), "\n",
      paste("Updog:", packageVersion("updog")), "\n",
      paste("PolyRAD:", packageVersion("polyRAD")), "\n",
      paste("vcfR:", packageVersion("vcfR")), "\n",
      paste("dplyr:", packageVersion("dplyr")), "\n",
      sep = ""
    )
  }

  # Popup for analysis summary
  observeEvent(input$dosage_summary, {
    showModal(modalDialog(
      title = "Summary Information",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close"),
        downloadButton("download_dosage_info", "Download")
      ),
      pre(
        paste(capture.output(dosage_summary_info()), collapse = "\n")
      )
    ))
  })


  # Download Summary Info
  output$download_dosage_info <- downloadHandler(
    filename = function() {
      paste("DosageCalling_summary_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      # Write the summary info to a file
      writeLines(paste(capture.output(dosage_summary_info()), collapse = "\n"), file)
    }
  )
}

## To be copied in the UI
# mod_DosageCall_ui("DosageCall_1")

## To be copied in the server
# mod_DosageCall_server("DosageCall_1")
