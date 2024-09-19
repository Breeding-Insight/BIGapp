#' DosageCall UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shinyjs enable disable useShinyjs
#' @importFrom shiny NS tagList
#' @importFrom future availableCores
#' @importFrom bs4Dash renderValueBox
#'
#'
mod_DosageCall_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        box(
          title = "Inputs", status = "info", solidHeader = TRUE, collapsible = FALSE, collapsed = FALSE,
          fileInput(ns("madc_file"), "Choose MADC or VCF File", accept = c(".csv",".vcf",".gz")),
          fileInput(ns("madc_passport"), "Choose Passport File (optional)", accept = c(".csv")),
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
          textInput(ns("output_name"), "Output File Name"),
          numericInput(ns("ploidy"), "Species Ploidy", min = 1, value = NULL),
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
          numericInput(ns("cores"), "Number of CPU Cores", min = 1, max = (future::availableCores() - 1), value = 1),
          actionButton(ns("run_analysis"), "Run Analysis"),
          useShinyjs(),
          downloadButton(ns('download_updog_vcf'), "Download VCF File", class = "butt"),

          div(style="display:inline-block; float:right",dropdownButton(
            HTML("<b>Input files</b>"),
            p(downloadButton(ns('download_vcf'),""), "VCF Example File"),
            p(downloadButton(ns('download_madc'),""), "MADC Example File"), hr(),
            p(HTML("<b>Parameters description:</b>"), actionButton(ns("goPar"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
            p(HTML("<b>Results description:</b>"), actionButton(ns("goRes"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
            p(HTML("<b>How to cite:</b>"), actionButton(ns("goCite"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
            p(HTML("<b>Updog tutorial:</b>"), actionButton(ns("goUpdog"), icon("arrow-up-right-from-square", verify_fa = FALSE), onclick ="window.open('https://dcgerard.github.io/updog/', '_blank')" )),
            circle = FALSE,
            status = "warning",
            icon = icon("info"), width = "500px",
            tooltip = tooltipOptions(title = "Click to see info!")
          ))
        ),
        valueBoxOutput(ns("MADCsnps"))
      ),

      fluidRow(
        box(title = "Status", width = 3, collapsible = TRUE, status = "info",
            progressBar(id = ns("pb_madc"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
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
#' @importFrom shinyjs enable disable
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
      need(!is.null(input$madc_passport), "Upload passport file to access results in this section."),
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

  disable("download_updog_vcf")

  ##This is for performing Updog Dosage Calling
  updog_out <- eventReactive(input$run_analysis,{

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
  })

  # Only make available the download button when analysis is finished
  observe({
    if (!is.null(updog_out())) {
      Sys.sleep(1)
      # enable the download button
      enable("download_updog_vcf")
    } else {
      disable("download_updog_vcf")
    }
  })

  output$download_updog_vcf <- downloadHandler(
    filename = function() {
      paste0(input$output_name, ".vcf.gz")
    },
    content = function(file) {
      #Save Updog output as VCF file
      temp <- tempfile()
      updog2vcf(
        multidog.object = updog_out(),
        output.file = temp,
        updog_version = packageVersion("updog"),
        compress = TRUE
      )

      # Move the file to the path specified by 'file'
      file.copy(paste0(temp, ".vcf.gz"), file)

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
}

## To be copied in the UI
# mod_DosageCall_ui("DosageCall_1")

## To be copied in the server
# mod_DosageCall_server("DosageCall_1")
