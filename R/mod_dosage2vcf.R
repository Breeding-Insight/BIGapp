#' dosage2vcf UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom shinyjs enable disable useShinyjs
#' @import shinydisconnect
#'
mod_dosage2vcf_ui <- function(id){
  ns <- NS(id)
  tagList(
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
        column(width = 5,
          box(
            title = "Inputs", width=12, status = "info", solidHeader = TRUE, collapsible = FALSE, collapsed = FALSE,
            selectInput(ns('file_type'), label = 'Select File Format', choices = c("DArT Dosage Reports","AgriSeq")),
            conditionalPanel(condition = "input.file_type == 'DArT Dosage Reports'",
                             ns = ns,
                             fileInput(ns("report_file"), "Choose DArT Dose Report File", accept = c(".csv")),
                             fileInput(ns("counts_file"), "Choose DArT Counts File", accept = c(".csv")),
            ),
            conditionalPanel(condition = "input.file_type == 'AgriSeq'",
                             ns = ns,
                             "Support for this file type is in-progress",
                             "",
                             #fileInput(ns("agriseq_file"), "Choose Input File", accept = c(".csv"))
            ),
            textInput(ns("d2v_output_name"), "Output File Name"),
            numericInput(ns("dosage2vcf_ploidy"), "Species Ploidy", min = 1, value = NULL),
            #actionButton(ns("run_analysis"), "Run Analysis"),
            useShinyjs(),
            downloadButton(ns('download_d2vcf'), "Download VCF File", class = "butt"),
            div(style="display:inline-block; float:right",dropdownButton(
              HTML("<b>Input files</b>"),
              p(downloadButton(ns('download_dose'), ""), "Dose Report Example File"),
              p(downloadButton(ns('download_counts'), ""), "Counts Example File"), hr(),
              p(HTML("<b>Parameters description:</b>"), actionButton(ns("goPar"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
              p(HTML("<b>Graphics description:</b>"), actionButton(ns("goRes"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
              p(HTML("<b>How to cite:</b>"), actionButton(ns("goCite"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
              actionButton(ns("d2vcf_summary"), "Summary"),
              circle = FALSE,
              status = "warning",
              icon = icon("info"), width = "500px",
              tooltip = tooltipOptions(title = "Click to see info!")
            ))
          )
      ),
      column(width = 4,  
        valueBoxOutput(ns("ReportSnps"), width=12),
        box(title = "Status", width = 12, collapsible = TRUE, status = "info",
            progressBar(id = ns("dosage2vcf_pb"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
          )
          ),
      column(width = 1),
        )
      )
    )
}

#' dosage2vcf Server Functions
#'
#' @import BIGr
#' @importFrom shinyjs enable disable
#'
#' @noRd
mod_dosage2vcf_server <- function(input, output, session, parent_session){

  ns <- session$ns

  # Help links
  observeEvent(input$goPar, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "DArT_Report2VCF_tabset",
                      selected = "DArT_Report2VCF_par")
    # expand specific box
    updateBox(id = "DArT_Report2VCF_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goRes, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "DArT_Report2VCF_tabset",
                      selected = "DArT_Report2VCF_results")
    # expand specific box
    updateBox(id = "DArT_Report2VCF_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goCite, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "DArT_Report2VCF_tabset",
                      selected = "DArT_Report2VCF_cite")
    # expand specific box
    updateBox(id = "DArT_Report2VCF_box", action = "toggle", session = parent_session)
  })

  snp_number <- reactiveVal(0)
  disable("download_d2vcf")

  #SNP counts value box
  output$ReportSnps <- renderValueBox({
    valueBox(snp_number(), "Number of Markers", icon = icon("dna"), color = "info")
  })

  observeEvent(input$run_analysis, {
    # Missing input with red border and alerts
    toggleClass(id = "d2v_output_name", class = "borderred", condition = (is.na(input$d2v_output_name) | is.null(input$d2v_output_name) | input$d2v_output_name == ""))
    toggleClass(id = "dosage2vcf_ploidy", class = "borderred", condition = (is.na(input$dosage2vcf_ploidy) | is.null(input$dosage2vcf_ploidy) | input$dosage2vcf_ploidy == ""))

    if (is.null(input$report_file$datapath) | is.null(input$counts_file$datapath)) {
      shinyalert(
        title = "Missing input!",
        text = "Upload Dose Report and Counts Files",
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
    req(input$report_file, input$counts_file, input$d2v_output_name, input$dosage2vcf_ploidy)

    updateProgressBar(session = session, id = "dosage2vcf_pb", value = 10, title = "Input files evaluated.")

    dosage_file_df <- read.csv(input$report_file$datapath)
    snp_number <- length(dosage_file_df$X.[-c(1:7)])

    #SNP counts value box
    output$ReportSnps <- renderValueBox({
      valueBox(snp_number, "Number of Markers", icon = icon("dna"), color = "info")
    })

    enable("download_d2vcf")

    updateProgressBar(session = session, id = "dosage2vcf_pb", value = 100, title = "Click in Download to continue.")

  })

  output$download_dose <- downloadHandler(
    filename = function() {
      paste0("BIGapp_Dose_Report_Example_file.csv")
    },
    content = function(file) {
      ex <- system.file("iris_DArT_Allele_Dose_Report.csv", package = "BIGapp")
      file.copy(ex, file)
    })

  output$download_counts <- downloadHandler(
    filename = function() {
      paste0("BIGapp_Counts_Example_file.csv")
    },
    content = function(file) {
      ex <- system.file("iris_DArT_Counts.csv", package = "BIGapp")
      file.copy(ex, file)
    })

  ##This is for the DArT files conversion to VCF
  output$download_d2vcf <- downloadHandler(
    filename = function() {
      output_name <- gsub("\\.vcf$", "", input$d2v_output_name)
      paste0(output_name, ".vcf.gz")
    },
    content = function(file) {
      # Ensure the files are uploaded
      req(input$report_file, input$counts_file, input$d2v_output_name, input$dosage2vcf_ploidy)

      # Get the uploaded file paths
      dosage_file <- input$report_file$datapath
      counts_file <- input$counts_file$datapath
      ploidy <- input$dosage2vcf_ploidy

      # Use a temporary file path without appending .vcf
      temp_base <- tempfile()

      #Status
      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 10, title = "Converting DArT files to VCF")

      # Convert to VCF using the BIGr package
      cat("Running BIGr::dosage2vcf...\n")
      dosage2vcf(
        dart.report = dosage_file,
        dart.counts = counts_file,
        output.file = temp_base,
        ploidy = as.numeric(ploidy)
      )

      # The output file should be temp_base.vcf
      output_name <- paste0(temp_base, ".vcf")

      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 50, title = "Writting vcf.")

      # Check if the VCF file was created
      if (file.exists(output_name)) {
        cat("VCF file created successfully.\n")

        # Compress the VCF file using gzip
        gzip_file <- paste0(output_name, ".gz")
        gz <- gzfile(gzip_file, "w")
        writeLines(readLines(output_name), gz)
        close(gz)

        # Check if the gzip file was created
        if (file.exists(gzip_file)) {
          cat("Gzip file created successfully.\n")

          # Move the compressed file to the path specified by 'file'
          file.copy(gzip_file, file)

          # Delete the temporary files
          unlink(gzip_file)
          unlink(output_name)

          cat("Temporary files deleted successfully.\n")
        } else {
          stop("Error: Failed to create the gzip file.")
        }
      } else {
        stop("Error: Failed to create the VCF file.")
      }

      #Status
      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 100, title = "Complete! - Downloading VCF")
    }
  )
  
  ##Summary Info
  d2vcf_summary_info <- function() {
    #Handle possible NULL values for inputs
    report_file_name <- if (!is.null(input$report_file$name)) input$report_file$name else "No file selected"
    counts_file_name <- if (!is.null(input$counts_file$name)) input$counts_file$name else "No file selected"
    selected_ploidy <- if (!is.null(input$dosage2vcf_ploidy)) as.character(input$dosage2vcf_ploidy) else "Not selected"
    
    #Print the summary information
    cat(
      "BIGapp Dosage2VCF Summary\n",
      "\n",
      paste0("Date: ", Sys.Date()), "\n",
      paste("R Version:", R.Version()$version.string), "\n",
      "\n",
      "### Input Files ###\n",
      "\n",
      paste("Input Dosage Report File:", report_file_name), "\n",
      paste("Input Counts File:", counts_file_name), "\n",
      "\n",
      "### User Selected Parameters ###\n",
      "\n",
      paste("Selected Ploidy:", selected_ploidy), "\n",
      "\n",
      "### R Packages Used ###\n",
      "\n",
      paste("BIGapp:", packageVersion("BIGapp")), "\n",
      paste("BIGr:", packageVersion("BIGr")), "\n",
      sep = ""
    )
  }
  
  # Popup for analysis summary
  observeEvent(input$d2vcf_summary, {
    showModal(modalDialog(
      title = "Summary Information",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close"),
        downloadButton("download_d2vcf_info", "Download")
      ),
      pre(
        paste(capture.output(d2vcf_summary_info()), collapse = "\n")
      )
    ))
  })
  
  
  # Download Summary Info
  output$download_d2vcf_info <- downloadHandler(
    filename = function() {
      paste("Dosage2VCF_summary_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      # Write the summary info to a file
      writeLines(paste(capture.output(d2vcf_summary_info()), collapse = "\n"), file)
    }
  )
}

## To be copied in the UI
# mod_dosage2vcf_ui("dosage2vcf_1")

## To be copied in the server
# mod_dosage2vcf_server("dosage2vcf_1")
