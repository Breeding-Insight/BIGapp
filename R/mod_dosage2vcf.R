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
#'
mod_dosage2vcf_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        box(
          title = "Inputs", status = "info", solidHeader = TRUE, collapsible = FALSE, collapsed = FALSE,
          fileInput(ns("report_file"), "Choose DArT Dose Report File", accept = c(".csv")),
          fileInput(ns("counts_file"), "Choose DArT Counts File", accept = c(".csv")),
          textInput(ns("d2v_output_name"), "Output File Name"),
          numericInput(ns("dosage2vcf_ploidy"), "Species Ploidy", min = 1, value = NULL),
          actionButton(ns("run_analysis"), "Run Analysis"),
          useShinyjs(),
          downloadButton(ns('download_d2vcf'), "Download VCF File", class = "butt"),
          div(style="display:inline-block; float:right",dropdownButton(

            tags$h3("DArT File Conversion"),
            "Converting DArT report files to VCF format. \n",
            "You can download examples of the expected files here: \n",
            downloadButton(ns('download_dose'), "Download Dose Report Example File"),
            downloadButton(ns('download_counts'), "Download Counts Example File"),
            circle = FALSE,
            status = "warning",
            icon = icon("info"), width = "500px",
            tooltip = tooltipOptions(title = "Click to see info!")
          ))
        ),
        valueBoxOutput(ns("ReportSnps"))
      ),
      fluidRow(
        box(title = "Status", width = 3, collapsible = TRUE, status = "info",
            progressBar(id = ns("dosage2vcf_pb"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
        )
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

    dosage_file_df <- read.csv(input$report_file$datapath)
    snp_number <- length(dosage_file_df$X.[-c(1:7)])

    #SNP counts value box
    output$ReportSnps <- renderValueBox({
      valueBox(snp_number, "Number of Markers", icon = icon("dna"), color = "info")
    })

    enable("download_d2vcf")
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
      paste0(input$d2v_output_name, ".vcf.gz")
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
      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 50, title = "Converting DArT files to VCF")

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
}

## To be copied in the UI
# mod_dosage2vcf_ui("dosage2vcf_1")

## To be copied in the server
# mod_dosage2vcf_server("dosage2vcf_1")
