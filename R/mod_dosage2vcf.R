#' dosage2vcf UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_dosage2vcf_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        box(
          title = "Inputs", status = "info", solidHeader = TRUE, collapsible = FALSE, collapsed = FALSE,
          fileInput(ns("report_file"), "Choose DArT Dose Report File", accept = c(".csv")),
          fileInput(ns("counts_file"), "Choose DArT Counts File", accept = c(".csv")),
          #checkboxInput("off-targets","Include off-target loci?"),
          #fileInput("sample_file", "Optional: Choose Sample List (disabled)", accept = c(".csv")),
          textInput(ns("d2v_output_name"), "Output File Name"),
          numericInput(ns("dosage2vcf_ploidy"), "Species Ploidy", min = 1, value = NULL),
          downloadButton(ns("download_d2vcf"), "Download VCF File"),
          div(style="display:inline-block; float:right",dropdownButton(

            tags$h3("DArT File Converstion"),
            "Converting DArT report files to VCF format. The VCF file will automatically
                    download when complete.",
            circle = FALSE,
            status = "warning",
            icon = icon("info"), width = "300px",
            tooltip = tooltipOptions(title = "Click to see info!")
          ))
        ),
        valueBoxOutput(ns("ReportSnps"))
        #valueBox("Help","Updog Manual", icon = icon("globe"), color = "warning")
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
#'
#' @noRd
mod_dosage2vcf_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    snp_number <- reactiveVal(0)

    #SNP counts value box
    output$ReportSnps <- renderValueBox({
      valueBox(snp_number(), "Markers in VCF File", icon = icon("dna"), color = "info")
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
  })
}

## To be copied in the UI
# mod_dosage2vcf_ui("dosage2vcf_1")

## To be copied in the server
# mod_dosage2vcf_server("dosage2vcf_1")
