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
                 "* Required",
                 selectInput(ns('file_type'),
                             label = 'Select File Format',
                             choices = c("DArT MADC file","DArT Dosage Reports", "Dosage Matrix"), selected = "DArT MADC file"),
                 conditionalPanel(condition = "input.file_type == 'DArT Dosage Reports'",
                                  ns = ns,
                                  fileInput(ns("report_file"), "Choose DArT Dose Report File", accept = c(".csv")),
                                  fileInput(ns("counts_file"), "Choose DArT Counts File", accept = c(".csv")),
                                  numericInput(ns("dosage2vcf_ploidy"), "Species Ploidy", min = 1, value = NULL)
                 ),
                 conditionalPanel(condition = "input.file_type == 'Dosage Matrix'",
                                  ns = ns,
                                  fileInput(ns("matrix_file"), "Choose Dosage Matrix File", accept = c(".csv")),
                                  selectInput(ns("dosage_counts"), "Dosage Allele Count", choices = c("Reference","Alternate"), selected = "Reference"),
                                  numericInput(ns("dosage2vcf_ploidy"), "Species Ploidy", min = 1, value = NULL)
                 ),
                 conditionalPanel(condition = "input.file_type == 'DArT MADC file'",
                                  ns = ns,
                                  fileInput(ns("madc_file"), "Choose DArT MADC File*", accept = c(".csv"), multiple = TRUE),
                                  radioButtons(ns("snp_type"),
                                               label = "Select Marker Type",
                                               choices = list("Target"= "target", "Target + Off-Target" = "target_off"),
                                               selected = "target"),
                                  conditionalPanel(condition = "input.snp_type == 'target_off'",
                                                   ns = ns,
                                                   fileInput(ns("hapDB_file"), "Upload haplotype database file (fasta)*"),
                                                   fileInput(ns("botloci_file"), "Upload bottom strand probes file (.botloci)*"),
                                                   sliderInput(ns("cores"), "Number of CPU Cores*", min = 1, max = (availableCores() - 1), value = 1, step = 1)
                                  ),
                                  conditionalPanel(condition = "input.snp_type == 'target'",
                                                   ns = ns,
                                                   radioButtons(ns("ref_alt"),
                                                                label = "Extract REF and ALT info:",
                                                                choices = list("Yes"= "yes", "No" = "no"),
                                                                selected = "no")                                  )
                 ),
                 textInput(ns("d2v_output_name"), "Output File Name"),
                 actionButton(ns("run_analysis"), "Convert File"),
                 uiOutput(ns('mybutton')),

                 div(style="display:inline-block; float:right",dropdownButton(
                   HTML("<b>Input files</b>"),
                   p(downloadButton(ns('download_dose'), ""), "Dose Report Example File"),
                   p(downloadButton(ns('download_counts'), ""), "Counts Example File"),
                   p(downloadButton(ns('download_madc'), ""), "MADC Example File"),hr(),
                   p(HTML("<b>Parameters description:</b>"), actionButton(ns("goPar"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>Results description:</b>"), actionButton(ns("goRes"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
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
#' @importFrom Rsamtools bgzip
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

  output$download_madc <- downloadHandler(
    filename = function() {
      paste0("BIGapp_MADC_Example_file.csv")
    },
    content = function(file) {
      ex <- system.file("iris_DArT_MADC.csv", package = "BIGapp")
      file.copy(ex, file)
    })


  vcf_out <- eventReactive(input$run_analysis,{
    # Ensure the files are uploaded
    # Missing input with red border and alerts
    if(input$file_type == "DArT Dosage Reports"){
      if (is.null(input$report_file$datapath) | is.null(input$counts_file$datapath) | input$d2v_output_name == "" | input$dosage2vcf_ploidy == "") {
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
      
      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 50, title = "Writing VCF")
      
      dosage2vcf(
        dart.report = dosage_file,
        dart.counts = counts_file,
        output.file = temp_base,
        ploidy = as.numeric(ploidy)
      )

      # The output file should be temp_base.vcf
      output_name <- paste0(temp_base, ".vcf")

      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 100, title = "Complete!")
      
      return(output_name)

    } else if(input$file_type == "DArT MADC file"){
      req(input$madc_file)
      # First check if the MADC file is valid (a non-fixedAlleleID MADC file)

      #Read only the first column for the first seven rows
      first_seven_rows <- read.csv(input$madc_file$datapath, header = FALSE, nrows = 7, colClasses = c(NA, "NULL"))

      #Check if all entries in the first column are either blank or "*"
      check_entries <- all(first_seven_rows[, 1] %in% c("", "*"))

      #Check if the MADC file has the filler rows or is processed from updated fixed allele ID pipeline
      if (check_entries) {
        #Note: This assumes that the first 7 rows are placeholder info from DArT processing
        shinyalert(
          title = "Raw MADC!",
          text = "This MADC file has not been processed by the updated Breeding Insight fixed allele ID pipeline.",
          size = "m",
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

        #Exit the analysis
        return()
      }

      #Now perform conversion depending on user options

      if(input$snp_type == "target_off"){
        if (is.null(input$madc_file$datapath) | input$hapDB_file$datapath == "" | input$botloci_file$datapath == "" | input$d2v_output_name == "") {
          shinyalert(
            title = "Missing input!",
            text = "Upload MADC, HaplotypeDB and BOTLOCI files",
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
        req(input$madc_file, input$hapDB_file, input$botloci_file)

        # Use a temporary file path without appending .vcf
        temp_base <- tempfile()

        # merge MADC if multiple
        if(length(input$madc_file$datapath) >1){
          updateProgressBar(session = session, id = "dosage2vcf_pb", value = 15, title = "Merging MADC files")

          merged_madc <- paste0(temp_base, ".csv")

          run_ids <- sapply(strsplit(input$madc_file$name, "_"), "[[",1)
          if(length(run_ids) == 0) run_ids <- NULL

          merge_MADCs(madc_list = as.list(input$madc_file$datapath), out_madc = merged_madc,run_ids = run_ids)
          read_madc <- merged_madc
        } else read_madc <- input$madc_file$datapath

        # The output file should be temp_base.vcf
        output_name <- paste0(temp_base, ".vcf")

        updateProgressBar(session = session, id = "dosage2vcf_pb", value = 30, title = "Writing VCF")

        get_OffTargets(madc = read_madc,
                       botloci = input$botloci_file$datapath,
                       hap_seq = input$hapDB_file$datapath,
                       n.cores= input$cores,
                       rm_multiallelic_SNP = TRUE,
                       out_vcf = output_name,
                       verbose = FALSE)

        updateProgressBar(session = session, id = "dosage2vcf_pb", value = 100, title = "Complete!")

        return(output_name)

      } else if(input$snp_type == "target"){

        # Use a temporary file path without appending .vcf
        temp_base <- tempfile()

        # merge MADC if multiple
        if(length(input$madc_file$datapath) >1){
          updateProgressBar(session = session, id = "dosage2vcf_pb", value = 15, title = "Merging MADC files")

          merged_madc <- paste0(temp_base, ".csv")

          run_ids <- sapply(strsplit(input$madc_file$name, "_"), "[[",1)
          if(length(run_ids) == 0) run_ids <- NULL

          merge_MADCs(madc_list = as.list(input$madc_file$datapath), out_madc = merged_madc,run_ids = run_ids)
          read_madc <- merged_madc
        } else read_madc <- input$madc_file$datapath


        # The output file should be temp_base.vcf
        output_name <- paste0(temp_base, ".vcf")

        updateProgressBar(session = session, id = "dosage2vcf_pb", value = 30, title = "Writing VCF")
        madc2vcf(read_madc, output_name, get_REF_ALT = input$ref_alt)

        updateProgressBar(session = session, id = "dosage2vcf_pb", value = 100, title = "Complete!")
        return(output_name)
      }
      
    } else if (input$file_type == "Dosage Matrix"){
      
      if (is.null(input$matrix_file$datapath) | input$d2v_output_name == "" | input$dosage2vcf_ploidy == "") {
        shinyalert(
          title = "Missing input!",
          text = "Upload Dosage Matrix",
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
      req(input$matrix_file, input$d2v_output_name, input$dosage2vcf_ploidy)
      #Status
      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 10, title = "Converting matrix to VCF")
      
      # Get the uploaded file paths
      matrix_file <- input$matrix_file$datapath
      ploidy <- input$dosage2vcf_ploidy
      
      # Use a temporary file path without appending .vcf
      temp_base <- tempfile()
      
      # The output file should be temp_base.vcf
      output_name <- paste0(temp_base, ".vcf")
      
      #Status
      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 50, title = "Writing VCF")
      
      # Convert to VCF using the BIGr package
      gmatrix2vcf(Gmat.file = matrix_file,
                  ploidy = ploidy,
                  output.file = output_name,
                  dosageCount = input$dosage_counts)
      
      #Status
      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 100, title = "Complete!")
      
      return(output_name)

    }
    
  })

  # Only make available the download button when analysis is finished
  output$mybutton <- renderUI({
    if(isTruthy(vcf_out()))
      downloadButton(ns("download_d2vcf"), "Download VCF file", class = "butt")
  })


  ##This is for the DArT files conversion to VCF
  output$download_d2vcf <- downloadHandler(
    filename = function() {
      output_name <- gsub("\\.vcf$", "", input$d2v_output_name)
      paste0(output_name, ".vcf.gz")
    },
    content = function(file) {
      bgzip_compress(vcf_out(), file)
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
