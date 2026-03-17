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
                 selectInput(ns('file_type'),
                             label = 'Select File Format',
                             choices = c("DArT MADC file","DArT Dosage Reports", "Dosage Matrix"), 
                             selected = "DArT MADC file"),
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
                                  fileInput(ns("madc_file"), "Choose DArT MADC File", accept = c(".csv"), multiple = TRUE),

                                  # Error warning for multiple raw MADC files
                                  conditionalPanel(
                                    condition = "output.raw_madc_multiple_files == 'true'",
                                    ns = ns,
                                    div(style = "background-color: #f8d7da; border: 1px solid #f5c6cb; padding: 10px; margin: 10px 0; border-radius: 4px; color: #721c24;",
                                        HTML("<b>⚠ Error:</b> Cannot upload multiple MADC files when any file is raw (unprocessed).<br>
                                             Raw MADC files cannot be merged. Please upload a single raw MADC file only.<br>
                                             <i>Note: Multiple files can only be merged if ALL files are processed MADC files with fixed allele IDs.</i>")
                                    )
                                  ),

                                  # RAW MADC SPECIFIC: marker file upload and info
                                  conditionalPanel(
                                    condition = "output.is_raw_madc == 'true'",
                                    ns = ns,
                                    hr(),
                                    div(style = "background-color: #e8f4f8; border-left: 4px solid #3498db; padding: 15px; margin: 10px 0;",
                                        h5("Raw MADC Detected", style = "color: #2c3e50; margin-top: 0;"),
                                        p("This file requires preprocessing before conversion to VCF. Upload a marker mapping file to continue.")
                                    ),
                                    fileInput(ns("marker_file"),
                                              "Upload Marker File (CSV)*",
                                              accept = c(".csv")),
                                    p(style = "font-size: 0.9em; color: #666;",
                                      HTML("<b>Marker file format:</b><br>
                                           • Column 1: CloneID (must match MADC CloneIDs exactly)<br>
                                           • Column 2: Chr (e.g., Chr01)<br>
                                           • Column 3: Pos (numeric position only)<br>
                                           • Column 4 (optional): Bottom strand (Y/N)")),
                                    # Warning for missing Ref/Alt
                                    conditionalPanel(
                                      condition = "output.show_ref_alt_warning == 'true'",
                                      ns = ns,
                                      div(style = "background-color: #fff3cd; border: 1px solid #ffc107; padding: 10px; margin: 10px 0; border-radius: 4px;",
                                          HTML("<b>⚠ Warning:</b> Some markers are missing |Ref or |Alt allele IDs and will be excluded from REF/ALT extraction.")
                                      )
                                    ),
                                    conditionalPanel(
                                      condition = "output.show_raw_ref_alt_unavailable == 'true'",
                                      ns = ns,
                                      div(style = "background-color: #fff3cd; border: 1px solid #ffc107; padding: 10px; margin: 10px 0; border-radius: 4px;",
                                          HTML("<b>Note:</b> REF/ALT extraction for raw MADC files is only available when the marker file includes the optional 4th bottom-strand column.")
                                      )
                                    )
                                  ),

                                  hr(),

                                  # SHARED: SNP type selection (used by both raw and processed)
                                  radioButtons(ns("snp_type"),
                                               label = "Select Marker Type",
                                               choices = list("Target"= "target", "Target + Off-Target" = "target_off"),
                                               selected = "target"),

                                  # PROCESSED MADC SPECIFIC: species selection
                                  conditionalPanel(
                                    condition = "output.is_raw_madc == 'false'",
                                    ns = ns,
                                    selectInput(ns('species'),
                                                label = 'Species',
                                                choices = c("alfalfa","blueberry", "cranberry", "cucumber", "lettuce", "pecan", "sweetpotato", "other"),
                                                selected = "alfalfa")
                                  ),

                                  # SHARED: target_off options (used by both raw and processed)
                                  conditionalPanel(condition = "input.snp_type == 'target_off'",
                                                   ns = ns,
                                                   # Processed MADC with "other" species needs botloci upload
                                                   conditionalPanel(condition = "output.is_raw_madc == 'false' && input.species == 'other'",
                                                                    ns = ns,
                                                                    fileInput(ns("botloci_file_target_off"), "Upload bottom strand probes file (.botloci)")
                                                   ),
                                                   conditionalPanel(condition = "output.show_raw_target_off_unavailable == 'true'",
                                                                    ns = ns,
                                                                    div(style = "background-color: #fff3cd; border: 1px solid #ffc107; padding: 10px; margin: 10px 0; border-radius: 4px;",
                                                                        HTML("<b>Note:</b> Raw Target + Off-Target conversion requires the marker file's optional 4th bottom-strand column.")
                                                                    )
                                                   ),
                                                   conditionalPanel(condition = "output.is_raw_madc == 'false'",
                                                                    ns = ns,
                                                                    fileInput(ns("hapDB_file"), "Upload haplotype database file (fasta) (optional)")
                                                   ),
                                                   sliderInput(ns("cores"), "Number of CPU Cores", min = 1, max = (availableCores() - 1), value = 1, step = 1)
                                  ),

                                  # SHARED: target options (used by both raw and processed)
                                  conditionalPanel(condition = "input.snp_type == 'target'",
                                                   ns = ns,
                                                   # For processed MADC, only show ref_alt for "other" species.
                                                   # For raw MADC, only show ref_alt once the marker file proves a 4th column exists.
                                                   conditionalPanel(condition = "output.show_ref_alt_controls == 'true'",
                                                                    ns = ns,
                                                                    radioButtons(ns("ref_alt"),
                                                                                 label = "Extract REF and ALT info:",
                                                                                 choices = list("Yes"= "TRUE", "No" = "FALSE"),
                                                                                 selected = "TRUE"),
                                                                    # Processed MADC with "other" species + ref_alt needs botloci
                                                                    conditionalPanel(condition = "output.is_raw_madc == 'false' && input.species == 'other' && input.ref_alt == 'TRUE'",
                                                                                     ns = ns,
                                                                                     fileInput(ns("botloci_file_target"), "Upload bottom strand probes file (.botloci)")
                                                                    )
                                                   )
                                  )
                 ),
                 hr(),
                 textInput(ns("d2v_output_name"), "Output File Name"),
                 hr(),
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

  # Reactive values for raw MADC handling
  raw_madc_state <- reactiveValues(
    is_raw = FALSE,
    fixed_madc_path = NULL,
    marker_file_valid = FALSE,
    marker_data = NULL,
    marker_has_bottom_strand_col = FALSE,
    missing_ref_alt_markers = NULL,
    botloci_from_marker = NULL,
    multiple_files_error = FALSE
  )

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
  
  disable("download_d2vcf")
  
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


  # Detect if uploaded MADC is raw and check for multiple files
  observeEvent(input$madc_file, {
    # Check ALL uploaded files to see if any are raw
    num_files <- length(input$madc_file$datapath)
    any_raw <- FALSE

    for (i in 1:num_files) {
      if (is_raw_madc(input$madc_file$datapath[i])) {
        any_raw <- TRUE
        break
      }
    }

    raw_madc_state$is_raw <- any_raw

    # If ANY file is raw AND multiple files uploaded, this is an error
    if (any_raw && num_files > 1) {
      shinyalert(
        title = "Multiple Files Not Supported",
        text = "Cannot upload multiple MADC files when any file is raw (unprocessed). Raw MADC files cannot be merged. Please upload a single raw MADC file, or ensure all files are processed MADC files.",
        type = "error",
        confirmButtonCol = "#004192"
      )
      # Mark as invalid so conversion can't proceed
      raw_madc_state$multiple_files_error <- TRUE
    } else {
      raw_madc_state$multiple_files_error <- FALSE
    }

    # Reset dependent states when new MADC uploaded
    raw_madc_state$marker_file_valid <- FALSE
    raw_madc_state$marker_data <- NULL
    raw_madc_state$marker_has_bottom_strand_col <- FALSE
    raw_madc_state$fixed_madc_path <- NULL
    raw_madc_state$missing_ref_alt_markers <- NULL
    raw_madc_state$botloci_from_marker <- NULL
  }, ignoreNULL = TRUE)

  # Output for conditional UI (use renderText with strings)
  output$is_raw_madc <- renderText({
    if (isTRUE(raw_madc_state$is_raw)) "true" else "false"
  })
  outputOptions(output, "is_raw_madc", suspendWhenHidden = FALSE)

  # Output for multiple files error warning
  output$raw_madc_multiple_files <- renderText({
    if (isTRUE(raw_madc_state$is_raw) && isTRUE(raw_madc_state$multiple_files_error)) "true" else "false"
  })
  outputOptions(output, "raw_madc_multiple_files", suspendWhenHidden = FALSE)

  # Validate marker file when uploaded
  observeEvent(input$marker_file, {
    req(input$marker_file, raw_madc_state$is_raw)

    # Short-circuit if multiple files error
    if (isTRUE(raw_madc_state$multiple_files_error)) {
      raw_madc_state$marker_file_valid <- FALSE
      raw_madc_state$marker_data <- NULL
      raw_madc_state$marker_has_bottom_strand_col <- FALSE
      raw_madc_state$missing_ref_alt_markers <- NULL
      raw_madc_state$botloci_from_marker <- NULL
      return(NULL)
    }

    # Validate marker file format
    validation <- validate_marker_file(input$marker_file$datapath)

    if (!validation$valid) {
      shinyalert(
        title = "Invalid Marker File",
        text = validation$message,
        type = "error",
        confirmButtonCol = "#004192"
      )
      raw_madc_state$marker_file_valid <- FALSE
      raw_madc_state$marker_data <- NULL
      raw_madc_state$marker_has_bottom_strand_col <- FALSE
      raw_madc_state$missing_ref_alt_markers <- NULL
      raw_madc_state$botloci_from_marker <- NULL
      return()
    }

    # Process marker file (create lookup and extract botloci)
    processed <- process_marker_file(validation$data)
    raw_madc_state$marker_data <- validation$data
    raw_madc_state$marker_has_bottom_strand_col <- processed$has_bottom_strand_col
    raw_madc_state$botloci_from_marker <- processed$botloci_ids

    # Check for missing Ref/Alt markers (only check first file since single file required)
    marker_clones <- validation$data[,1]
    raw_madc_state$missing_ref_alt_markers <- check_missing_ref_alt(
      input$madc_file$datapath[1],
      marker_clones
    )

    raw_madc_state$marker_file_valid <- TRUE

    shinyalert(
      title = "Marker File Valid",
      text = sprintf("Marker file validated successfully. %d markers ready for processing.",
                     nrow(validation$data)),
      type = "success",
      confirmButtonCol = "#004192"
    )
  }, ignoreNULL = TRUE)

  # Output for conditional warning display
  output$show_ref_alt_warning <- renderText({
    show_warn <- isTRUE(raw_madc_state$is_raw) &&
      identical(input$snp_type, "target") &&
      isTRUE(raw_madc_state$marker_has_bottom_strand_col) &&
      identical(input$ref_alt, "TRUE") &&
      !is.null(raw_madc_state$missing_ref_alt_markers)

    if (show_warn) "true" else "false"
  })
  outputOptions(output, "show_ref_alt_warning", suspendWhenHidden = FALSE)

  output$show_ref_alt_controls <- renderText({
    show_controls <- identical(input$snp_type, "target") && (
      (isTRUE(raw_madc_state$is_raw) &&
         isTRUE(raw_madc_state$marker_file_valid) &&
         isTRUE(raw_madc_state$marker_has_bottom_strand_col)) ||
        (!isTRUE(raw_madc_state$is_raw) && identical(input$species, "other"))
    )

    if (show_controls) "true" else "false"
  })
  outputOptions(output, "show_ref_alt_controls", suspendWhenHidden = FALSE)

  output$show_raw_ref_alt_unavailable <- renderText({
    show_note <- isTRUE(raw_madc_state$is_raw) &&
      identical(input$snp_type, "target") &&
      isTRUE(raw_madc_state$marker_file_valid) &&
      !isTRUE(raw_madc_state$marker_has_bottom_strand_col)

    if (show_note) "true" else "false"
  })
  outputOptions(output, "show_raw_ref_alt_unavailable", suspendWhenHidden = FALSE)

  output$show_raw_target_off_unavailable <- renderText({
    show_note <- isTRUE(raw_madc_state$is_raw) &&
      identical(input$snp_type, "target_off") &&
      isTRUE(raw_madc_state$marker_file_valid) &&
      !isTRUE(raw_madc_state$marker_has_bottom_strand_col)

    if (show_note) "true" else "false"
  })
  outputOptions(output, "show_raw_target_off_unavailable", suspendWhenHidden = FALSE)

  observe({
    if (isTRUE(raw_madc_state$is_raw) &&
        isTRUE(raw_madc_state$marker_file_valid) &&
        !isTRUE(raw_madc_state$marker_has_bottom_strand_col) &&
        identical(input$ref_alt, "TRUE")) {
      updateRadioButtons(session, "ref_alt", selected = "FALSE")
    }
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
      get_ref_alt <- identical(input$ref_alt, "TRUE")

      # Check if raw MADC - if so, preprocess with fixMADC
      if (raw_madc_state$is_raw) {

        # Block if multiple files uploaded
        if (raw_madc_state$multiple_files_error) {
          shinyalert(
            title = "Multiple Files Not Allowed",
            text = "Raw MADC files cannot be merged. Please upload a single MADC file only.",
            type = "error",
            confirmButtonCol = "#004192"
          )
          return(NULL)
        }

        # Require marker file to be uploaded and validated
        if (!raw_madc_state$marker_file_valid) {
          shinyalert(
            title = "Marker File Required",
            text = "Please upload and validate a marker file before converting raw MADC files.",
            type = "warning",
            confirmButtonCol = "#004192"
          )
          return(NULL)
        }

        req(input$marker_file)
        req(raw_madc_state$marker_data)

        if (identical(input$snp_type, "target") &&
            get_ref_alt &&
            !isTRUE(raw_madc_state$marker_has_bottom_strand_col)) {
          shinyalert(
            title = "Bottom Strand Column Required",
            text = "Raw REF/ALT extraction requires the optional 4th marker-file column with bottom-strand values.",
            type = "warning",
            confirmButtonCol = "#004192"
          )
          return(NULL)
        }

        if (identical(input$snp_type, "target_off") &&
            !isTRUE(raw_madc_state$marker_has_bottom_strand_col)) {
          shinyalert(
            title = "Bottom Strand Column Required",
            text = "Raw Target + Off-Target conversion requires the optional 4th marker-file column with bottom-strand values.",
            type = "warning",
            confirmButtonCol = "#004192"
          )
          return(NULL)
        }

        marker_data_for_run <- raw_madc_state$marker_data
        if (identical(input$snp_type, "target") &&
            get_ref_alt &&
            !is.null(raw_madc_state$missing_ref_alt_markers)) {
          marker_data_for_run <- marker_data_for_run[
            !(marker_data_for_run[,1] %in% raw_madc_state$missing_ref_alt_markers),
            ,
            drop = FALSE
          ]
        }

        if (nrow(marker_data_for_run) == 0) {
          shinyalert(
            title = "No Valid REF/ALT Markers",
            text = "No markers with complete |Ref and |Alt allele IDs remain after filtering.",
            type = "error",
            confirmButtonCol = "#004192"
          )
          return(NULL)
        }

        marker_file_for_run <- write_marker_file(marker_data_for_run)
        processed_marker_for_run <- process_marker_file(marker_data_for_run)

        # Update progress
        updateProgressBar(session = session, id = "dosage2vcf_pb",
                          value = 5, title = "Preprocessing raw MADC file...")

        # Preprocess single MADC file using helper function
        fixed_madc_path <- tryCatch({
          preprocess_raw_madc(
            input$madc_file$datapath[1],
            marker_file_for_run
          )
        }, error = function(e) {
          shinyalert(
            title = "Preprocessing Error",
            text = paste("Error preprocessing MADC file:", e$message),
            type = "error",
            confirmButtonCol = "#004192"
          )
          NULL
        })

        if (is.null(fixed_madc_path)) {
          return(NULL)
        }

        raw_madc_state$fixed_madc_path <- fixed_madc_path

        updateProgressBar(session = session, id = "dosage2vcf_pb",
                          value = 15, title = "Raw MADC preprocessing complete")

      } else {
        # Not raw MADC - use original file paths
        raw_madc_state$fixed_madc_path <- NULL
      }

      # Determine which MADC file path to use and select botloci
      if (raw_madc_state$is_raw) {
        # RAW MADC: use preprocessed file and botloci from 4th column
        madc_to_process <- raw_madc_state$fixed_madc_path

        # Determine botloci: derive only from the marker file's optional 4th column.
        botloci <- if (isTRUE(processed_marker_for_run$has_bottom_strand_col)) {
          write_botloci_file(processed_marker_for_run$botloci_ids)
        } else {
          NULL
        }
      } else {
        # PROCESSED MADC: use original files and species-based botloci
        madc_to_process <- input$madc_file$datapath

        # Resolve botloci_file from either target_off or target branch
        botloci_input_path <- if (input$snp_type == "target_off") {
          if (!is.null(input$botloci_file_target_off)) input$botloci_file_target_off$datapath else NULL
        } else if (input$snp_type == "target") {
          if (!is.null(input$botloci_file_target)) input$botloci_file_target$datapath else NULL
        } else {
          NULL
        }

        # Select botloci: use manual upload if provided, otherwise species defaults
        if (!is.null(botloci_input_path)) {
          botloci <- botloci_input_path
        } else {
          botloci <- switch(input$species,
                            "alfalfa" = "https://raw.githubusercontent.com/Breeding-Insight/BIGapp-PanelHub/refs/heads/main/alfalfa/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_f180bp.botloci",
                            "blueberry" = "https://raw.githubusercontent.com/Breeding-Insight/BIGapp-PanelHub/refs/heads/main/blueberry/20200819-BI-Blueberry_10K_SNPs_forDArT_3K_ref_alt.botloci",
                            "cranberry" = "https://raw.githubusercontent.com/Breeding-Insight/BIGapp-PanelHub/refs/heads/main/cranberry/Cranberry_unique_alignment_126MAS_3K_54BB_rmDupTags_f180bp.botloci",
                            "cucumber" = "https://raw.githubusercontent.com/Breeding-Insight/BIGapp-PanelHub/refs/heads/main/cucumber/Cucumber_DArT3K_10192022_f180bp.botloci",
                            "lettuce" = "https://raw.githubusercontent.com/Breeding-Insight/BIGapp-PanelHub/refs/heads/main/lettuce/Lettuce_DArT3K_08172022_bait_f180bp.botloci",
                            "pecan" = "https://raw.githubusercontent.com/Breeding-Insight/BIGapp-PanelHub/refs/heads/main/pecan/Pecan_unique_alignment_top48_MAS_14K_3K_f180bp.botloci",
                            "sweetpotato" = "https://raw.githubusercontent.com/Breeding-Insight/BIGapp-PanelHub/refs/heads/main/sweetpotato/sweetpotato_20K_SNPset_f180bp_forDArT_3K_f180bp.botloci",
                            "other" = NULL)
        }
      }
      
      #Now perform conversion depending on user options
      if(input$snp_type == "target_off"){
        if (is.null(input$madc_file$datapath) | input$d2v_output_name == "") {
          shinyalert(
            title = "Missing input!",
            text = "Upload MADC and/or define a output file name",
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
        req(input$madc_file)

        # Use a temporary file path without appending .vcf
        temp_base <- tempfile()

        # For raw MADC: use preprocessed file (already in madc_to_process)
        # For processed MADC: merge if multiple files, then use
        if (raw_madc_state$is_raw) {
          # Raw MADC: single file already preprocessed
          read_madc <- madc_to_process
          updateProgressBar(session = session, id = "dosage2vcf_pb", value = 20, title = "Preparing conversion...")
        } else {
          # Processed MADC: merge if multiple
          if(length(madc_to_process) > 1){
            updateProgressBar(session = session, id = "dosage2vcf_pb", value = 15, title = "Merging MADC files")

            merged_madc <- paste0(temp_base, ".csv")

            run_ids <- sapply(strsplit(input$madc_file$name, "_"), "[[", 1)
            if(length(run_ids) == 0) run_ids <- NULL

            merge_MADCs(madc_list = as.list(madc_to_process), out_madc = merged_madc, run_ids = run_ids)
            read_madc <- merged_madc
          } else {
            read_madc <- madc_to_process
          }
        }

        # The output file should be temp_base.vcf
        output_name <- paste0(temp_base, ".vcf")

        updateProgressBar(session = session, id = "dosage2vcf_pb", value = 30, title = "Writing VCF")

        # Safe hapDB_file handling
        hap_seq_file <- if (raw_madc_state$is_raw || is.null(input$hapDB_file)) {
          NULL
        } else {
          input$hapDB_file$datapath
        }

        madc2vcf_all(madc = read_madc,
                     botloci_file = botloci,
                     hap_seq_file = hap_seq_file,
                     n.cores= input$cores,
                     rm_multiallelic_SNP = TRUE,
                     out_vcf = output_name,
                     verbose = FALSE)

        updateProgressBar(session = session, id = "dosage2vcf_pb", value = 100, title = "Complete!")

        return(output_name)
        
      } else if(input$snp_type == "target"){

        # Use a temporary file path without appending .vcf
        temp_base <- tempfile()

        # For raw MADC: use preprocessed file (already in madc_to_process)
        # For processed MADC: merge if multiple files, then use
        if (raw_madc_state$is_raw) {
          # Raw MADC: single file already preprocessed
          read_madc <- madc_to_process
          updateProgressBar(session = session, id = "dosage2vcf_pb", value = 20, title = "Preparing conversion...")
        } else {
          # Processed MADC: merge if multiple
          if(length(madc_to_process) > 1){
            updateProgressBar(session = session, id = "dosage2vcf_pb", value = 15, title = "Merging MADC files")

            merged_madc <- paste0(temp_base, ".csv")

            run_ids <- sapply(strsplit(input$madc_file$name, "_"), "[[", 1)
            if(length(run_ids) == 0) run_ids <- NULL

            merge_MADCs(madc_list = as.list(madc_to_process), out_madc = merged_madc, run_ids = run_ids)
            read_madc <- merged_madc
          } else {
            read_madc <- madc_to_process
          }
        }

        # The output file should be temp_base.vcf
        output_name <- paste0(temp_base, ".vcf")

        updateProgressBar(session = session, id = "dosage2vcf_pb", value = 30, title = "Writing VCF")
        madc2vcf_targets(read_madc, output_name, get_REF_ALT = get_ref_alt, botloci_file = botloci)

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
