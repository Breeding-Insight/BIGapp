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
                             choices = c("DArT MADC file", "Dosage Matrix"), 
                             selected = "DArT MADC file"),
                 conditionalPanel(condition = "input.file_type == 'Dosage Matrix'",
                                  ns = ns,
                                  fileInput(ns("matrix_file"), "Choose Dosage Matrix File", accept = c(".csv")),
                                  selectInput(ns("dosage_counts"), "Dosage Allele Count", choices = c("Reference","Alternate"), selected = "Reference"),
                                  numericInput(ns("dosage2vcf_ploidy"), "Species Ploidy", min = 1, value = NULL)
                 ),
                 conditionalPanel(condition = "input.file_type == 'DArT MADC file'",
                                  ns = ns,
                                  fileInput(ns("madc_file"), "Choose DArT MADC File", accept = c(".csv"), multiple = TRUE),
                                  hr(),
                                  radioButtons(ns("snp_type"),
                                               label = "Select Marker Type",
                                               choices = list("Target"= "target", "Target + Off-Target" = "target_off", "Microhaplotypes" = "multiallelic"),
                                               selected = "target"),
                                  selectInput(ns('species'),
                                              label = 'Species',
                                              choices = c("alfalfa", "blueberry", "cranberry", "cucumber", "pecan", "potato", "strawberry", "sweetpotato", "other"), 
                                              selected = "alfalfa"),
                                  conditionalPanel(condition = "input.snp_type == 'target_off'",
                                                   ns = ns,
                                                   conditionalPanel(condition = "input.species == 'other'",
                                                                    ns = ns,
                                                                    fileInput(ns("botloci_file"), "Upload bottom strand probes file (.botloci)"),
                                                                    fileInput(ns("hapDB_file"), "Upload haplotype database file (fasta) (optional)"),
                                                                    fileInput(ns("markers_info_file"), "Upload markers information (_lut.csv from HapApp) (optional)"),
                                                   ),
                                                   sliderInput(ns("cores"), "Number of CPU Cores", min = 1, max = (availableCores() - 1), value = 1, step = 1), br(),
                                                   div(
                                                       style = "text-align: left; margin-top: 10px;",
                                                        actionButton(ns("advanced_options_all"),
                                                                      label = HTML(paste(icon("cog", style = "color: #007bff;"), "Advanced Options")),
                                                                      style = "background-color: transparent; border: none; color: #007bff; font-size: smaller; text-decoration: underline; padding: 0;"
                                                    )
                                                  )
                                  ),
                                  conditionalPanel(condition = "input.snp_type == 'target'",
                                                   ns = ns,
                                                   radioButtons(ns("collapse_matches_counts"),
                                                                label = "Collapse Matches Counts:",
                                                                choices = list("Yes"= TRUE, "No" = FALSE),
                                                                selected = FALSE),  
                                                   conditionalPanel(condition = "input.species == 'other'",
                                                                    ns = ns,
                                                                    radioButtons(ns("ref_alt"),
                                                                                 label = "Extract REF and ALT info:",
                                                                                 choices = list("Yes"= TRUE, "No" = FALSE),
                                                                                 selected = TRUE),                                                                  
                                                                    conditionalPanel(condition = "input.ref_alt == 'TRUE'",
                                                                                     ns = ns,
                                                                                     fileInput(ns("botloci_file"), "Upload bottom strand probes file (.botloci)"),
                                                                                     fileInput(ns("hapDB_file"), "Upload haplotype database file (fasta) (optional)"),
                                                                                     fileInput(ns("markers_info_file"), "Upload markers information (_lut.csv from HapApp) (optional)"),
                                                                    )
                                                   )
                                  )
                 ),
                 hr(),
                 textInput(ns("d2v_output_name"), "Output File Name", placeholder = "e.g. my_output"),
                 hr(),
                 actionButton(ns("run_analysis"), "Convert File", width = "100%"), br(), br(),
                 uiOutput(ns('mybutton')),
                 
                 div(style="display:inline-block; float:right",dropdownButton(
                   HTML("<b>Input files</b>"),
                  p(downloadButton(ns('download_madc_fixed'), ""), "Fixed Allele MADC Example File"),
                   p(downloadButton(ns('download_madc'), ""), "Raw MADC Example File"),hr(),
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
        column(width = 5,
               box(title = "Status & Log", width = 12, collapsible = TRUE, status = "info",
                   progressBar(id = ns("dosage2vcf_pb"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " "),
                   hr(),
                   tags$style(HTML(paste0(
                     "#", ns("d2vcf_log"), " { max-height: 300px; overflow-y: auto; background: #1e1e1e !important; color: #d4d4d4; padding: 10px; border-radius: 4px; font-size: 12px; }"
                   ))),
                   verbatimTextOutput(ns("d2vcf_log")),
                   uiOutput(ns("download_d2vcf_log_btn"))
               )
        )
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

  # Reactive value to accumulate log messages
  d2vcf_log <- reactiveVal("")

  # Clear log when a new analysis is triggered
  observeEvent(input$run_analysis, {
    d2vcf_log("")
  }, priority = 10)

  output$d2vcf_log <- renderText({
    d2vcf_log()
  })

  output$download_d2vcf_log_btn <- renderUI({
    req(nchar(d2vcf_log()) > 0)
    downloadButton(ns("download_d2vcf_log"), "Download Log", style = "margin-top: 8px;")
  })

  output$download_d2vcf_log <- downloadHandler(
    filename = function() {
      paste0("dosage2vcf_log_", Sys.Date(), ".txt")
    },
    content = function(file) {
      writeLines(d2vcf_log(), file)
    }
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
  
  output$download_madc <- downloadHandler(
    filename = function() {
      paste0("BIGapp_MADC_Example_file.csv")
    },
    content = function(file) {
      ex <- system.file("iris_DArT_MADC.csv", package = "BIGapp")
      file.copy(ex, file)
    })
  
  output$download_madc_fixed <- downloadHandler(
    filename = function() {
      paste0("BIGapp_MADC_Fixed_Example_file.csv")
    },
    content = function(file) {
      github_path <- "https://raw.githubusercontent.com/Breeding-Insight/BIGapp-PanelHub/refs/heads/long_seq/"
      alfalfa_madc <- paste0(github_path, "test_madcs/alfalfa_madc.csv")
      download.file(alfalfa_madc, destfile = file, mode = "wb")
    })
  
  # Default model choices
  advanced_options_all <- reactiveValues(
    rm_multiallelic_SNP = TRUE,
    multiallelic_SNP_dp_thr = 0,
    multiallelic_SNP_sample_thr = 0,
    alignment_score_thr = 40,
    add_others = FALSE,
    others_max_snps = 5,
    others_rm_with_indels = TRUE
  )

    # Close popup window when user "saves options"
  observeEvent(input$save_advanced_options_all, {
    advanced_options_all$rm_multiallelic_SNP <- input$rm_multiallelic_SNP
    advanced_options_all$multiallelic_SNP_dp_thr <- input$multiallelic_SNP_dp_thr
    advanced_options_all$multiallelic_SNP_sample_thr <- input$multiallelic_SNP_sample_thr
    advanced_options_all$alignment_score_thr <- input$alignment_score_thr
    advanced_options_all$add_others <- input$add_others
    advanced_options_all$others_max_snps <- input$others_max_snps
    advanced_options_all$others_rm_with_indels <- as.logical(input$others_rm_with_indels)
    # Save other inputs as needed
    
    removeModal() # Close the modal after saving
  })

  # UI popup window for input
  observeEvent(input$advanced_options_all, {
    showModal(modalDialog(
      title = "Advanced Options",

      numericInput(
         inputId = ns("alignment_score_thr"),
         label = "Alignment Score Threshold",
         value = advanced_options_all$alignment_score_thr,
         min = 0
       ), hr(),
      # rm_multiallelic_SNP
      selectInput(
        inputId = ns("rm_multiallelic_SNP"),
        label = "Remove Multiallelic SNPs",
        choices = c("TRUE" = TRUE, "FALSE" = FALSE),
        selected = advanced_options_all$rm_multiallelic_SNP
      ),

      # Conditionally show thresholds when rm_multiallelic_SNP is FALSE
      conditionalPanel(
        condition = paste0("input['", ns("rm_multiallelic_SNP"), "'] == 'FALSE'"),
        numericInput(
          inputId = ns("multiallelic_SNP_dp_thr"),
          label = "Multiallelic SNP Depth Threshold",
          value = advanced_options_all$multiallelic_SNP_dp_thr,
          min = 0
        ),
        numericInput(
          inputId = ns("multiallelic_SNP_sample_thr"),
          label = "Multiallelic SNP Sample Threshold",
          value = advanced_options_all$multiallelic_SNP_sample_thr,
          min = 0
        )
      ), hr(),

      # add_others
      selectInput(
        inputId = ns("add_others"),
        label = "Add Others",
        choices = c("TRUE" = TRUE, "FALSE" = FALSE),
        selected = advanced_options_all$add_others
      ),

      # Conditionally show others options when add_others is TRUE
      conditionalPanel(
        condition = paste0("input['", ns("add_others"), "'] == 'TRUE'"),
        numericInput(
          inputId = ns("others_max_snps"),
          label = "Others Max SNPs",
          value = advanced_options_all$others_max_snps,
          min = 0
        ),
        selectInput(
          inputId = ns("others_rm_with_indels"),
          label = "Remove Others with Indels",
          choices = c("TRUE" = TRUE, "FALSE" = FALSE),
          selected = advanced_options_all$others_rm_with_indels
        )
      ),

      footer = tagList(
        modalButton("Close"),
        actionButton(ns("save_advanced_options_all"), "Save")
      )
    ))
  })

  vcf_out <- eventReactive(input$run_analysis,{
    # Ensure the files are uploaded
    # Missing input with red border and alerts
    if(input$file_type == "DArT MADC file"){
      req(input$madc_file)
      # First check if the MADC file is valid (a non-fixedAlleleID MADC file)
            
      # Select species botloci
      github_path <- "https://raw.githubusercontent.com/Breeding-Insight/BIGapp-PanelHub/refs/heads/long_seq/"

      botloci <- switch(input$species,
                        "alfalfa" = paste0(github_path, "alfalfa/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_f180bp.botloci"),
                        "blueberry" = paste0(github_path, "blueberry/20200819-BI-Blueberry_10K_SNPs_forDArT_3K_ref_alt.botloci"),
                        "cranberry" = paste0(github_path, "cranberry/Cranberry_unique_alignment_126MAS_3K_54BB_rmDupTags_f180bp.botloci"),
                        "cucumber" = paste0(github_path, "cucumber/Cucumber_DArT3K_10192022_f180bp.botloci"),
                        "pecan" = paste0(github_path, "pecan/Pecan_unique_alignment_top48_MAS_14K_3K_f180bp.botloci"),
                        "potato" = paste0(github_path, "potato/potato_20K_SNPset_f180bp_forDArT_3K_f180bp.botloci"),
                        "strawberry" = paste0(github_path, "strawberry/strawberry_20K_SNPset_f180bp_forDArT_3K_f180bp.botloci"),
                        "sweetpotato" = paste0(github_path, "sweetpotato/sweetpotato_20K_SNPset_f180bp_forDArT_3K_f180bp.botloci"),
                        "other" = input$botloci_file$datapath)

      microhapDB <- switch(input$species,
                        "alfalfa" = paste0(github_path, "alfalfa/alfalfa_allele_db_v001.fa"),
                        "blueberry" = paste0(github_path, "blueberry/blueberry_allele_db_v001.fa"),
                        "cranberry" = paste0(github_path, "cranberry/cranberry_allele_db_v001.fa"),
                        "cucumber" = paste0(github_path, "cucumber/cucumber_allele_db_v001.fa"),
                        "pecan" = paste0(github_path, "pecan/pecan_allele_db_v001.fa"),
                        "potato" = paste0(github_path, "potato/potato_allele_db_v001.fa"),
                        "strawberry" = paste0(github_path, "strawberry/strawberry_allele_db_v001.fa"),
                        "sweetpotato" = paste0(github_path, "sweetpotato/sweetpotato_allele_db_v000_refAlt109bp.fa"),
                        "other" = input$hapDB_file$datapath)
      
      markers_info <- switch(input$species,
                        "alfalfa" = paste0(github_path, "alfalfa/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_snpID_lut.csv"),
                        "blueberry" = paste0(github_path, "blueberry/20200819-BI-Blueberry_10K_SNPs_forDArT_3K_chrID_snpID_lut.csv"),
                        "cranberry" = paste0(github_path, "cranberry/Cranberry_unique_alignment_126MAS_3K_54BB_rmDupTags_lut.csv"),
                        "cucumber" = paste0(github_path, "cucumber/20221019-BI-Cucumber_DArTag-probe-desgin_snpID_lut.csv"),
                        "pecan" = paste0(github_path, "pecan/Pecan_unique_alignment_top48_MAS_14K_3K_snpID_lut.csv"),
                        "potato" = paste0(github_path, "potato/potato_dartag_v2_3915markers_rm7dupTags_6traitMarkers_rm1dup_snpID_lut.csv"),
                        "strawberry" = paste0(github_path, "strawberry/strawberry_9k_probe_3K_snpID_lut.csv"),
                        "sweetpotato" = paste0(github_path, "sweetpotato/sweetpotato_20K_SNPset_f180bp_forDArT_3K_snpID_lut.csv"),
                        "other" = if (!is.null(input$markers_info_file)) input$markers_info_file$datapath else NULL)

      # Guard optional inputs (may be NULL if not uploaded)
      microhapDB  <- if (length(microhapDB) == 0)  NULL else microhapDB
      markers_info <- if (length(markers_info) == 0) NULL else markers_info
      
      #Now perform conversion depending on user options

      # Shared setup for all MADC snp_type branches
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

      # Merge MADC if multiple
      if(length(input$madc_file$datapath) > 1){
        updateProgressBar(session = session, id = "dosage2vcf_pb", value = 15, title = "Merging MADC files")

        merged_madc <- paste0(temp_base, ".csv")

        run_ids <- sapply(strsplit(input$madc_file$name, "_"), "[[", 1)
        if(length(run_ids) == 0) run_ids <- NULL

        merge_MADCs(madc_list = as.list(input$madc_file$datapath), out_madc = merged_madc, run_ids = run_ids)
        read_madc <- merged_madc
      } else read_madc <- input$madc_file$datapath

      # The output file should be temp_base.vcf
      output_name <- paste0(temp_base, ".vcf")

      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 30, title = "Writing VCF")

      log_lines <- character(0)
      tryCatch(
        withCallingHandlers(
          {
            if(input$snp_type == "target_off"){
              madc2vcf_all(madc = read_madc,
                           botloci_file = botloci,
                           hap_seq_file = microhapDB,
                           markers_info = markers_info,
                           n.cores= input$cores,
                           rm_multiallelic_SNP = as.logical(advanced_options_all$rm_multiallelic_SNP),
                           multiallelic_SNP_dp_thr = advanced_options_all$multiallelic_SNP_dp_thr,
                           multiallelic_SNP_sample_thr = advanced_options_all$multiallelic_SNP_sample_thr,
                           alignment_score_thr = advanced_options_all$alignment_score_thr,
                           add_others = as.logical(advanced_options_all$add_others),
                           others_max_snps = advanced_options_all$others_max_snps,
                           others_rm_with_indels = as.logical(advanced_options_all$others_rm_with_indels),
                           out_vcf = output_name,
                           verbose = TRUE)
            } else if(input$snp_type == "target"){
              madc2vcf_targets(read_madc, output_name, get_REF_ALT = as.logical(input$ref_alt), botloci_file = botloci,
                               markers_info = markers_info, collapse_matches_counts = input$collapse_matches_counts)
            } else if(input$snp_type == "multiallelic"){
              madc2vcf_multi(
                madc_file    = read_madc,
                botloci_file = botloci,
                outfile      = output_name,
                markers_info = markers_info,
                ploidy       = 4L,
                verbose      = TRUE
              )
            }
          },
          message = function(m) {
            log_lines <<- c(log_lines, conditionMessage(m))
            invokeRestart("muffleMessage")
          },
          warning = function(w) {
            log_lines <<- c(log_lines, paste0("Warning: ", conditionMessage(w), "\n"))
            shinyalert(
              title = "Warning",
              text = conditionMessage(w),
              size = "s",
              type = "warning",
              showConfirmButton = TRUE,
              confirmButtonText = "OK",
              confirmButtonCol = "#004192"
            )
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          log_lines <<- c(log_lines, paste0("Error: ", conditionMessage(e), "\n"))
          shinyalert(
            title = "Error",
            text = conditionMessage(e),
            size = "s",
            type = "error",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192"
          )
        }
      )

      d2vcf_log(paste(log_lines, collapse = ""))
      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 100, title = "Complete!")
      return(output_name)
      
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
