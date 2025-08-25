#' Filtering UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom bs4Dash tabBox valueBoxOutput
#' @importFrom shiny NS tagList
#' @importFrom purrr map set_names
#' @importFrom stringr str_split
#' @importFrom shinyjs enable disable useShinyjs
#'
#' @import dplyr
#' @import shinydisconnect
#'
#'
mod_Filtering_ui <- function(id){
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
             box(width = 12,
                 title = "Quality Filtering", status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                 fileInput(ns("updog_rdata"),"Choose VCF File", accept = c(".vcf",".gz")),
                 textInput(ns("filter_output_name"), "Output File Name"),
                 numericInput(ns("filter_ploidy"),"Species Ploidy", min = 0, value = NULL),
                 numericInput(ns("filter_maf"),"MAF filter", min = 0, max=1, value = 0.05, step = 0.01),
                 numericInput(ns("size_depth"),"Min Read Depth (Marker per Sample)", min = 0, max = 300, value = 10, step = 1),
                 numericInput(ns("snp_miss"),"Remove SNPs with >= % missing data", min = 0, max = 100, value = 50, step = 1),
                 numericInput(ns("sample_miss"),"Remove Samples with >= % missing data", min = 0, max = 100, value = 50, step = 1),
                 "Updog Filtering Parameters",
                 checkboxInput(ns("use_updog"), "Use Updog Filtering Parameters?", value = FALSE),
                 conditionalPanel(
                   condition = "input.use_updog == true", ns = ns,
                   div(
                     numericInput(ns("OD_filter"), "Max OD (Updog filter)", min = 0, value = 0.05),
                     sliderInput(ns("Bias"), "Bias (Updog filter)", min = 0, max = 10, value = c(0.5, 2), step = 0.1),
                     numericInput(ns("Prop_mis"), "Max Prop_mis (Updog filter)", min = 0, max = 1, value = 0.05, step = 0.05),
                     numericInput(ns("maxpostprob_filter"), "Minimum maxpostprob (Updog filter)", min = 0, value = 0.5, step = 0.1)
                   )
                 ),
                 actionButton(ns("run_filters"), "Apply Filters"),
                 uiOutput(ns("mybutton")),
                 div(style="display:inline-block; float:right",dropdownButton(
                   HTML("<b>Input files</b>"),
                   p(downloadButton(ns('download_vcf'),""), "VCF Example File"),
                   p(HTML("<b>Parameters description:</b>"), actionButton(ns("goPar"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>Results description:</b>"), actionButton(ns("goRes"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>How to cite:</b>"), actionButton(ns("goCite"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   actionButton(ns("filtering_summary"), "Summary"),
                   circle = FALSE,
                   status = "warning",
                   icon = icon("info"), width = "300px",
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
             )
      ),
      column(width = 6,
             uiOutput(ns("din_tabs")),
      ),
      column(width = 3,
             valueBoxOutput(ns("snp_retained_box"), width = NULL),
             valueBoxOutput(ns("snp_removed_box"), width = NULL),
             box(title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE,
                 selectInput(ns("vcf_type"), "Histogram Data", choices = c("Unfiltered VCF", "Filtered VCF"), selected = "Unfiltered VCF"),
                 sliderInput(ns("hist_bins"),"Histogram Bins", min = 1, max = 1200, value = c(50), step = 1), width = NULL,
                 div(style="display:inline-block; float:left", dropdownButton(
                   selectInput(inputId = ns('filter_hist'), label = 'Figure', choices = c("Bias Histogram",
                                                                                          "OD Histogram",
                                                                                          "Prop_mis Histogram",
                                                                                          "SNP_mis",
                                                                                          "Sample_mis")),
                   selectInput(inputId = ns('image_type'), label = 'File Type', choices = c("jpeg","tiff","png","svg"), selected = "jpeg"),
                   sliderInput(inputId = ns('image_res'), label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                   sliderInput(inputId = ns('image_width'), label = 'Width', value = 8, min = 1, max = 20, step=0.5),
                   sliderInput(inputId = ns('image_height'), label = 'Height', value = 5, min = 1, max = 20, step = 0.5),
                   downloadButton(ns("download_filter_hist"), "Save"),
                   circle = FALSE,
                   status = "danger", label = "Save",
                   icon = icon("floppy-disk"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see inputs!")
                 ))
             ),
             box(title = "Status", width =12, collapsible = TRUE, status = "info",
                 progressBar(id = ns("pb_filter"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
             ),
             # A placeholder for the download button. It will be rendered in the shinyalert modal.
             uiOutput(ns("download_ui_placeholder")) 
      )
    )
  )
}

#' Filtering Server Functions
#'
#' @import vcfR
#' @import BIGr
#' @importFrom shinyWidgets virtualSelectInput
#' @importFrom shinyjs enable disable useShinyjs
#' @importFrom graphics abline axis hist
#'
#' @noRd
mod_Filtering_server <- function(input, output, session, parent_session){

  ns <- session$ns

  # Help links
  observeEvent(input$goPar, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "VCF_Filtering_tabset",
                      selected = "VCF_Filtering_par")
    # expand specific box
    updateBox(id = "VCF_Filtering_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goRes, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "VCF_Filtering_tabset",
                      selected = "VCF_Filtering_results")
    # expand specific box
    updateBox(id = "VCF_Filtering_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goCite, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "VCF_Filtering_tabset",
                      selected = "VCF_Filtering_cite")
    # expand specific box
    updateBox(id = "VCF_Filtering_box", action = "toggle", session = parent_session)
  })
  
  
  ## Advanced options popup
  #Default model choices
  advanced_options <- reactiveValues(
    sample_list = NULL,
    remove_list = NULL,
    remove_file = NULL
  )
  
  #List the ped file name if previously uploaded
  output$uploaded_file_name <- renderText({
    if (!is.null(advanced_options$remove_file)) {
      paste("Previously uploaded file:", advanced_options$remove_file$name)
    } else {
      ""  # Return an empty string if no file has been uploaded
    }
  })
  
  #Get list of sample names from VCF file
  observeEvent(input$updog_rdata, {
    #### VCF sanity check
    checks <- vcf_sanity_check(input$updog_rdata$datapath, 
                               max_markers = 16000, 
                               depth_support_fields = c("DP", "AD", "RA"))
    
    error_if_false <- c(
      "VCF_header", "VCF_columns", "unique_FORMAT", "GT",
      "samples", "chrom_info", "pos_info", "VCF_compressed", "allele_counts"
    )
    
    error_if_true <- c(
      "multiallelics", "phased_GT",  
      "duplicated_samples", "duplicated_markers"
    )
    
    warning_if_false <- c("ref_alt","max_markers")
    
    checks_result <- vcf_sanity_messages(checks, 
                                         error_if_false, 
                                         error_if_true, 
                                         warning_if_false = warning_if_false, 
                                         warning_if_true = NULL)
    
    print(checks)
    print(checks_result)
    if(checks_result) return() # Stop the analysis if checks fail
    #########
    
    
    #populate preview_data
    preview_vcf <- read.vcfR(input$updog_rdata$datapath, verbose = FALSE, nrows = 1)
    
    #Get names
    advanced_options$sample_list <- names(data.frame(preview_vcf@gt, check.names=FALSE)[,-1])
    
    rm(preview_vcf)
  })
  
  #UI popup window for input
  observeEvent(input$advanced_options, {
    showModal(modalDialog(
      title = "Advanced Options",
      size = "l",
      div(
        h4(
          "Remove Samples From VCF",
          style = "font-size: 18px; color: black;" # Smaller and purple
        ),
      ),
      fluidRow(
        column(
          width = 5,
          virtualSelectInput(
            inputId = ns("remove_list"),
            label = "Select Samples to Remove",
            choices = advanced_options$sample_list,
            selected = advanced_options$remove_list,
            showValueAsTags = TRUE,
            search = TRUE,
            multiple = TRUE
          )
        ),
        column(
          width = 1,
          h4(
            "or",
            style = "font-size: 20px; color: blue;"
          )
        ),
        column(
          width = 5,
          div(
            fileInput(ns("remove_file"), "Upload Sample File", accept = ".txt"),
            conditionalPanel(
              condition = "output.uploaded_file_name !== ''",
              textOutput(ns("uploaded_file_name"))
            )
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
    #Only close the window if one of the options has been selected
    if (!is.null(input$remove_list) && !is.null(input$remove_file)) {
      showNotification("Please select only one method (list or file) to remove samples. Please refresh BIGapp", type = "warning")
    } else {
      advanced_options$remove_list <- input$remove_list
      advanced_options$remove_file <- input$remove_file
      # Save other inputs as needed
      
      removeModal()
    }
    
  })

  #vcf
  filtering_files <- reactiveValues(
    raw_vcf_df = NULL,
    sample_miss_df = NULL,
    snp_miss_df = NULL,
    raw_snp_miss_df = NULL,
    raw_sample_miss_df = NULL,
    maf_df = NULL,
    raw_maf_df = NULL,
    format_fields = NULL,
    removed_names = NULL

  )

  # Function to choose user selected dataset
  current_hist_data <- reactive({
    req(input$vcf_type)

    # Switch between 'pre-filtered' or 'filtered' data based on user choice
    if (input$vcf_type == "Unfiltered VCF") {
      list(
        snp_miss = filtering_files$raw_snp_miss_df,
        sample_miss = filtering_files$raw_sample_miss_df,
        maf_data = filtering_files$raw_maf_df
      )
    } else {
      list(
        snp_miss = filtering_files$snp_miss_df,
        sample_miss = filtering_files$sample_miss_df,
        maf_data = filtering_files$maf_df
      )
    }
  })

  #Reactive boxes
  output$snp_retained_box <- renderValueBox({
    valueBox(
      value = 0,
      subtitle = "SNPs Retained",
      icon = icon("dna"),
      color = "info"
    )
  })

  output$snp_removed_box <- renderValueBox({
    valueBox(
      value = 0,
      subtitle = "Percent SNPs Removed",
      icon = icon("filter-circle-xmark"),
      color = "info"
    )
  })

  vcf <- eventReactive(input$run_filters, {

    # Ensure the files are uploaded
    # Missing input with red border and alerts
    toggleClass(id = "filter_ploidy", class = "borderred", condition = (is.na(input$filter_ploidy) | is.null(input$filter_ploidy) | input$filter_ploidy == ""))
    toggleClass(id = "filter_output_name", class = "borderred", condition = (is.na(input$filter_output_name) | is.null(input$filter_output_name) | input$filter_output_name == ""))

    if (is.null(input$updog_rdata$datapath)) {
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

    req(input$filter_ploidy, input$filter_output_name,input$updog_rdata)

    #Status
    updateProgressBar(session = session, id = "pb_filter", value = 10, title = "Processing VCF file")

    #Input file
    vcf <- read.vcfR(input$updog_rdata$datapath, verbose = FALSE)

    # Identify if have updog parameters
    format_fields <- unique(vcf@gt[,1])
    info_fields <- vcf@fix[1,8]
    updog_par <- grepl("MPP", format_fields) & grepl("PMC", info_fields) & grepl("BIAS", info_fields) & grepl("OD", info_fields)
    filtering_files$format_fields <- updog_par

    if(length(updog_par) > 1) {
      shinyalert(
        title = "Malformed VCF",
        text = "Make sure all your markers have the same information in the FORMAT field.",
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
    
    if (input$use_updog & updog_par) {
      # Use Updog filtering parameters
      OD_filter <- as.numeric(input$OD_filter)
      Prop_mis <- as.numeric(input$Prop_mis)
      Bias_min <- as.numeric(input$Bias[1])
      Bias_max <- as.numeric(input$Bias[2])
      max_post <- as.numeric(input$maxpostprob_filter)

      # Perform filtering with Updog parameters
      # (insert your filtering code here)
    } else {
      # Do not use Updog filtering parameters
      OD_filter = NULL
      Prop_mis = NULL
      Bias_min = NULL
      Bias_max = NULL
      max_post = NULL
    }

    #Variables
    size_depth <- input$size_depth
    output_name <- input$filter_output_name
    snp_miss <- input$snp_miss / 100
    sample_miss <- input$sample_miss / 100
    ploidy <- as.numeric(input$filter_ploidy)
    maf_filter <- input$filter_maf

    #Starting SNPs
    starting_snps <- nrow(vcf)
    output$snp_removed_box <- renderValueBox({
      valueBox(
        value = round(((starting_snps - final_snps)/starting_snps*100),1),
        subtitle = "Percent SNPs Removed",
        icon = icon("dna"),
        color = "info"
      )
    })

    #export INFO dataframe
    filtering_files$raw_vcf_df <- data.frame(vcf@fix)

    #Pb
    updateProgressBar(session = session, id = "pb_filter", value = 40, title = "Filtering VCF file")

    #Filtering
    #Raw VCF info
    gt_matrix <- extract.gt(filterVCF(vcf, ploidy = ploidy,filter.DP = as.numeric(size_depth), output.file = NULL), element = "GT", as.numeric = FALSE)
    starting_samples <- colnames(gt_matrix)
    filtering_files$raw_snp_miss_df <- rowMeans(is.na(gt_matrix)) #SNP missing values
    filtering_files$raw_sample_miss_df <- as.numeric(colMeans(is.na(gt_matrix))) #Sample missing values

    rm(gt_matrix) #Remove gt matrix
    
    #Remove the samples if any are manually selected from advanced options
    if (!is.null(advanced_options$remove_list)) {
      advanced_options$remove_file <- NULL #Prioritize manually selected samples if a file was also uploaded (add a user warning if both are uploaded in model)
      
      vcf_temp <- subset_vcf(vcf, remove.sample.list = advanced_options$remove_list)
      vcf <- vcf_temp[[1]]
      removed_samples <- vcf_temp[[2]]
      rm(vcf_temp)
    } else if (!is.null(advanced_options$remove_file)) {
      
      #Remove the samples
      vcf_temp <- subset_vcf(vcf, remove.sample.file = advanced_options$remove_file$datapath)
      vcf <- vcf_temp[[1]]
      removed_samples <- vcf_temp[[2]]
      rm(vcf_temp)
    } else {
      removed_samples <- 0
    }

    # Filtered VCF
    vcf <- filterVCF(vcf.file = vcf,
                     ploidy=ploidy,
                     output.file=NULL,
                     filter.OD = OD_filter,
                     filter.BIAS.min = Bias_min,
                     filter.BIAS.max = Bias_max,
                     filter.DP = as.numeric(size_depth),
                     filter.PMC = Prop_mis,
                     filter.SAMPLE.miss = as.numeric(sample_miss),
                     filter.SNP.miss = as.numeric(snp_miss),
                     filter.MAF = as.numeric(maf_filter),
                     filter.MPP = max_post)

    if (length(vcf@gt) == 0) {
      shinyalert(
        title = "All markers were filtered out",
        text = "Loose the parameters to access results in this tab",
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

    #Getting missing data information
    #Add support for genotype matrix filtering?
    #Pb
    updateProgressBar(session = session, id = "pb_filter", value = 50, title = "Calculating Missing Data")

    gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
    filtering_files$snp_miss_df <- rowMeans(is.na(gt_matrix)) #SNP missing values
    filtering_files$sample_miss_df <- as.numeric(colMeans(is.na(gt_matrix))) #Sample missing values
    final_samples <- colnames(gt_matrix)
    rm(gt_matrix) #Remove gt matrix

    #Pb
    updateProgressBar(session = session, id = "pb_filter", value = 80, title = "Exporting Filtered VCF")

    #Get final_snps
    final_snps <- nrow(vcf)
    #Updating value boxes
    output$snp_retained_box <- renderValueBox({
      valueBox(
        value = final_snps,
        subtitle = "SNPs Retained",
        icon = icon("dna"),
        color = "info"
      )
    })

    #User warning if samples were removed during filtering
    sample_removed <- length(starting_samples) - length(final_samples)
    removed_names <- setdiff(starting_samples, final_samples)
    filtering_files$removed_names <- removed_names
    
    # Define the download handler
    output$download_removed_samples <- downloadHandler(
      filename = function() {
        "removed_samples.txt"
      },
      content = function(file) {
        if (!is.null(filtering_files$removed_names)) {
          writeLines(filtering_files$removed_names, file)
        }
      }
    )
    
    if (sample_removed > 0 && removed_samples == 0) {
      showModal(modalDialog(
        title = "Samples Filtered",
        footer = tagList(
          modalButton("OK"),
          # Use a downloadButton instead of a manual <a> tag
          downloadButton(ns("download_removed_samples"), "Save Sample List")
        ),
        paste(sample_removed, "samples were removed during filtering.")
      ))
    } else if (sample_removed > 0 && removed_samples > 0) {
      shinyalert(
        title = "Samples Filtered",
        text = paste(sample_removed, "samples were removed during filtering.\n",removed_samples,"of",sample_removed,"were manually removed by user input."),
        size = "s",
        closeOnEsc = TRUE,
        closeOnClickOutside = FALSE,
        html = TRUE,
        type = "info",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
        confirmButtonCol = "#004192",
        showCancelButton = FALSE,
        animation = TRUE
      )
    }

    # Status
    updateProgressBar(session = session, id = "pb_filter", value = 100, title = "Finished!")


    vcf
  })

  #Update plots
  output$din_tabs <- renderUI({

    if (input$run_filters == 0) {
      tabBox(width =12, collapsible = FALSE, status = "info",
             id = "updog_tab", height = "600px",
             tabPanel("Results", p("Upload VCF file to access results in this section."))
      )
    } else {
      if (!any(is.null(filtering_files$format_fields)) && any(filtering_files$format_fields == TRUE) && input$vcf_type == "Unfiltered VCF") {
        # Include "Bias Histogram", "OD Histogram", and "Prop_mis Histogram" for Pre-Filtered VCF
        tabBox(
          width = 12, collapsible = FALSE, status = "info",
          id = "updog_tab", height = "600px",
          tabPanel("Bias Histogram", icon = icon("image"), plotOutput(ns("bias_hist"), height = '550px')),
          tabPanel("OD Histogram", icon = icon("image"), plotOutput(ns("od_hist"), height = '550px')),
          tabPanel("Prop_mis Histogram", icon = icon("image"), plotOutput(ns("maxpostprob_hist"), height = '550px')),
          tabPanel("SNP_miss", icon = icon("image"), plotOutput(ns("missing_snp_hist"), height = '550px')),
          tabPanel("Sample_miss", icon = icon("image"), plotOutput(ns("missing_sample_hist"), height = '550px'))
        )
      } else {
        # For Filtered VCF or non-updog VCF, only include "SNP_miss" and "Sample_miss"
        tabBox(
          width = 12, collapsible = FALSE, status = "info",
          id = "updog_tab", height = "600px",
          tabPanel("SNP_miss", icon = icon("image"), plotOutput(ns("missing_snp_hist"), height = '550px')),
          tabPanel("Sample_miss", icon = icon("image"), plotOutput(ns("missing_sample_hist"), height = '550px'))
        )
      }
    }
  })

  # Only make available the download button when analysis is finished
  output$mybutton <- renderUI({
    if(isTruthy(vcf()))
      downloadButton(ns("start_updog_filter"), "Download VCF file", class = "butt")
  })

  #Updog filtering
  output$start_updog_filter <- downloadHandler(
    filename = function() {
      output_name <- gsub("\\.vcf$", "", input$filter_output_name)
      paste0(output_name, ".vcf.gz")
    },
    content = function(file) {

      #Writing file
      temp_file <- tempfile(fileext = ".vcf.gz")
      write.vcf(vcf(), file = temp_file)

      # Avoid exporting gziped instead of bgziped
      gunzip(temp_file)
      temp_file <- gsub("\\.gz$", "", temp_file)

      # Check if the VCF file was created
      if (file.exists(temp_file)) {
        cat("VCF file created successfully.\n")

        # Move the file to the path specified by 'file'
        bgzip_compress(temp_file, file)

        # Delete the temporary file
        unlink(temp_file)
      } else {
        stop("Error: Failed to create the VCF file.")
      }

    }
  )

  #Download figures for VCF Filtering
  output$download_filter_hist <- downloadHandler(

    filename = function() {
      if (input$image_type == "jpeg") {
        paste("VCF-histogram-", Sys.Date(), ".jpg", sep="")
      } else if (input$image_type == "png") {
        paste("VCF-histogram-", Sys.Date(), ".png", sep="")
      } else if (input$image_type == "tiff") {
        paste("VCF-histogram-", Sys.Date(), ".tiff", sep="")
      } else {
        paste("VCF-histogram-", Sys.Date(), ".svg", sep="")
      }
    },
    content = function(file) {
      req(input$image_type)

      if (input$image_type == "jpeg") {
        jpeg(file, width = as.numeric(input$image_width), height = as.numeric(input$image_height), res= as.numeric(input$image_res), units = "in")
      } else if (input$image_type == "png") {
        png(file, width = as.numeric(input$image_width), height = as.numeric(input$image_height), res= as.numeric(input$image_res), units = "in")
      } else if (input$image_type == "tiff") {
        tiff(file, width = as.numeric(input$image_width), height = as.numeric(input$image_height), res= as.numeric(input$image_res), units = "in")
      } else {
        svg(file, width = as.numeric(input$image_width), height = as.numeric(input$image_height))
      }

      # Conditional plotting based on input selection
      req(filtering_output$df, filtering_files)
      if (input$filter_hist == "Bias Histogram") {

        hist(as.numeric(filtering_output$df$BIAS),
             main = "Unfiltered SNP bias histogram",
             xlab = "bias",
             ylab = "SNPs",
             col = "lightblue",
             border = "black",
             xlim = c(0,5),
             breaks = as.numeric(input$hist_bins))
        axis(1, at = seq(0, 5, by = .2), labels = rep("", length(seq(0, 5, by = 0.2))))  # Add ticks
        abline(v = mean(as.numeric(filtering_output$df$BIAS)), col = "red", lty = 2)  # Mean line
        abline(v = median(as.numeric(filtering_output$df$BIAS)), col = "green", lty = 2)  # Median line
        abline(v = 0.5, col = "black", lty = 2)  # proposed lower line
        abline(v = 2, col = "black", lty = 2)  # proposed upper line
        legend("topright", legend=c("mean", "median", "suggested threshold"),
               col=c("red", "green","black"), lty=2, cex=0.8)

      } else if (input$filter_hist == "OD Histogram") {

        #Plot
        hist(as.numeric(filtering_output$df$OD),
             main = "Unfiltered SNP overdispersion parameter histogram",
             xlab = "OD",
             ylab = "SNPs",
             col = "lightblue",
             border = "black",
             xlim = c(0,0.6),
             breaks = as.numeric(input$hist_bins))
        axis(1, at = seq(0, 0.6, by = .01), labels = rep("", length(seq(0, 0.6, by = 0.01))))  # Add ticks
        abline(v = 0.05, col = "black", lty = 2)  # proposed filter by updog

        # Add vertical lines
        abline(v = mean(as.numeric(filtering_output$df$OD)), col = "red", lty = 2)  # Mean line
        abline(v = median(as.numeric(filtering_output$df$OD)), col = "green", lty = 2)  # Median line
        abline(v = 0.05, col = "black", lty = 2)  # proposed filter by updog
        legend("topright", legend=c("mean", "median", "suggested threshold"),
               col=c("red", "green","black"), lty=2, cex=0.8)

      } else if (input$filter_hist == "Prop_mis Histogram") {

        hist(as.numeric(filtering_output$df$PMC),
             main = "The estimated proportion of individuals misclassified in the SNP from updog",
             xlab = "Proportion of Misclassified Genotypes per SNP",
             ylab = "Number of SNPs",
             col = "lightblue",
             border = "black",
             xlim = c(0,1),
             breaks = as.numeric(input$hist_bins))
        axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

        # Add vertical lines
        abline(v = mean(as.numeric(filtering_output$df$PMC)), col = "red", lty = 2)  # Mean line
        abline(v = median(as.numeric(filtering_output$df$PMC)), col = "green", lty = 2)  # Median line
        abline(v = quantile(as.numeric(filtering_output$df$PMC), 0.95), col = "blue", lty = 2)
        legend("topright", legend=c("mean", "median", "quantile"),
               col=c("red", "green","blue"), lty=2, cex=0.8)

      } else if (input$filter_hist == "SNP_mis") {

        #Histogram
        hist(
          as.numeric(current_hist_data()$snp_miss),
          main = paste("SNP Missing Data -", input$vcf_type),
          xlab = "Proportion of Missing Data per SNP",
          ylab = "Number of SNPs",
          col = "lightblue",
          border = "black",
          xlim = c(0,1),
          breaks = as.numeric(input$hist_bins))
        axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

        # Add vertical lines
        abline(v = mean(as.numeric(current_hist_data()$snp_miss)), col = "red", lty = 2)  # Mean line
        abline(v = median(as.numeric(current_hist_data()$snp_miss)), col = "green", lty = 2)  # Median line
        abline(v = quantile(as.numeric(current_hist_data()$snp_miss), 0.95), col = "blue", lty = 2)
        legend("topright", legend=c("mean", "median", "quantile"),
               col=c("red", "green","blue"), lty=2, cex=0.8)

      } else if (input$filter_hist == "Sample_mis") {

        hist(
          as.numeric(current_hist_data()$sample_miss),
          main = paste("Sample Missing Data -", input$vcf_type),
          xlab = "Proportion of Missing Data per Sample",
          ylab = "Number of Samples",
          col = "lightblue",
          border = "black",
          xlim = c(0,1),
          breaks = as.numeric(input$hist_bins))
        axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

        # Add vertical lines
        abline(v = mean(as.numeric(current_hist_data()$sample_miss)), col = "red", lty = 2)  # Mean line
        abline(v = median(as.numeric(current_hist_data()$sample_miss)), col = "green", lty = 2)  # Median line
        abline(v = quantile(as.numeric(current_hist_data()$sample_miss), 0.95), col = "blue", lty = 2)
        legend("topright", legend=c("mean", "median", "quantile"),
               col=c("red", "green","blue"), lty=2, cex=0.8)
      }
      dev.off()
    }
  )

  # Commented code
  ##Updog file stats
  #Consider Extracting the GT info or UD info if present as a datafrfame,
  #Obtaining the info in the INFO column as it's own dataframe with a column for each value
  #Then remove the VCF file and use the remaining dataframes for producing the figures

  filtering_output <- reactiveValues(df = NULL)

  observeEvent(filtering_files$raw_vcf_df, {

    # Apply the function to each row and bind the results into a new dataframe
    new_df <- data.frame(filtering_files$raw_vcf_df) %>%
      mutate(INFO_list = map(INFO, split_info_column)) %>%
      unnest_wider(INFO_list)

    #Save df to reactive value
    filtering_output$df <- new_df


    ##Make plots
    #Number of SNPs
    nrow(filtering_files$raw_vcf_df)

    ###Bias

    #Histogram
    if(any(grepl("BIAS", colnames(new_df)))){
      output$bias_hist <- renderPlot({
        hist(as.numeric(new_df$BIAS),
             main = "Unfiltered SNP bias histogram",
             xlab = "bias",
             ylab = "SNPs",
             col = "lightblue",
             border = "black",
             xlim = c(0,5),
             breaks = as.numeric(input$hist_bins))
        axis(1, at = seq(0, 5, by = .2), labels = rep("", length(seq(0, 5, by = 0.2))))  # Add ticks
        abline(v = mean(as.numeric(new_df$BIAS)), col = "red", lty = 2)  # Mean line
        abline(v = median(as.numeric(new_df$BIAS)), col = "green", lty = 2)  # Median line
        abline(v = 0.5, col = "black", lty = 2)  # proposed lower line
        abline(v = 2, col = "black", lty = 2)  # proposed upper line
        legend("topright", legend=c("mean", "median", "suggested threshold"),
               col=c("red", "green","black"), lty=2, cex=0.8)
      })
    }

    ###OD
    if(any(grepl("OD", colnames(new_df)))){

      quantile(as.numeric(new_df$OD), 0.95)
      #Histogram
      output$od_hist <- renderPlot({
        hist(as.numeric(new_df$OD),
             main = "Unfiltered SNP overdispersion parameter histogram",
             xlab = "OD",
             ylab = "SNPs",
             col = "lightblue",
             border = "black",
             xlim = c(0,0.6),
             breaks = as.numeric(input$hist_bins))
        axis(1, at = seq(0, 0.6, by = .01), labels = rep("", length(seq(0, 0.6, by = 0.01))))  # Add ticks
        abline(v = 0.05, col = "black", lty = 2)  # proposed filter by updog

        # Add vertical lines
        abline(v = mean(as.numeric(new_df$OD)), col = "red", lty = 2)  # Mean line
        abline(v = median(as.numeric(new_df$OD)), col = "green", lty = 2)  # Median line
        abline(v = 0.05, col = "black", lty = 2)  # proposed filter by updog
        legend("topright", legend=c("mean", "median", "suggested threshold"),
               col=c("red", "green","black"), lty=2, cex=0.8)

      })
    }

    ##MAXPOSTPROB

    #Histogram
    if(any(grepl("PMC", colnames(new_df)))){

      output$maxpostprob_hist <- renderPlot({

        #Histogram
        hist(as.numeric(new_df$PMC),
             main = "The estimated proportion of individuals misclassified in the SNP from updog",
             xlab = "Proportion of Misclassified Genotypes per SNP",
             ylab = "Number of SNPs",
             col = "lightblue",
             border = "black",
             xlim = c(0,1),
             breaks = as.numeric(input$hist_bins))
        axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

        # Add vertical lines
        abline(v = mean(as.numeric(new_df$PMC)), col = "red", lty = 2)  # Mean line
        abline(v = median(as.numeric(new_df$PMC)), col = "green", lty = 2)  # Median line
        abline(v = quantile(as.numeric(new_df$PMC), 0.95), col = "blue", lty = 2)
        legend("topright", legend=c("mean", "median", "quantile"),
               col=c("red", "green","blue"), lty=2, cex=0.8)

      })
    }

    #Missing data
    output$missing_snp_hist <- renderPlot({
      req(current_hist_data()$snp_miss)

      #Histogram
      hist(
        as.numeric(current_hist_data()$snp_miss),
           main = paste("SNP Missing Data -", input$vcf_type),
           xlab = "Proportion of Missing Data per SNP",
           ylab = "Number of SNPs",
           col = "lightblue",
           border = "black",
           xlim = c(0,1),
           breaks = as.numeric(input$hist_bins))
      axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

      # Add vertical lines
      abline(v = mean(as.numeric(current_hist_data()$snp_miss)), col = "red", lty = 2)  # Mean line
      abline(v = median(as.numeric(current_hist_data()$snp_miss)), col = "green", lty = 2)  # Median line
      abline(v = quantile(as.numeric(current_hist_data()$snp_miss), 0.95), col = "blue", lty = 2)
      legend("topright", legend=c("mean", "median", "quantile"),
             col=c("red", "green","blue"), lty=2, cex=0.8)
    })

    output$missing_sample_hist <- renderPlot({

      #Histogram
      hist(
        as.numeric(current_hist_data()$sample_miss),
           main = paste("Sample Missing Data -", input$vcf_type),
           xlab = "Proportion of Missing Data per Sample",
           ylab = "Number of Samples",
           col = "lightblue",
           border = "black",
           xlim = c(0,1),
           breaks = as.numeric(input$hist_bins))
      axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

      # Add vertical lines
      abline(v = mean(as.numeric(current_hist_data()$sample_miss)), col = "red", lty = 2)  # Mean line
      abline(v = median(as.numeric(current_hist_data()$sample_miss)), col = "green", lty = 2)  # Median line
      abline(v = quantile(as.numeric(current_hist_data()$sample_miss), 0.95), col = "blue", lty = 2)
      legend("topright", legend=c("mean", "median", "quantile"),
             col=c("red", "green","blue"), lty=2, cex=0.8)
    })

    ##Read Depth (I would prefer that this show the mean depth for SNPs or Samples instead of all loci/sample cells)
    #quantile(as.numeric(new_df$DP), 0.95)
  })

  output$download_vcf <- downloadHandler(
    filename = function() {
      paste0("BIGapp_VCF_Example_file.vcf.gz")
    },
    content = function(file) {
      ex <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")
      file.copy(ex, file)
    })

  ##Summary Info
  filtering_summary_info <- function() {
    #Handle possible NULL values for inputs
    genotype_file_name <- if (!is.null(input$updog_rdata$name)) input$updog_rdata$name else "No file selected"
    selected_ploidy <- if (!is.null(input$filter_ploidy)) as.character(input$filter_ploidy) else "Not selected"

    #Print the summary information
    cat(
      "BIGapp VCF Filtering Summary\n",
      "\n",
      paste0("Date: ", Sys.Date()), "\n",
      paste(R.Version()$version.string), "\n",
      "\n",
      "### Input Files ###\n",
      "\n",
      paste("Input Genotype File:", genotype_file_name), "\n",
      "\n",
      "### User Selected Parameters ###\n",
      "\n",
      paste("Selected Ploidy:", selected_ploidy), "\n",
      paste("MAF Filter:", input$filter_maf), "\n",
      paste("Min Read Depth (Marker per Sample):", input$size_depth), "\n",
      paste("Remove SNPs with >= % missing data:", input$snp_miss), "\n",
      paste("Remove Samples with >= % missing data:", input$sample_miss), "\n",
      paste("Use Updog Filtering Parameters?:", input$use_updog), "\n",
      paste("Max OD (Updog filter):", ifelse(input$use_updog,input$OD_filter, "NA")), "\n",
      paste("Bias Minimum (Updog filter):", ifelse(input$use_updog,input$Bias[1], "NA")), "\n",
      paste("Bias Maximum (Updog filter):", ifelse(input$use_updog,input$Bias[2], "NA")), "\n",
      paste("Max Prop_mis (Updog filter):", ifelse(input$use_updog,input$Prop_mis,"NA")), "\n",
      paste("Minimum maxpostprob (Updog filter):", ifelse(input$use_updog,input$maxpostprob_filter,"NA")), "\n",
      "\n",
      "### R Packages Used ###\n",
      "\n",
      paste("BIGapp:", packageVersion("BIGapp")), "\n",
      paste("BIGr:", packageVersion("BIGr")), "\n",
      paste("Updog:", packageVersion("updog")), "\n",
      sep = ""
    )
  }

  # Popup for analysis summary
  observeEvent(input$filtering_summary, {
    showModal(modalDialog(
      title = "Summary Information",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close"),
        downloadButton("download_filtering_info", "Download")
      ),
      pre(
        paste(capture.output(filtering_summary_info()), collapse = "\n")
      )
    ))
  })


  # Download Summary Info
  output$download_filtering_info <- downloadHandler(
    filename = function() {
      paste("Filtering_summary_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      # Write the summary info to a file
      writeLines(paste(capture.output(filtering_summary_info()), collapse = "\n"), file)
    }
  )
  
}

## To be copied in the UI
# mod_Filtering_ui("Filtering_1")

## To be copied in the server
# mod_Filtering_server("Filtering_1")
