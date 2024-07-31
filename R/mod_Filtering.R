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
#' @import dplyr
#'
#'
mod_Filtering_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 3,
             box(width = 12,
                 title = "Quality Filtering", status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                 fileInput(ns("updog_rdata"),"Choose VCF File", accept = c(".vcf",".gz")),
                 textInput(ns("filter_output_name"), "Output File Name"),
                 numericInput(ns("filter_ploidy"),"Ploidy", min = 0, value = NULL),
                 numericInput(ns("filter_maf"),"MAF filter", min = 0, max=1, value = 0.05, step = 0.01),
                 sliderInput(ns("size_depth"),"Min Read Depth (Marker per Sample)", min = 0, max = 300, value = 10, step = 1),
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
                 downloadButton(ns("start_updog_filter"), "Download Filtered VCF", icon = icon("download")),
                 div(style="display:inline-block; float:right",dropdownButton(
                   tags$h3("Updog Filter Parameters"),
                   "Add description of each filter. Presently, all filtering parameters that are typically used for processing
                    a VCF file from Updog dosage calling are included. If a VCF file does not contain these values, it will only be
                    filtered for read depth, missing data, and maf.",
                   circle = FALSE,
                   status = "warning",
                   icon = icon("info"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see info!")
                 ))
             )
      ),
      column(width = 6,
             tabBox(width =12, collapsible = FALSE, status = "info",
                    id = "updog_tab", height = "600px",
                    tabPanel("Bias Histogram", icon = icon("image"), plotOutput(ns("bias_hist"), height = '550px')),
                    tabPanel("OD Histogram", icon = icon("image"), plotOutput(ns("od_hist"), height = '550px')),
                    tabPanel("Prop_mis Histogram", icon = icon("image"), plotOutput(ns("maxpostprob_hist"), height = '550px')),
                    tabPanel("SNP_miss", icon = icon("image"), plotOutput(ns("missing_snp_hist"), height = '550px')),
                    tabPanel("Sample_miss", icon = icon("image"), plotOutput(ns("missing_sample_hist"), height = '550px'))
             )
      ),
      column(width = 3,
             valueBoxOutput(ns("snp_retained_box"), width = NULL),
             valueBoxOutput(ns("snp_removed_box"), width = NULL),
             box(title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE,
                 sliderInput(ns("hist_bins"),"Histogram Bins", min = 1, max = 1200, value = c(50), step = 1), width = NULL,
                 div(style="display:inline-block; float:left",dropdownButton(
                   tags$h3("Save Image"),
                   selectInput(inputId = ns('filter_hist'), label = 'Figure', choices = c("Bias Histogram",
                                                                                          "OD Histogram",
                                                                                          "Prop_mis Histogram",
                                                                                          "SNP_mis",
                                                                                          "Sample_mis")),
                   selectInput(inputId = ns('image_type'), label = 'File Type', choices = c("jpeg","tiff","png"), selected = "jpeg"),
                   sliderInput(inputId = ns('image_res'), label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                   sliderInput(inputId = ns('image_width'), label = 'Width', value = 8, min = 1, max = 20, step=0.5),
                   sliderInput(inputId = ns('image_height'), label = 'Height', value = 5, min = 1, max = 20, step = 0.5),
                   downloadButton(ns("download_filter_hist"), "Save"),
                   circle = FALSE,
                   status = "danger",
                   icon = icon("floppy-disk"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see inputs!")
                 ))
             ),
             box(title = "Status", width =12, collapsible = TRUE, status = "info",
                 progressBar(id = ns("pb_filter"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
             )
      )
    )
  )
}

#' Filtering Server Functions
#'
#' @import vcfR
#' @import BIGr
#' @importFrom graphics abline axis hist
#'
#' @noRd
mod_Filtering_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    #vcf
    filtering_files <- reactiveValues(
      raw_vcf_df = NULL,
      sample_miss_df = NULL,
      snp_miss_df = NULL

    )

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

    #Updog filtering
    output$start_updog_filter <- downloadHandler(
      filename = function() {
        paste0(input$filter_output_name, ".vcf.gz")
      },
      content = function(file) {
        # Ensure the files are uploaded
        req(input$filter_ploidy, input$filter_output_name,input$updog_rdata)

        if (input$use_updog) {
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

        temp_file <- tempfile(fileext = ".vcf.gz")

        updateProgressBar(session = session, id = "pb_filter", value = 10, title = "Processing VCF file")
        #Input file
        vcf <- read.vcfR(input$updog_rdata$datapath, verbose = FALSE)
        #Starting SNPs
        starting_snps <- nrow(vcf)
        #export INFO dataframe
        filtering_files$raw_vcf_df <- data.frame(vcf@fix)

        #Pb
        updateProgressBar(session = session, id = "pb_filter", value = 40, title = "Filtering VCF file")

        #Filtering
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

        #Getting missing data information
        #Add support for genotype matrix filtering?
        #Pb
        updateProgressBar(session = session, id = "pb_filter", value = 50, title = "Calculating Missing Data")

        gt_matrix <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
        filtering_files$snp_miss_df <- rowMeans(is.na(gt_matrix)) #SNP missing values
        filtering_files$sample_miss_df <- as.numeric(colMeans(is.na(gt_matrix))) #Sample missing values
        rm(gt_matrix) #Remove gt matrix

        #Pb
        updateProgressBar(session = session, id = "pb_filter", value = 80, title = "Exporting Filtered VCF")

        #Writing file
        write.vcf(vcf, file = temp_file)

        #Get final_snps
        final_snps <- nrow(vcf)

        # Check if the VCF file was created
        if (file.exists(temp_file)) {
          cat("VCF file created successfully.\n")

          # Move the file to the path specified by 'file'
          file.copy(temp_file, file, overwrite = TRUE)

          # Delete the temporary file
          unlink(temp_file)
        } else {
          stop("Error: Failed to create the VCF file.")
        }

        # Status
        updateProgressBar(session = session, id = "pb_filter", value = 100, title = "Finished!")

        #Updating value boxes
        output$snp_retained_box <- renderValueBox({
          valueBox(
            value = final_snps,
            subtitle = "SNPs Retained",
            icon = icon("dna"),
            color = "info"
          )
        })
        output$snp_removed_box <- renderValueBox({
          valueBox(
            value = round(((starting_snps - final_snps)/starting_snps*100),1),
            subtitle = "Percent SNPs Removed",
            icon = icon("dna"),
            color = "info"
          )
        })

        #Unload vcf
        rm(vcf)
      }
    )

    #Download figures for VCF Filtering
    output$download_filter_hist <- downloadHandler(

      filename = function() {
        if (input$image_type == "jpeg") {
          paste("VCF-histogram-", Sys.Date(), ".jpg", sep="")
        } else if (input$image_type == "png") {
          paste("VCF-histogram-", Sys.Date(), ".png", sep="")
        } else {
          paste("VCF-histogram-", Sys.Date(), ".tiff", sep="")
        }
      },
      content = function(file) {
        req(input$image_type)

        if (input$image_type == "jpeg") {
          jpeg(file, width = as.numeric(input$image_width), height = as.numeric(input$image_height), res= as.numeric(input$image_res), units = "in")
        } else if (input$image_type == "png") {
          png(file, width = as.numeric(input$image_width), height = as.numeric(input$image_height), res= as.numeric(input$image_res), units = "in")
        } else {
          tiff(file, width = as.numeric(input$image_width), height = as.numeric(input$image_height), res= as.numeric(input$image_res), units = "in")
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

        } else if (input$filter_hist == "SNP_mis") {

          hist(as.numeric(filtering_files$snp_miss_df),
               main = "Ratio of Missing Data per SNP After Filtering",
               xlab = "Proportion of Missing Data per SNP",
               ylab = "Number of SNPs",
               col = "lightblue",
               border = "black",
               xlim = c(0,1),
               breaks = as.numeric(input$hist_bins))
          axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

          # Add vertical lines
          abline(v = mean(as.numeric(filtering_files$snp_miss_df)), col = "red", lty = 2)  # Mean line
          abline(v = median(as.numeric(filtering_files$snp_miss_df)), col = "green", lty = 2)  # Median line
          abline(v = quantile(as.numeric(filtering_files$snp_miss_df), 0.95), col = "blue", lty = 2)

        } else if (input$filter_hist == "Sample_mis") {

          hist(as.numeric(filtering_files$sample_miss_df),
               main = "Ratio of Missing Data per Sample After Filtering",
               xlab = "Proportion of Missing Data per Sample",
               ylab = "Number of Samples",
               col = "lightblue",
               border = "black",
               xlim = c(0,1),
               breaks = as.numeric(input$hist_bins))
          axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

          # Add vertical lines
          abline(v = mean(as.numeric(filtering_files$sample_miss_df)), col = "red", lty = 2)  # Mean line
          abline(v = median(as.numeric(filtering_files$sample_miss_df)), col = "green", lty = 2)  # Median line
          abline(v = quantile(as.numeric(filtering_files$sample_miss_df), 0.95), col = "blue", lty = 2)
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


      # Function to split INFO column and expand it into multiple columns
      split_info_column <- function(info) {
        # Split the INFO column by semicolon
        info_split <- str_split(info, ";")[[1]]

        # Create a named list by splitting each element by equals sign
        info_list <- set_names(map(info_split, ~ str_split(.x, "=")[[1]][2]),
                               map(info_split, ~ str_split(.x, "=")[[1]][1]))

        return(info_list)
      }

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
      })

      ###OD
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

      })

      ##MAXPOSTPROB

      #Histogram

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

      })

      #Missing data
      output$missing_snp_hist <- renderPlot({

        #Histogram
        hist(as.numeric(filtering_files$snp_miss_df),
             main = "Ratio of Missing Data per SNP After Filtering",
             xlab = "Proportion of Missing Data per SNP",
             ylab = "Number of SNPs",
             col = "lightblue",
             border = "black",
             xlim = c(0,1),
             breaks = as.numeric(input$hist_bins))
        axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

        # Add vertical lines
        abline(v = mean(as.numeric(filtering_files$snp_miss_df)), col = "red", lty = 2)  # Mean line
        abline(v = median(as.numeric(filtering_files$snp_miss_df)), col = "green", lty = 2)  # Median line
        abline(v = quantile(as.numeric(filtering_files$snp_miss_df), 0.95), col = "blue", lty = 2)

      })

      output$missing_sample_hist <- renderPlot({

        #Histogram
        hist(as.numeric(filtering_files$sample_miss_df),
             main = "Ratio of Missing Data per Sample After Filtering",
             xlab = "Proportion of Missing Data per Sample",
             ylab = "Number of Samples",
             col = "lightblue",
             border = "black",
             xlim = c(0,1),
             breaks = as.numeric(input$hist_bins))
        axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

        # Add vertical lines
        abline(v = mean(as.numeric(filtering_files$sample_miss_df)), col = "red", lty = 2)  # Mean line
        abline(v = median(as.numeric(filtering_files$sample_miss_df)), col = "green", lty = 2)  # Median line
        abline(v = quantile(as.numeric(filtering_files$sample_miss_df), 0.95), col = "blue", lty = 2)

      })

      ##Read Depth (I would prefer that this show the mean depth for SNPs or Samples instead of all loci/sample cells)
      quantile(as.numeric(new_df$DP), 0.95)
    })
  })
}

## To be copied in the UI
# mod_Filtering_ui("Filtering_1")

## To be copied in the server
# mod_Filtering_server("Filtering_1")
