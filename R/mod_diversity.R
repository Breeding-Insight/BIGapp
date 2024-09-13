#' diversity UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_diversity_ui <- function(id){
  ns <- NS(id)
  tagList(
    # Add GWAS content here
    fluidRow(
      column(width = 3,
             box(title="Inputs", width = 12, collapsible = TRUE, collapsed = FALSE, status = "info", solidHeader = TRUE,
                 fileInput(ns("diversity_file"), "Choose VCF File", accept = c(".csv",".vcf",".gz")),
                 numericInput(ns("diversity_ploidy"), "Species Ploidy", min = 1, value = NULL),
                 actionButton(ns("diversity_start"), "Run Analysis"),
                 div(style="display:inline-block; float:right",dropdownButton(
                   tags$h3("Diversity Parameters"),
                   "You can download an examples of the expected input file here: \n",
                   downloadButton(ns('download_vcf'), "Download VCF Example File"),
                   circle = FALSE,
                   status = "warning",
                   icon = icon("info"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see info!")
                 ))
             ),
             box(title = "Plot Controls", width=12, status = "warning", solidHeader = TRUE, collapsible = TRUE,
                 sliderInput(ns("hist_bins"),"Histogram Bins", min = 1, max = 200, value = c(20), step = 1),
                 div(style="display:inline-block; float:left",dropdownButton(
                   tags$h3("Save Image"),
                   selectInput(inputId = ns('div_figure'), label = 'Figure', choices = c("Dosage Plot",
                                                                                         "MAF Histogram",
                                                                                         "OHet Histogram",
                                                                                         "Marker Plot")),
                   selectInput(inputId = ns('div_image_type'), label = 'File Type', choices = c("jpeg","pdf","tiff","png"), selected = "jpeg"),
                   sliderInput(inputId = ns('div_image_res'), label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                   sliderInput(inputId = ns('div_image_width'), label = 'Width', value = 8, min = 1, max = 20, step=0.5),
                   sliderInput(inputId = ns('div_image_height'), label = 'Height', value = 5, min = 1, max = 20, step = 0.5),
                   fluidRow(
                     downloadButton(ns("download_div_figure"), "Save Image"),
                     downloadButton(ns("download_div_file"), "Save Files")),
                   circle = FALSE,
                   status = "danger",
                   icon = icon("floppy-disk"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see inputs!")
                 ))
             )
      ),
      column(width = 6,
             box(
               title = "Plots", status = "info", solidHeader = FALSE, width = 12, height = 550,
               bs4Dash::tabsetPanel(
                 tabPanel("Dosage Plot", plotOutput(ns('dosage_plot')),style = "overflow-y: auto; height: 500px"),
                 tabPanel("MAF Plot", plotOutput(ns('maf_plot')),style = "overflow-y: auto; height: 500px"),
                 tabPanel("OHet Plot", plotOutput(ns('het_plot')),style = "overflow-y: auto; height: 500px"),
                 tabPanel("Marker Plot", plotOutput(ns('marker_plot')),style = "overflow-y: auto; height: 500px"), #Can this be an interactive plotly?
                 tabPanel("Sample Table", DTOutput(ns('sample_table')),style = "overflow-y: auto; height: 470px"),
                 tabPanel("SNP Table", DTOutput(ns('snp_table')),style = "overflow-y: auto; height: 470px")
               )
             )
      ),
      column(width = 3,
             valueBoxOutput(ns("mean_het_box"), width = NULL),
             valueBoxOutput(ns("mean_maf_box"), width = NULL),
             box(title = "Status", width = 12, collapsible = TRUE, status = "info",
                 progressBar(id = ns("pb_diversity"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
             )
      )
    )
  )
}

#' diversity Server Functions
#'
#' @importFrom graphics axis hist points
#' @import ggplot2
#' @importFrom scales comma_format
#'
#' @noRd
mod_diversity_server <- function(input, output, session, parent_session){

    ns <- session$ns
    #######Genomic Diversity analysis

    #Genomic Diversity output files
    diversity_items <- reactiveValues(
      diversity_df = NULL,
      dosage_df = NULL,
      het_df = NULL,
      maf_df = NULL,
      pos_df = NULL,
      markerPlot = NULL
    )

    #Reactive boxes
    output$mean_het_box <- renderValueBox({
      valueBox(
        value = 0,
        subtitle = "Mean Heterozygosity",
        icon = icon("dna"),
        color = "info"
      )
    })

    output$mean_maf_box <- renderValueBox({
      valueBox(
        value = 0,
        subtitle = "Mean MAF",
        icon = icon("dna"),
        color = "info"
      )
    })

    observeEvent(input$diversity_start, {
      toggleClass(id = "diversity_ploidy", class = "borderred", condition = (is.na(input$diversity_ploidy) | is.null(input$diversity_ploidy)))
      #toggleClass(id = "zero_value", class = "borderred", condition = (is.na(input$zero_value) | is.null(input$zero_value)))

      if (is.null(input$diversity_file$datapath)) {
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
      req(input$diversity_file, input$diversity_ploidy)

      #Input variables (need to add support for VCF file)
      ploidy <- as.numeric(input$diversity_ploidy)
      geno <- input$diversity_file$datapath

      #Status
      updateProgressBar(session = session, id = "pb_diversity", value = 20, title = "Importing VCF")

      #Import genotype information if in VCF format
      vcf <- read.vcfR(geno)

      #Save position information
      diversity_items$pos_df <- data.frame(vcf@fix[, 1:2])

      #Get items in FORMAT column
      info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT

      # Apply the function to the first INFO string
      info_ids <- extract_info_ids(info[1])

      #Status
      updateProgressBar(session = session, id = "pb_diversity", value = 40, title = "Converting to Numeric")

      #Get the genotype values and convert to numeric format
      #Extract GT and convert to numeric calls
      geno_mat <- extract.gt(vcf, element = "GT")
      geno_mat <- apply(geno_mat, 2, convert_to_dosage)
      rm(vcf) #Remove VCF
    

      print(class(geno_mat))
      #Convert genotypes to alternate counts if they are the reference allele counts
      #Importantly, the dosage plot is based on the input format NOT the converted genotypes
      is_reference <- TRUE #(input$zero_value == "Reference Allele Counts")

      print("Genotype file successfully imported")
      ######Get MAF plot (Need to remember that the VCF genotypes are likely set as 0 = homozygous reference, where the dosage report is 0 = homozygous alternate)

      print("Starting percentage calc")
      #Status
      updateProgressBar(session = session, id = "pb_diversity", value = 70, title = "Calculating...")
      # Calculate percentages for both genotype matrices
      percentages1 <- calculate_percentages(geno_mat, ploidy)
      # Combine the data matrices into a single data frame
      percentages1_df <- as.data.frame(t(percentages1))
      percentages1_df$Data <- "Dosages"
      # Assuming my_data is your dataframe
      print("Percentage Complete: melting dataframe")
      melted_data <- percentages1_df %>%
        pivot_longer(cols = -(Data),names_to = "Dosage", values_to = "Percentage")

      diversity_items$dosage_df <- melted_data

      print("Dosage calculations worked")

      #Convert the genotype calls prior to het,af, and maf calculation
      geno_mat <- data.frame(convert_genotype_counts(df = geno_mat, ploidy = ploidy, is_reference),
                             check.names = FALSE)

      # Calculating heterozygosity
      diversity_items$het_df <- calculate_heterozygosity(geno_mat, ploidy = ploidy)

      print("Heterozygosity success")
      diversity_items$maf_df <- calculateMAF(geno_mat, ploidy = ploidy)
      diversity_items$maf_df <- diversity_items$maf_df[, c(1,3)]

      print("MAF success")

      #Updating value boxes
      output$mean_het_box <- renderValueBox({
        valueBox(
          value = round(mean(diversity_items$het_df$ObservedHeterozygosity),3),
          subtitle = "Mean Heterozygosity",
          icon = icon("dna"),
          color = "info"
        )
      })
      output$mean_maf_box <- renderValueBox({
        valueBox(
          value = round(mean(diversity_items$maf_df$MAF),3),
          subtitle = "Mean MAF",
          icon = icon("dna"),
          color = "info"
        )
      })

      #Status
      updateProgressBar(session = session, id = "pb_diversity", value = 100, title = "Complete!")
    })

    box_plot <- reactive({
      validate(
        need(!is.null(diversity_items$dosage_df), "Input VCF, define parameters and click `run analysis` to access results in this session.")
      )

      #Plotting
      box <- ggplot(diversity_items$dosage_df, aes(x=Dosage, y=Percentage, fill=Data)) +
        #geom_point(aes(color = Data), position = position_dodge(width = 0.8), width = 0.2, alpha = 0.5) +  # Add jittered points
        geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.9) +
        labs(x = "\nDosage", y = "Percentage\n", title = "Genotype Distribution by Sample") +
        theme_bw() +
        theme(
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14)
        )

      box
    })

    output$dosage_plot <- renderPlot({
      box_plot()
    })

    #Het plot
    het_plot <- reactive({
      validate(
        need(!is.null(diversity_items$het_df) & !is.null(input$hist_bins), "Input VCF, define parameters and click `run analysis` to access results in this session.")
      )
      hist(diversity_items$het_df$ObservedHeterozygosity, breaks = as.numeric(input$hist_bins), col = "tan3", border = "black", xlim= c(0,1),
           xlab = "Observed Heterozygosity",
           ylab = "Number of Samples",
           main = "Sample Observed Heterozygosity")

      axis(1, at = seq(0, 1, by = 0.1), labels = TRUE)
    })

    output$het_plot <- renderPlot({
      het_plot()
    })

    #AF Plot
    #af_plot <- reactive({
    #  validate(
    #    need(!is.null(diversity_items$maf_df) & !is.null(input$hist_bins), "Input VCF, define parameters and click `run analysis` to access results in this session.")
    #  )
    #  hist(diversity_items$maf_df$AF, breaks = as.numeric(input$hist_bins), col = "grey", border = "black", xlab = "Alternate Allele Frequency",
    #       ylab = "Frequency", main = "Alternate Allele Frequency Distribution")
    #})

    #output$af_plot <- renderPlot({
    #  af_plot()
    #})

    #MAF plot
    maf_plot <- reactive({
      validate(
        need(!is.null(diversity_items$maf_df) & !is.null(input$hist_bins), "Input VCF, define parameters and click `run analysis` to access results in this session.")
      )

      hist(diversity_items$maf_df$MAF, breaks = as.numeric(input$hist_bins), col = "grey", border = "black", xlab = "Minor Allele Frequency (MAF)",
           ylab = "Frequency", main = "Minor Allele Frequency Distribution")
    })

    #Marker plot
    marker_plot <- reactive({
      validate(
        need(!is.null(diversity_items$pos_df), "Input VCF, define parameters and click `run analysis` to access results in this session.")
      )
      #Order the Chr column
      diversity_items$pos_df$POS <- as.numeric(diversity_items$pos_df$POS)
      # Sort the dataframe and pad with a 0 if only a single digit is provided
      diversity_items$pos_df$CHROM <- ifelse(
        nchar(diversity_items$pos_df$CHROM) == 1, 
        paste0("0", diversity_items$pos_df$CHROM), 
        diversity_items$pos_df$CHROM
      )
      diversity_items$pos_df <- diversity_items$pos_df[order(diversity_items$pos_df$CHROM), ]

      #Plot

      # Create custom breaks for the x-axis labels (every 13Mb)
      x_breaks <- seq(0, max(diversity_items$pos_df$POS), by = (max(diversity_items$pos_df$POS)/5))
      x_breaks <- c(x_breaks, max(diversity_items$pos_df$POS))  # Add 114Mb as a custom break

      # Create custom labels for the x-axis using the 'Mb' suffix
      x_labels <- comma_format()(x_breaks / 1000000)
      x_labels <- paste0(x_labels, "Mb")

      suppressWarnings({
        markerPlot <- ggplot(diversity_items$pos_df, aes(x = as.numeric(POS), y = CHROM, group = as.factor(CHROM))) +
          geom_point(aes(color = as.factor(CHROM)), shape = 108, size = 5, show.legend = FALSE) +
          xlab("Position") +
          #ylab("Markers\n") +
          theme(axis.text = element_text(size = 11, color = "black"),
                axis.text.x.top = element_text(size = 11, color = "black"),
                axis.title = element_blank(),
                panel.grid = element_blank(),
                axis.ticks.length.x = unit(-0.15, "cm"),
                axis.ticks.margin = unit(0.1, "cm"),
                axis.ticks.y = element_blank(),
                axis.line.x.top = element_line(color="black"),
                panel.background = element_rect(fill="white"),
                plot.margin = margin(10, 25, 10, 10)
          ) +
          scale_x_continuous(
            breaks = x_breaks,     # Set custom breaks for x-axis labels
            labels = x_labels,     # Set custom labels with "Mb" suffixes
            position = "top",       # Move x-axis labels and ticks to the top
            expand = c(0,0),
            limits = c(0,max(diversity_items$pos_df$POS))
          )
      })
      #Display plot
      markerPlot
    })

    output$marker_plot <- renderPlot({
      marker_plot()
    })

    output$maf_plot <- renderPlot({
      maf_plot()
    })

    sample_table <- reactive({
      validate(
        need(!is.null(diversity_items$het_df), "Input VCF, define parameters and click `run analysis` to access results in this session.")
      )
      diversity_items$het_df
    })

    output$sample_table <- renderDT({sample_table()}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

    snp_table <- reactive({
      validate(
        need(!is.null(diversity_items$maf_df), "Input VCF, define parameters and click `run analysis` to access results in this session.")
      )
      diversity_items$maf_df
    })

    output$snp_table <- renderDT({snp_table()}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

    #Download Figures for Diversity Tab (Need to convert figures to ggplot)
    output$download_div_figure <- downloadHandler(

      filename = function() {
        if (input$div_image_type == "jpeg") {
          paste("genomic-diversity-", Sys.Date(), ".jpg", sep="")
        } else if (input$div_image_type == "png") {
          paste("genomic-diversity-", Sys.Date(), ".png", sep="")
        } else {
          paste("genomic-diversity-", Sys.Date(), ".tiff", sep="")
        }
      },
      content = function(file) {
        req(input$div_figure)

        if (input$div_image_type == "jpeg") {
          jpeg(file, width = as.numeric(input$div_image_width), height = as.numeric(input$div_image_height), res= as.numeric(input$div_image_res), units = "in")
        } else if (input$div_image_type == "png") {
          png(file, width = as.numeric(input$div_image_width), height = as.numeric(input$div_image_height), res= as.numeric(input$div_image_res), units = "in")
        } else {
          tiff(file, width = as.numeric(input$div_image_width), height = as.numeric(input$div_image_height), res= as.numeric(input$div_image_res), units = "in")
        }

        # Conditional plotting based on input selection
        if (input$div_figure == "Dosage Plot") {
          print(box_plot())
        } else if (input$div_figure == "AF Histogram") {
          af_plot()
        } else if (input$div_figure == "MAF Histogram") {
          maf_plot()
        } else if (input$div_figure == "OHet Histogram") {
          het_plot()
        } else if (input$div_figure == "Marker Plot") {
          print(marker_plot())
        }

        dev.off()
      }

    )

    #Download files for Genotype Diversity
    output$download_div_file <- downloadHandler(
      filename = function() {
        paste0("genomic-diversity-results-", Sys.Date(), ".zip")
      },
      content = function(file) {
        # Temporary files list
        temp_dir <- tempdir()
        temp_files <- c()

        if (!is.null(diversity_items$het_df)) {
          # Create a temporary file for assignments
          het_file <- file.path(temp_dir, paste0("Sample-statistics-", Sys.Date(), ".csv"))
          write.csv(diversity_items$het_df, het_file, row.names = FALSE)
          temp_files <- c(temp_files, het_file)
        }

        if (!is.null(diversity_items$maf_df)) {
          # Create a temporary file for BIC data frame
          maf_file <- file.path(temp_dir, paste0("SNP-statistics-", Sys.Date(), ".csv"))
          write.csv(diversity_items$maf_df, maf_file, row.names = FALSE)
          temp_files <- c(temp_files, maf_file)
        }

        # Zip files only if there's something to zip
        if (length(temp_files) > 0) {
          zip(file, files = temp_files, extras = "-j") # Using -j to junk paths
        }

        # Optionally clean up
        file.remove(temp_files)
      }
    )

    output$download_vcf <- downloadHandler(
      filename = function() {
        paste0("BIGapp_VCF_Example_file.vcf.gz")
      },
      content = function(file) {
        ex <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")
        file.copy(ex, file)
      })
}

## To be copied in the UI
# mod_diversity_ui("diversity_1")

## To be copied in the server
# mod_diversity_server("diversity_1")