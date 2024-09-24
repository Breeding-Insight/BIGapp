#' dapc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @importFrom shinycssloaders withSpinner
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_dapc_ui <- function(id){
  ns <- NS(id)
  tagList(
    # Add PCA content here
    fluidRow(
      column(width = 3,
             bs4Dash::box(
               title = "Inputs", width = 12, solidHeader = TRUE, status = "info",
               bs4Dash::tabsetPanel(
                 tabPanel("Step 1:(K)", width = 12,
                          fileInput(ns("dosage_file"), "Choose VCF File", accept = c(".csv",".vcf",".gz")),
                          numericInput(ns("dapc_kmax"), "Maximum K", min = 1, value = 5),
                          numericInput(ns("dapc_ploidy"), "Species Ploidy", min = 1, value = NULL),
                          actionButton(ns("K_start"), "Run Step 1"),
                          div(style="display:inline-block; float:right",dropdownButton(
                            tags$h3("DAPC Inputs"),
                            "You can download an examples of the expected input file here: \n",
                            downloadButton(ns('download_vcf'), "Download VCF Example File"),hr(),
                            actionButton(ns("dapc_summary"), "Summary"),
                            circle = FALSE,
                            status = "warning",
                            icon = icon("info"), width = "300px",
                            tooltip = tooltipOptions(title = "Click to see info!")
                          )),
                 ),
                 tabPanel("Step 2:(DAPC)", width = 12,
                          fileInput(ns("dosage_file"), "Choose VCF File", accept = c(".csv",".vcf",".gz")),
                          numericInput(ns("dapc_k"), "Number of Clusters (K)", min = 1, value = NULL),
                          numericInput(ns("dapc_ploidy"), "Species Ploidy", min = 1, value = NULL),
                          actionButton(ns("dapc_start"), "Run Step 2"),
                          div(style="display:inline-block; float:right",dropdownButton(
                            tags$h3("DAPC Inputs"),
                            "You can download an examples of the expected input file here: \n",
                            downloadButton(ns('download_vcf'), "Download VCF Example File"),
                            circle = FALSE,
                            status = "warning",
                            icon = icon("info"), width = "300px",
                            tooltip = tooltipOptions(title = "Click to see info!")
                          )),
                 )
               )
             ),
             bs4Dash::box(
               title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE,collapsed = FALSE, width = 12,
               "Change the plot parameters", br(),
               selectInput(ns("color_choice"), "Color Palette", choices = c("YlOrRd","YlOrBr","YlGnBu","YlGn",
                                                                            "Reds","RdPu","Purples","PuRd","PuBuGn","PuBu",
                                                                            "OrRd","Oranges","Greys","Greens","GnBu","BuPu",
                                                                            "BuGn","Blues","Set3","Set2","Set1","Pastel2",
                                                                            "Pastel1","Paired","Dark2","Accent","Spectral",
                                                                            "RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn",
                                                                            "PiYG","BrBG"), selected = "Paired"),
               materialSwitch(ns('plot_BICX'), label = "Label Suggested K?", status = "success", value = TRUE),
               div(style="display:inline-block; float:right",
                   dropdownButton(
                     tags$h3("Save Output"),
                     selectInput(inputId = ns('dapc_figure'), label = 'Figure', choices = c("BIC Plot", "DAPC Plot"), selected = "BIC Plot"),
                     selectInput(inputId = ns('dapc_image_type'), label = 'File', choices = c("jpeg","tiff","png"), selected = "jpeg"),
                     sliderInput(inputId = ns('dapc_image_res'), label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                     sliderInput(inputId = ns('dapc_image_width'), label = 'Width', value = 3, min = 1, max = 10, step=0.5),
                     sliderInput(inputId = ns('dapc_image_height'), label = 'Height', value = 3, min = 1, max = 10, step = 0.5),
                     fluidRow(
                       downloadButton(ns("download_dapc_image"), "Save Image"),
                       downloadButton(ns("download_dapc_file"), "Save Files")),
                     circle = FALSE,
                     status = "danger",
                     icon = icon("floppy-disk"), width = "300px",
                     tooltip = tooltipOptions(title = "Click to see inputs!")
                   ))
             ),
             bs4Dash::box(title = "Status", width = 12, collapsible = TRUE, status = "info",
                          progressBar(id = ns("pb_dapc"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
             )
      ),
      column(width = 8,
             bs4Dash::box(title = "DAPC Data", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "info", collapsed = FALSE, maximizable = T,
                          bs4Dash::tabsetPanel(
                            tabPanel("BIC Values",DTOutput(ns('BIC_table'))),
                            tabPanel("DAPC Values", DTOutput(ns('DAPC_table'))), # Placeholder for plot outputs
                            br(), br()
                          )),
             bs4Dash::box(title = "DAPC Plots", status = "info", solidHeader = FALSE, width = 12, height = 550, maximizable = T,
                          bs4Dash::tabsetPanel(
                            tabPanel("BIC Plot",withSpinner(plotOutput(ns("BIC_plot"), height = '460px'))),
                            tabPanel("DAPC Plot", withSpinner(plotOutput(ns("DAPC_plot"), height = '460px'))),
                            tabPanel("STRUCTURE Plot", "Not yet supported"))
             )
      ),
      column(width = 1)
    )
  )
}

#' dapc Server Functions
#' @import grDevices
#' @importFrom graphics axis
#' @importClassesFrom adegenet genlight
#' @importFrom adegenet find.clusters dapc optim.a.score pop<- nInd scatter.dapc
#' @importFrom vcfR read.vcfR extract.gt
#' @importFrom stats BIC as.formula lm logLik median model.matrix na.omit prcomp qbeta quantile runif sd setNames
#' @noRd
mod_dapc_server <- function(input, output, session, parent_session){

  ns <- session$ns


  dapc_items <- reactiveValues(
    grp = NULL,
    bestK = NULL,
    BIC = NULL,
    assignments = NULL,
    dapc = NULL
  )

  ##DAPC analysis
  #Make it a two step process 1) estimate K, and 2) perform DAPC
  observeEvent(input$K_start, {

    toggleClass(id = "dapc_ploidy", class = "borderred", condition = (is.na(input$dapc_ploidy) | is.null(input$dapc_ploidy)))
    if (is.null(input$dosage_file$datapath)) {
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
    req(input$dosage_file$datapath, input$dapc_ploidy)

    ploidy <- as.numeric(input$dapc_ploidy)
    maxK <- as.numeric(input$dapc_kmax)
    geno <- input$dosage_file$datapath

    ##Add in VCF with the vcfR package (input VCF, then convert to genlight using vcf2genlight function)

    #Import genotype information if in VCF format
    vcf <- read.vcfR(geno)

    #Get items in FORMAT column
    info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT

    # Apply the function to the first INFO string
    info_ids <- extract_info_ids(info[1])

    #Get the genotype values if the updog dosage calls are present
    if ("UD" %in% info_ids) {
      genotypeMatrix <- extract.gt(vcf, element = "UD")
      class(genotypeMatrix) <- "numeric"
      rm(vcf) #Remove vcf
    }else{
      #Extract GT and convert to numeric calls
      genotypeMatrix <- extract.gt(vcf, element = "GT")
      genotypeMatrix <- apply(genotypeMatrix, 2, convert_to_dosage)
      rm(vcf) #Remove VCF
    }

    #Perform analysis
    get_k <- findK(genotypeMatrix, maxK, ploidy)

    #Assign results to reactive values
    dapc_items$grp <- get_k$grp
    dapc_items$bestK <- get_k$bestK
    dapc_items$BIC <- get_k$BIC
  })

  observeEvent(input$dapc_start, {

    toggleClass(id = "dapc_ploidy", class = "borderred", condition = (is.na(input$dapc_ploidy) | is.null(input$dapc_ploidy)))
    toggleClass(id = "dapc_k", class = "borderred", condition = (is.na(input$dapc_k) | is.null(input$dapc_k)))

    if (is.null(input$dosage_file$datapath)) {
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
    req(input$dosage_file$datapath, input$dapc_ploidy, input$dapc_k)

    geno <- input$dosage_file$datapath
    ploidy <- as.numeric(input$dapc_ploidy)
    selected_K <- as.numeric(input$dapc_k)

    #Import genotype information if in VCF format
    vcf <- read.vcfR(geno)

    #Get items in FORMAT column
    info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT

    # Apply the function to the first INFO string
    info_ids <- extract_info_ids(info[1])

    #Get the genotype values if the updog dosage calls are present
    if ("UD" %in% info_ids) {
      genotypeMatrix <- extract.gt(vcf, element = "UD")
      class(genotypeMatrix) <- "numeric"
      rm(vcf) #Remove vcf
    }else{
      #Extract GT and convert to numeric calls
      genotypeMatrix <- extract.gt(vcf, element = "GT")
      genotypeMatrix <- apply(genotypeMatrix, 2, convert_to_dosage)
      rm(vcf) #Remove VCF
    }

    #Perform analysis
    clusters <- performDAPC(genotypeMatrix, selected_K, ploidy)

    #Assign results to reactive value
    dapc_items$assignments <- clusters$Q
    dapc_items$dapc <- clusters$dapc
  })

  ###Outputs from DAPC
  #Output the BIC plot
  BIC_plot <- reactive({
    validate(
      need(!is.null(dapc_items$BIC), "Input VCF, define parameters and click `run analysis` in Step 1:(K) to access results in this session.")
    )

    BIC <- dapc_items$BIC
    selected_K <- as.numeric(dapc_items$bestK)
    plot(BIC, type = "o", xaxt = 'n')
    axis(1, at = seq(1, nrow(BIC), 1), labels = TRUE)

    if (input$plot_BICX) {
      plot(BIC, type = "o", xaxt = 'n')
      axis(1, at = seq(1, nrow(BIC), 1), labels = TRUE)
      points(selected_K, BIC[selected_K,2], pch = "x", col = "red", cex = 2)
    } else {
      plot(BIC, type = "o", xaxt = 'n')
      axis(1, at = seq(1, nrow(BIC), 1), labels = TRUE)
    }
  })

  output$BIC_plot <- renderPlot({
    BIC_plot()
  })

  # #Output the DAPC scatter plot
  DAPC_plot <- reactive({
    validate(
      need(!is.null(dapc_items$dapc), "Input VCF, define parameters and click `run analysis` in Step 2:(DAPC) to access results in this session.")
    )

    #Get colors
    palette <- brewer.pal(as.numeric(input$dapc_k), input$color_choice)
    my_palette <- colorRampPalette(palette)(as.numeric(input$dapc_k))

    sc1 <- scatter.dapc(dapc_items$dapc,
                        bg = "white", solid = 1, cex = 1, # cex circle size
                        col = my_palette,
                        pch = 20, # shapes
                        cstar = 1, # 0 or 1, arrows from center of cluster
                        cell = 2, # size of elipse
                        scree.da = T, # plot da
                        scree.pca = T, # plot pca
                        posi.da = "topright",
                        posi.pca="bottomright",
                        mstree = F, # lines connecting clusters
                        lwd = 1, lty = 2,
                        leg = F, clab = 1) # legend and label of legend clusters. clab 0 or 1
  })

  output$DAPC_plot <- renderPlot({
    DAPC_plot()
  })

  # #Output datatables

  BIC_table <-   reactive({
    validate(
      need(!is.null(dapc_items$BIC), "Input VCF, define parameters and click `run analysis` in Step 1:(K) to access results in this session.")
    )
    dapc_items$BIC
  })

  output$BIC_table <- renderDT({
    BIC_table()
  }, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

  assignments_table <- reactive({
    validate(
      need(!is.null(dapc_items$assignments), "Input VCF, define parameters and click `run analysis` in Step 2:(DAPC) to access results in this session.")
    )
    dapc_items$assignments
  })

  output$DAPC_table <- renderDT({
    assignments_table()
  }, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

  #Download figures for DAPC (change this so that the figures were already saved as reactive values to then print() here)
  output$download_dapc_image <- downloadHandler(

    filename = function() {
      if (input$dapc_image_type == "jpeg") {
        paste("dapc-", Sys.Date(), ".jpg", sep="")
      } else if (input$dapc_image_type == "png") {
        paste("dapc-", Sys.Date(), ".png", sep="")
      } else {
        paste("dapc-", Sys.Date(), ".tiff", sep="")
      }
    },
    content = function(file) {
      #req(all_plots$pca_2d, all_plots$pca3d, all_plots$scree, input$pca_image_type, input$pca_image_res, input$pca_image_width, input$pca_image_height) #Get the plots
      req(input$dapc_figure)

      if (input$dapc_image_type == "jpeg") {
        jpeg(file, width = as.numeric(input$dapc_image_width), height = as.numeric(input$dapc_image_height), res= as.numeric(input$dapc_image_res), units = "in")
      } else if (input$dapc_image_type == "png") {
        png(file, width = as.numeric(input$dapc_image_width), height = as.numeric(input$dapc_image_height), res= as.numeric(input$dapc_image_res), units = "in")
      } else {
        tiff(file, width = as.numeric(input$dapc_image_width), height = as.numeric(input$dapc_image_height), res= as.numeric(input$dapc_image_res), units = "in")
      }

      # Conditional plotting based on input selection
      if (input$dapc_figure == "DAPC Plot") {
        req(dapc_items$dapc, input$dapc_k)

        #Get colors
        palette <- brewer.pal(as.numeric(input$dapc_k), input$color_choice)
        my_palette <- colorRampPalette(palette)(as.numeric(input$dapc_k))

        sc1 <- scatter.dapc(dapc_items$dapc,
                            bg = "white", solid = 1, cex = 1, # cex circle size
                            col = my_palette,
                            pch = 20, # shapes
                            cstar = 1, # 0 or 1, arrows from center of cluster
                            cell = 2, # size of elipse
                            scree.da = T, # plot da
                            scree.pca = T, # plot pca
                            posi.da = "topright",
                            posi.pca="bottomright",
                            mstree = F, # lines connecting clusters
                            lwd = 1, lty = 2,
                            leg = F, clab = 1) # legend and label of legend clusters. clab 0 or 1

      } else if (input$dapc_figure == "BIC Plot") {
        req(dapc_items$BIC, dapc_items$bestK)

        BIC <- dapc_items$BIC
        selected_K <- as.numeric(dapc_items$bestK)
        plot(BIC, type = "o", xaxt = 'n')
        axis(1, at = seq(1, nrow(BIC), 1), labels = TRUE)

        if (input$plot_BICX) {
          plot(BIC, type = "o", xaxt = 'n')
          axis(1, at = seq(1, nrow(BIC), 1), labels = TRUE)
          points(selected_K, BIC[selected_K,2], pch = "x", col = "red", cex = 2)
        } else {
          plot(BIC, type = "o", xaxt = 'n')
          axis(1, at = seq(1, nrow(BIC), 1), labels = TRUE)
        }
      }
      dev.off()
    }
  )

  #Download files for DAPC
  output$download_dapc_file <- downloadHandler(
    filename = function() {
      paste0("dapc-results-", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Temporary files list
      temp_dir <- tempdir()
      temp_files <- c()

      if (!is.null(dapc_items$assignments)) {
        # Create a temporary file for assignments
        assignments_file <- file.path(temp_dir, paste0("DAPC-values-", Sys.Date(), ".csv"))
        write.csv(dapc_items$assignments, assignments_file, row.names = TRUE)
        temp_files <- c(temp_files, assignments_file)
      }

      if (!is.null(dapc_items$BIC)) {
        # Create a temporary file for BIC data frame
        bicDF_file <- file.path(temp_dir, paste0("BIC-values-", Sys.Date(), ".csv"))
        write.csv(dapc_items$BIC, bicDF_file, row.names = FALSE)
        temp_files <- c(temp_files, bicDF_file)
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
  
  ##Summary Info
  dapc_summary_info <- function() {
    #Handle possible NULL values for inputs
    dosage_file_name <- if (!is.null(input$dosage_file$name)) input$dosage_file$name else "No file selected"
    #passport_file_name <- if (!is.null(input$passport_file$name)) input$passport_file$name else "No file selected"
    selected_ploidy <- if (!is.null(input$dapc_ploidy)) as.character(input$dapc_ploidy) else "Not selected"

    #Print the summary information
    cat(
      "BIGapp DAPC Summary\n",
      "\n",
      paste0("Date: ", Sys.Date()), "\n",
      paste("R Version:", R.Version()$version.string), "\n",
      "\n",
      "### Input Files ###\n",
      "\n",
      paste("Input Genotype File:", dosage_file_name), "\n",
      #paste("Input Passport File:", passport_file_name), "\n",
      "\n",
      "### User Selected Parameters ###\n",
      "\n",
      paste("Selected Ploidy:", selected_ploidy), "\n",
      paste("Maximum K:", input$dapc_kmax), "\n",
      paste("Number of Clusters (K):", input$dapc_k), "\n",
      "\n",
      "### R Packages Used ###\n",
      "\n",
      paste("BIGapp:", packageVersion("BIGapp")), "\n",
      paste("adegenet:", packageVersion("adegenet")), "\n",
      paste("ggplot2:", packageVersion("ggplot2")), "\n",
      paste("vcfR:", packageVersion("vcfR")), "\n",
      paste("RColorBrewer:", packageVersion("RColorBrewer")), "\n",
      sep = ""
    )
  }

  # Popup for analysis summary
  observeEvent(input$dapc_summary, {
    showModal(modalDialog(
      title = "Summary Information",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close"),
        downloadButton("download_dapc_info", "Download")
      ),
      pre(
        paste(capture.output(dapc_summary_info()), collapse = "\n")
      )
    ))
  })


  # Download Summary Info
  output$download_dapc_info <- downloadHandler(
    filename = function() {
      paste("dapc_summary_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      # Write the summary info to a file
      writeLines(paste(capture.output(dapc_summary_info()), collapse = "\n"), file)
    }
  )
}

## To be copied in the UI
# mod_dapc_ui("dapc_1")

## To be copied in the server
# mod_dapc_server("dapc_1")
