#' dapc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#'
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
                          fileInput(ns("dosage_file"), "Choose Genotypes File", accept = c(".csv",".vcf",".vcf.gz")),
                          #fileInput("passport_file", "Choose Passport File (Sample IDs in first column)", accept = c(".csv")),
                          #textInput("output_name", "Output File Name (disabled)"),
                          #Dropdown will update after pasport upload
                          numericInput(ns("dapc_kmax"), "Maximum K", min = 1, value = 5),
                          numericInput(ns("dapc_ploidy"), "Species Ploidy", min = 1, value = 2),
                          actionButton(ns("K_start"), "Run Analysis"),
                          div(style="display:inline-block; float:right",dropdownButton(
                            tags$h3("DAPC Inputs"),
                            "DAPC Input file and analysis info. The DAPC analysis is broken down into two steps. The first step (Step 1), uses Kmeans clustering to estimate the most likely number of clusters within the dataset. This is visualized in the BIC plot and is typically the minimum BIC value. Step 2 is the DAPC analysis where the most likely value for K (number of clusters) is input and the cluster memberships are determined in the DAPC results",
                            circle = FALSE,
                            status = "warning",
                            icon = icon("info"), width = "300px",
                            tooltip = tooltipOptions(title = "Click to see info!")
                          )),
                          #style = "overflow-y: auto; height: 420px"
                 ),
                 tabPanel("Step 2:(DAPC)", width = 12,
                          fileInput(ns("dosage_file"), "Choose Genotypes File", accept = c(".csv",".vcf",".vcf.gz")),
                          #fileInput("passport_file", "Choose Passport File (Sample IDs in first column)", accept = c(".csv")),
                          #textInput("output_name", "Output File Name (disabled)"),
                          #Dropdown will update after pasport upload
                          numericInput(ns("dapc_k"), "Number of Clusters (K)", min = 1, value = NULL),
                          numericInput(ns("dapc_ploidy"), "Species Ploidy", min = 1, value = 2),
                          actionButton(ns("dapc_start"), "Run Analysis"),
                          div(style="display:inline-block; float:right",dropdownButton(
                            tags$h3("DAPC Inputs"),
                            "DAPC Input file and analysis info",
                            circle = FALSE,
                            status = "warning",
                            icon = icon("info"), width = "300px",
                            tooltip = tooltipOptions(title = "Click to see info!")
                          )),
                          #style = "overflow-y: auto; height: 420px"
                 )
               )
               #style = "overflow-y: auto; height: 420px"
             ),
             bs4Dash::box(
               title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE,collapsed = FALSE, width = 12,
               "Change the plot parameters", br(),
               #selectInput('group_info', label = 'Color Variable (eg, Taxon)', choices = NULL),
               selectInput(ns("color_choice"), "Color Palette", choices = c("YlOrRd","YlOrBr","YlGnBu","YlGn",
                                                                            "Reds","RdPu","Purples","PuRd","PuBuGn","PuBu",
                                                                            "OrRd","Oranges","Greys","Greens","GnBu","BuPu",
                                                                            "BuGn","Blues","Set3","Set2","Set1","Pastel2",
                                                                            "Pastel1","Paired","Dark2","Accent","Spectral",
                                                                            "RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn",
                                                                            "PiYG","BrBG"), selected = "Paired"),
               materialSwitch(ns('plot_BICX'), label = "Label Suggested K?", status = "success", value = TRUE),
               #selectInput("pc_X", "X-Axis (2D-Plot only)", choices = c("PC1","PC2","PC3","PC4","PC5"), selected = "PC1"),
               #selectInput("pc_Y", "Y-Axis (2D-Plot only)", choices = c("PC1","PC2","PC3","PC4","PC5"), selected = "PC2"),
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
             bs4Dash::box(title = "DAPC Data", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "info", collapsed = FALSE,
                          bs4Dash::tabsetPanel(
                            tabPanel("BIC Values",DTOutput(ns('BIC_table'))),
                            tabPanel("DAPC Values", DTOutput(ns('DAPC_table'))) # Placeholder for plot outputs
                            #style = "overflow-y: auto; height: 420px"
                          )),
             bs4Dash::box(title = "DAPC Plots", status = "info", solidHeader = FALSE, width = 12, height = 550,
                          bs4Dash::tabsetPanel(
                            tabPanel("BIC Plot",shinycssloaders::withSpinner(plotOutput(ns("BIC_plot"), height = '460px'))),
                            tabPanel("DAPC Plot", shinycssloaders::withSpinner(plotOutput(ns("DAPC_plot"), height = '460px'))),
                            tabPanel("STRUCTURE Plot", "Not yet supported"))
                          #tabPanel("STRUCTURE Plot", plotOutput("STRUCTURE_plot", height = '460px')))
             )
      ),
      column(width = 1)
    )
  )
}

#' dapc Server Functions
#'
#' @noRd
mod_dapc_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    ##DAPC analysis
    #Make it a two step process 1) estimate K, and 2) perform DAPC

    dapc_items <- reactiveValues(
      grp = NULL,
      bestK = NULL,
      bicDF = NULL,
      assignments = NULL,
      dapc = NULL

    )

    observeEvent(input$K_start, {

      ploidy <- as.numeric(input$dapc_ploidy)
      maxK <- as.numeric(input$dapc_kmax)
      geno <- input$dosage_file$datapath

      ##Add in VCF with the vcfR package (input VCF, then convert to genlight using vcf2genlight function)

      #Import genotype data
      #genotypeMatrix <- read.csv(geno, header = TRUE, row.names = 1, check.names = FALSE)
      #Import genotype information if in VCF format
      vcf <- read.vcfR(geno)

      convert_to_dosage <- function(gt) {
        # Split the genotype string
        alleles <- strsplit(gt, "[|/]")
        # Sum the alleles, treating NA values appropriately
        sapply(alleles, function(x) {
          if (any(is.na(x))) {
            return(NA)
          } else {
            return(sum(as.numeric(x), na.rm = TRUE))
          }
        })
      }

      #Get items in FORMAT column
      info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT
      extract_info_ids <- function(info_string) {
        # Split the INFO string by ';'
        info_parts <- strsplit(info_string, ":")[[1]]
        # Extract the part before the '=' in each segment
        info_ids <- gsub("=.*", "", info_parts)
        return(info_ids)
      }

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

      findK <- function(genotypeMatrix, maxK, ploidy) {
        # Convert the genotype matrix to a genlight object
        genlight_new <- new("genlight", t(genotypeMatrix),
                            ind.names = row.names(t(genotypeMatrix)),
                            loc.names = colnames(t(genotypeMatrix)),
                            ploidy = ploidy,
                            NA.char = NA)

        #Assign the populations as the sample names since there is no assumed populations
        pop(genlight_new) <- genlight_new@ind.names

        #Estimate number of clusters
        #Retain all pca for the find.clusters step. Retain as few as possible while maximizing variance captured for DAPC step.
        #Choose is the option to allow adegenet to select the best cluster number based on the BIC minimum
        #The default criterion is "diffNgroup, which is not necessarily the minimum BIC, but based on the sharp decrease of the BIC value.
        #Either way, this is a suggestion, and the number of clusters to use should be made with biology considerations.
        graphics.off() #Prevent plot from automatically displaying
        grp <- find.clusters(genlight_new, max.n.clust = maxK,
                             n.pca = nInd(genlight_new),
                             stat = "BIC",
                             criterion = "diffNgroup",
                             parallel = FALSE,
                             choose = FALSE)

        # Identify the best K based on lowest BIC
        bestK <- length(grp$size)

        # Create a BIC dataframe
        bicDF <- data.frame(K = 1:maxK, BIC = as.data.frame(grp$Kstat)$`grp$Kstat`)

        return(list(bestK = as.numeric(bestK), grp = grp, BIC = bicDF))

      }

      #Perform analysis
      get_k <- findK(genotypeMatrix, maxK, ploidy)


      #Assign results to reactive values
      dapc_items$grp <- get_k$grp
      dapc_items$bestK <- get_k$bestK
      dapc_items$bicDF <- get_k$BIC

    })

    observeEvent(input$dapc_start, {

      #req()
      geno <- input$dosage_file$datapath
      ploidy <- as.numeric(input$dapc_ploidy)
      selected_K <- as.numeric(input$dapc_k)

      #Import genotype data
      #genotypeMatrix <- read.csv(geno, header = TRUE, row.names = 1, check.names = FALSE)
      #Import genotype information if in VCF format
      vcf <- read.vcfR(geno)

      convert_to_dosage <- function(gt) {
        # Split the genotype string
        alleles <- strsplit(gt, "[|/]")
        # Sum the alleles, treating NA values appropriately
        sapply(alleles, function(x) {
          if (any(is.na(x))) {
            return(NA)
          } else {
            return(sum(as.numeric(x), na.rm = TRUE))
          }
        })
      }

      #Get items in FORMAT column
      info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT
      extract_info_ids <- function(info_string) {
        # Split the INFO string by ';'
        info_parts <- strsplit(info_string, ":")[[1]]
        # Extract the part before the '=' in each segment
        info_ids <- gsub("=.*", "", info_parts)
        return(info_ids)
      }

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

      performDAPC <- function(genotypeMatrix, selected_K, ploidy) {

        #Convert matrix to genlight
        genlight_new <- new("genlight", t(genotypeMatrix),
                            ind.names = row.names(t(genotypeMatrix)),
                            loc.names = colnames(t(genotypeMatrix)),
                            ploidy = ploidy,
                            NA.char = NA)

        #Get groups based on specified cluster number (K)
        graphics.off() #Prevent plot from automatically displaying
        grp <- find.clusters(genlight_new, n.clust = selected_K,
                             n.pca = nInd(genlight_new),
                             stat = "BIC",
                             criterion = "diffNgroup",
                             parallel = FALSE,
                             choose = FALSE)

        # Find the optimal number of principal components
        #NOTE: The default n.da is K-1, but I have read previously to use #Samples - 1?
        dapc1 <- dapc(genlight_new, grp$grp,
                      n.pca = nInd(genlight_new),
                      n.da = nInd(genlight_new)-1,
                      parallel = FALSE)

        #xval <- xvalDapc(genind_obj, n.pca.max = 10, n.da = NULL) #xval does not accept NAs
        a.score <- optim.a.score(dapc1, plot = FALSE)
        n.pca <- a.score$best

        # Perform DAPC with the best K
        finalDapc <- dapc(genlight_new, grp$grp, n.pca = n.pca, n.da = selected_K-1, parallel= FALSE)

        # Extract the membership probabilities
        Q <- as.data.frame(finalDapc$posterior)

        # Add cluster assignments to Q dataframe
        Q$Cluster_Assignment <- finalDapc$assign

        #a data.frame giving the contributions of original variables (alleles in the case of genetic data) to the principal components of DAPC.
        #dapc$var.contr

        # Return list containing BIC dataframe, Q dataframe w/ dapc assignments
        return(list(Q = Q, dapc = finalDapc))
      }

      #Perform analysis
      clusters <- performDAPC(genotypeMatrix, selected_K, ploidy)

      #Assign results to reactive value
      dapc_items$assignments <- clusters$Q
      dapc_items$dapc <- clusters$dapc
    })

    ###Outputs from DAPC
    #Output the BIC plot
    output$BIC_plot <- renderPlot({
      req(dapc_items$bicDF, dapc_items$bestK)

      BIC <- dapc_items$bicDF
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

    #Output the DAPC scatter plot
    output$DAPC_plot <- renderPlot({
      req(dapc_items$dapc, input$dapc_k)

      #Get colors
      palette <- brewer.pal(as.numeric(input$dapc_k), input$color_choice)
      my_palette <- colorRampPalette(palette)(as.numeric(input$dapc_k))

      sc1 <- scatter(dapc_items$dapc,
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
    #Output datatables
    output$BIC_table <- renderDT({dapc_items$bicDF}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))
    output$DAPC_table <- renderDT({dapc_items$assignments}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

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

          sc1 <- scatter(dapc_items$dapc,
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
          req(dapc_items$bicDF, dapc_items$bestK)

          BIC <- dapc_items$bicDF
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
          #} else if (input$pca_figure == "3D Plot") {
          #print(all_plots$pca3d)  # Assuming you might want a 3D plot as well
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

        if (!is.null(dapc_items$bicDF)) {
          # Create a temporary file for BIC data frame
          bicDF_file <- file.path(temp_dir, paste0("BIC-values-", Sys.Date(), ".csv"))
          write.csv(dapc_items$bicDF, bicDF_file, row.names = FALSE)
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
  })
}

## To be copied in the UI
# mod_dapc_ui("dapc_1")

## To be copied in the server
# mod_dapc_server("dapc_1")
