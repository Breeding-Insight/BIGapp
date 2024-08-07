#' PCA UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom AGHmatrix Gmatrix
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyWidgets virtualSelectInput progressBar
#' @importFrom shiny NS tagList
#' @importFrom shinyWidgets materialSwitch
#' @importFrom shinyjs inlineCSS toggleClass
#'
mod_PCA_ui <- function(id){
  ns <- NS(id)
  tagList(
    # Add PCA content here
    fluidRow(
      useShinyjs(),
      inlineCSS(list(.borderred = "border-color: red", .redback = "background: red")),

      column(width = 3,
             box(
               title = "Inputs", width = 12, solidHeader = TRUE, status = "info",
               p("* Required"),
               fileInput(ns("dosage_file"), "Choose Genotypes File*", accept = c(".csv",".vcf",".gz")),
               fileInput(ns("passport_file"), "Choose Passport File (Sample IDs in first column)", accept = c(".csv")),
               #Dropdown will update after passport upload
               numericInput(ns("pca_ploidy"), "Species Ploidy*", min = 2, value = NULL),
               actionButton(ns("pca_start"), "Run Analysis"),
               div(style="display:inline-block; float:right",
                   dropdownButton(
                     tags$h3("PCA Inputs"),
                     "PCA Input file and analysis info",
                     circle = FALSE,
                     status = "warning",
                     icon = icon("info"), width = "300px",
                     tooltip = tooltipOptions(title = "Click to see info!")
                   )),
               style = "overflow-y: auto; height: 420px"
             ),
             box(
               title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE,collapsed = FALSE, width = 12,
               "Change the PCA plot parameters", br(),
               selectInput(ns('group_info'), label = 'Color Variable (eg, Taxon)', choices = NULL),
               materialSwitch(ns('use_cat'), label = "Color Specific Category Within Variable?", status = "success"),
               conditionalPanel(condition = "input.use_cat",
                                ns = ns,
                                virtualSelectInput(
                                  inputId = ns("cat_color"),
                                  label = "Select Category To Color:",
                                  choices = NULL,
                                  showValueAsTags = TRUE,
                                  search = TRUE,
                                  multiple = TRUE
                                ),
                                selectInput(ns("grey_choice"), "Select Grey", choices = c("Light Grey", "Grey", "Dark Grey", "Black"), selected = "Grey")
               ),
               selectInput(ns("color_choice"), "Color Palette", choices = list("Standard Palettes" = c("Set1","Set3","Pastel2",
                                                                                                      "Pastel1","Accent","Spectral",
                                                                                                      "RdYlGn","RdGy"),
                                                                            "Colorblind Friendly" = c("Set2","Paired","Dark2","YlOrRd","YlOrBr","YlGnBu","YlGn",
                                                                                                       "Reds","RdPu","Purples","PuRd","PuBuGn","PuBu",
                                                                                                       "OrRd","Oranges","Greys","Greens","GnBu","BuPu",
                                                                                                       "BuGn","Blues","RdYlBu",
                                                                                                       "RdBu", "PuOr","PRGn","PiYG","BrBG"
                                                                            )), 
                                                                            selected = "Set1"),
               selectInput(ns("pc_X"), "X-Axis (2D-Plot only)", choices = c("PC1","PC2","PC3","PC4","PC5"), selected = "PC1"),
               selectInput(ns("pc_Y"), "Y-Axis (2D-Plot only)", choices = c("PC1","PC2","PC3","PC4","PC5"), selected = "PC2"),
               div(style="display:inline-block; float:right",dropdownButton(
                 tags$h3("Save Image"),
                 selectInput(inputId = ns('pca_figure'), label = 'Figure', choices = c("2D Plot", "Scree Plot"), selected = "2D Plot"),
                 selectInput(inputId = ns('pca_image_type'), label = 'File', choices = c("jpeg","tiff","png"), selected = "jpeg"),
                 sliderInput(inputId = ns('pca_image_res'), label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                 sliderInput(inputId = ns('pca_image_width'), label = 'Width', value = 10, min = 1, max = 20, step=0.5),
                 sliderInput(inputId = ns('pca_image_height'), label = 'Height', value = 6, min = 1, max = 20, step = 0.5),
                 downloadButton(ns("download_pca"), "Save"),
                 circle = FALSE,
                 status = "danger",
                 icon = icon("floppy-disk"), width = "300px",
                 tooltip = tooltipOptions(title = "Click to see inputs!")
               ))
             )
      ),
      column(width = 8,
             box(title = "Passport Data", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "info", collapsed = FALSE,
                 DTOutput(ns('passport_table')),
                 style = "overflow-y: auto; height: 420px"
             ),
             box(
               title = "PCA Plots", status = "info", solidHeader = FALSE, width = 12, height = 550,
               bs4Dash::tabsetPanel(
                 tabPanel("3D-Plot",withSpinner(plotlyOutput(ns("pca_plot"), height = '460px'))),
                 tabPanel("2D-Plot", withSpinner(plotOutput(ns("pca_plot_ggplot"), height = '460px'))),
                 tabPanel("Scree Plot", withSpinner(plotOutput(ns("scree_plot"), height = '460px')))) # Placeholder for plot outputs
             )
      ),
      column(width = 1)
    )
  )
}

#' PCA Server Functions
#'
#' @importFrom factoextra get_eigenvalue
#' @import grDevices
#' @importFrom shinyWidgets updateVirtualSelect
#' @importFrom plotly layout plotlyOutput renderPlotly add_markers plot_ly
#' @importFrom RColorBrewer brewer.pal
#' @importFrom shinyjs toggleClass
#'
#' @noRd
mod_PCA_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    #PCA reactive values
    pca_data <- reactiveValues(
      pc_df_pop = NULL,
      variance_explained = NULL,
      my_palette = NULL
    )

    # Update dropdown menu choices based on uploaded passport file
    passport_table <- reactive({
      validate(
        need(!is.null(input$passport_file), "Upload passport file to access results in this section."),
      )
      info_df <- read.csv(input$passport_file$datapath, header = TRUE, check.names = FALSE)
      info_df[,1] <- as.character(info_df[,1]) #Makes sure that the sample names are characters instead of numeric

      updateSelectInput(session, "group_info", choices = colnames(info_df))
      info_df
    })

    output$passport_table <- renderDT({
      passport_table()},
      options = list(scrollX = TRUE,
                     autoWidth = FALSE,
                     pageLength = 4))

    #PCA specific category selection
    observeEvent(input$group_info, {
      #updateMaterialSwitch(session, inputId = "use_cat", status = "success")

      # Get selected column name
      selected_col <- input$group_info

      # Extract unique values from the selected column
      unique_values <- unique(passport_table()[[selected_col]])

      #Add category selection
      updateVirtualSelect("cat_color", choices = unique_values, session = session)

    })

    #PCA events
    observeEvent(input$pca_start, {

      # Missing input with red border and alerts
      toggleClass(id = "pca_ploidy", class = "borderred", condition = is.na(input$pca_ploidy))
      if (is.null(input$dosage_file)) {
        shinyalert(
          title = "Missing input!",
          text = "Upload Genotypes File",
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
      req(input$pca_ploidy, input$dosage_file$datapath)

      # Get inputs
      geno <- input$dosage_file$datapath
      g_info <- as.character(input$group_info)
      output_name <- input$output_name
      ploidy <- input$pca_ploidy

      #Notification
      showNotification("PCA analysis in progress...")

      #Import genotype info if genotype matrix format
      if (grepl("\\.csv$", geno)) {
        genomat <- read.csv(geno, header = TRUE, row.names = 1, check.names = FALSE)
      } else{

        #Import genotype information if in VCF format
        vcf <- read.vcfR(geno)

        #Get items in FORMAT column
        info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT

        # Apply the function to the first INFO string
        info_ids <- extract_info_ids(info[1])

        #Get the genotype values if the updog dosage calls are present
        if ("UD" %in% info_ids) {
          genomat <- extract.gt(vcf, element = "UD")
          class(genomat) <- "numeric"
          rm(vcf) #Remove vcf
        }else{
          #Extract GT and convert to numeric calls
          genomat <- extract.gt(vcf, element = "GT")
          genomat <- apply(genomat, 2, convert_to_dosage)
          rm(vcf) #Remove VCF
        }
      }

      #Start analysis

      # Passport info
      if (!is.null(input$passport_file$datapath) && input$passport_file$datapath != "") {
        info_df <- read.csv(input$passport_file$datapath, header = TRUE, check.names = FALSE)

        # Check for duplicates in the first column
        duplicated_samples <- info_df[duplicated(info_df[, 1]), 1]
        if (length(duplicated_samples) > 0) {
          shinyalert(
            title = "Duplicate Samples Detected in Passport File",
            text = paste("The following samples are duplicated:", paste(unique(duplicated_samples), collapse = ", ")),
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
          req(length(duplicated_samples) == 0) # Stop the analysis if duplicates are found
        }

      } else {
        info_df <- data.frame(SampleID = colnames(genomat))
      }

      # Print the modified dataframe
      row.names(info_df) <- info_df[,1]

      #Plotting
      #First build a relationship matrix using the genotype values
      G.mat.updog <- Gmatrix(t(genomat), method = "VanRaden", ploidy = as.numeric(ploidy), missingValue = "NA")

      #PCA
      prin_comp <- prcomp(G.mat.updog, scale = TRUE)
      eig <- get_eigenvalue(prin_comp)
      round(sum(eig$variance.percent[1:3]),1)

      ###Simple plots
      # Extract the PC scores
      pc_scores <- prin_comp$x

      # Create a data frame with PC scores
      pc_df <- data.frame(PC1 = pc_scores[, 1], PC2 = pc_scores[, 2],
                          PC3 = pc_scores[, 3], PC4 = pc_scores[, 4],
                          PC5 = pc_scores[, 5], PC6 = pc_scores[, 6],
                          PC7 = pc_scores[, 7], PC8 = pc_scores[, 8],
                          PC9 = pc_scores[, 9], PC10 = pc_scores[, 10])


      # Compute the percentage of variance explained for each PC
      variance_explained <- round(100 * prin_comp$sdev^2 / sum(prin_comp$sdev^2), 1)


      # Retain only samples in common
      row.names(info_df) <- info_df[,1]
      info_df <- info_df[row.names(pc_df),]

      #Add the information for each sample
      pc_df_pop <- merge(pc_df, info_df, by.x = "row.names", by.y = "row.names", all.x = TRUE)


      # Ignore color input if none is entered by user
      if (g_info != "") {
        pc_df_pop[[g_info]] <- as.factor(pc_df_pop[[g_info]])
      } else {
        g_info <- NULL
      }

      #Update global variable
      pca_dataframes <- pc_df_pop

      # Generate a distinct color palette if g_info is provided
      if (!is.null(g_info) && g_info != "") {
        unique_countries <- unique(pc_df_pop[[g_info]])
        palette <- brewer.pal(length(unique_countries), input$color_choice)
        my_palette <- colorRampPalette(palette)(length(unique_countries))
      } else {
        unique_countries <- NULL
        my_palette <- NULL
      }

      # Store processed data in reactive values
      pca_data$pc_df_pop <- pc_df_pop
      pca_data$variance_explained <- variance_explained
      pca_data$my_palette <- my_palette

      #End of PCA section
    })

    ##2D PCA plotting
    pca_2d <- reactive({
      validate(
        need(!is.null(pca_data$pc_df_pop), "Input Genotype file, Species ploidy, and run the analysis to access results in this section.")
      )

      # Generate colors
      if (!is.null(pca_data$my_palette)) {
        unique_countries <- unique(pca_data$pc_df_pop[[input$group_info]])
        palette <- brewer.pal(length(unique_countries), input$color_choice)
        my_palette <- colorRampPalette(palette)(length(unique_countries))
      } else {
        unique_countries <- NULL
        my_palette <- NULL
      }

      # Define a named vector to map input labels to grey values
      label_to_value <- c("Light Grey" = "grey80",
                          "Grey" = "grey60",
                          "Dark Grey" = "grey40",
                          "Black" = "black")

      # Get the corresponding value based on the selected grey
      selected_grey <- label_to_value[[input$grey_choice]]

      #Set factor
      if (!input$use_cat && is.null(my_palette)) {
        print("No Color Info")
      }else{
        pca_data$pc_df_pop[[input$group_info]] <- as.factor(pca_data$pc_df_pop[[input$group_info]])
      }

      # Similar plotting logic here

      cat_colors <- c(input$cat_color, "grey")
      plot <- {if(!is.null(input$group_info) & input$group_info != "")
        ggplot(pca_data$pc_df_pop, aes(x = pca_data$pc_df_pop[[input$pc_X]],
                                       y = pca_data$pc_df_pop[[input$pc_Y]],
                                       color = factor(pca_data$pc_df_pop[[input$group_info]]))) else
                                         ggplot(pca_data$pc_df_pop, aes(x = pca_data$pc_df_pop[[input$pc_X]],
                                                                        y = pca_data$pc_df_pop[[input$pc_Y]]))} +
        geom_point(size = 2, alpha = 0.8) +
        {if(input$use_cat) scale_color_manual(values = setNames(c(my_palette, "grey"), cat_colors), na.value = selected_grey) else
          if(!is.null(my_palette)) scale_color_manual(values = my_palette)} +
        guides(color = guide_legend(override.aes = list(size = 5.5), nrow = 17)) +
        theme_minimal() +
        theme(
          panel.border = element_rect(color = "black", fill = NA),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 16)
        ) +
        labs(
          x = paste0(input$pc_X, "(", pca_data$variance_explained[as.numeric(substr(input$pc_X, 3, 3))], "%)"),
          y = paste0(input$pc_Y, "(", pca_data$variance_explained[as.numeric(substr(input$pc_Y, 3, 3))], "%)"),
          color = input$group_info
        )

      plot  # Assign the plot to your reactiveValues
    })

    #Plot the 2d plot
    output$pca_plot_ggplot <- renderPlot({
      pca_2d()
    })

    #3D PCA plotting
    pca_plot <- reactive({
      #Plotly
      validate(
        need(!is.null(pca_data$pc_df_pop), "Input Genotype file, Species ploidy, and run the analysis to access results in this section.")
      )

      #Generate colors
      unique_countries <- unique(pca_data$pc_df_pop[[input$group_info]])
      palette <- brewer.pal(length(unique_countries),input$color_choice)
      my_palette <- colorRampPalette(palette)(length(unique_countries))

      tit = paste0('Total Explained Variance =', sum(pca_data$variance_explained[1:3]))

      fig <- plot_ly(pca_data$pc_df_pop, x = ~PC1, y = ~PC2, z = ~PC3, color = pca_data$pc_df_pop[[input$group_info]],
                     colors = my_palette) %>%
        add_markers(size = 12, text = paste0("Sample:",pca_data$pc_df_pop$Row.names))

      fig <- fig %>%
        layout(
          title = tit,
          scene = list(bgcolor = "white")
        )

      fig # Return the Plotly object here
    })

    output$pca_plot <- renderPlotly({
      pca_plot()
    })

    pca_scree <- reactive({
      #PCA scree plot
      validate(
        need(!is.null(pca_data$variance_explained), "Input Genotype file, Species ploidy, and run the analysis to access the results in this section.")
      )

      var_explained <- pca_data$variance_explained

      # Create a data frame for plotting
      plot_data <- data.frame(PC = 1:10, Variance_Explained = var_explained[1:10])

      # Use ggplot for plotting
      plot <- ggplot(plot_data, aes(x = PC, y = Variance_Explained)) +
        geom_bar(stat = "identity", fill = "lightblue", alpha = 0.9, color = "black") +  # Bars with some transparency
        geom_line(color = "black") +  # Connect points with a line
        geom_point(color = "black") +  # Add points on top of the line for emphasis
        scale_x_continuous(breaks = 1:10, limits = c(0.5, 10.5)) +
        xlab("Principal Component") +
        ylab("% Variance Explained") +
        ylim(0, 100) +
        theme_bw() +
        theme(
          panel.border = element_rect(color = "black", fill = NA),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 16)
        )
      plot
    })

    #Scree plot
    output$scree_plot <- renderPlot({
      pca_scree()
    })

    #Download figures for PCA
    output$download_pca <- downloadHandler(

      filename = function() {
        if (input$pca_image_type == "jpeg") {
          paste("pca-", Sys.Date(), ".jpg", sep="")
        } else if (input$pca_image_type == "png") {
          paste("pca-", Sys.Date(), ".png", sep="")
        } else {
          paste("pca-", Sys.Date(), ".tiff", sep="")
        }
      },
      content = function(file) {
        req(input$pca_figure)

        if (input$pca_image_type == "jpeg") {
          jpeg(file, width = as.numeric(input$pca_image_width), height = as.numeric(input$pca_image_height), res= as.numeric(input$pca_image_res), units = "in")
        } else if (input$pca_image_type == "png") {
          png(file, width = as.numeric(input$pca_image_width), height = as.numeric(input$pca_image_height), res= as.numeric(input$pca_image_res), units = "in")
        } else {
          tiff(file, width = as.numeric(input$pca_image_width), height = as.numeric(input$pca_image_height), res= as.numeric(input$pca_image_res), units = "in")
        }

        # Conditional plotting based on input selection
        if (input$pca_figure == "2D Plot") {
          pca_2d()
        } else if (input$pca_figure == "Scree Plot") {
          pca_scree()
        } else {
          plot(x = 1:10, y = 1:10, main = "Fallback Simple Test Plot")  # Fallback simple test plot
        }
        dev.off()
      }
    )
  })
}

## To be copied in the UI
# mod_PCA_ui("PCA_1")

## To be copied in the server
# mod_PCA_server("PCA_1")
