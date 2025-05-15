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
    fluidPage(
      # Add PCA content here
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
        useShinyjs(),
        inlineCSS(list(.borderred = "border-color: red", .redback = "background: red")),

        column(width = 3,
               box(
                 title = "Inputs", width = 12, solidHeader = TRUE, status = "info",
                 p("* Required"),
                 fileInput(ns("dosage_file"), "Choose VCF File*", accept = c(".csv",".vcf",".gz")),
                 fileInput(ns("passport_file"), "Choose Trait File", accept = c(".csv")),
                 #Dropdown will update after passport upload
                 numericInput(ns("pca_ploidy"), "Species Ploidy*", min = 2, value = NULL),
                 br(),
                 actionButton(ns("pca_start"), "Run Analysis"),
                 div(style="display:inline-block; float:right",dropdownButton(
                   HTML("<b>Input files</b>"),
                   p(downloadButton(ns('download_vcf'),""), "VCF Example File"),
                   p(downloadButton(ns('download_pheno'),""), "Trait Example File"), hr(),
                   p(HTML("<b>Parameters description:</b>"), actionButton(ns("goPar"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>Results description:</b>"), actionButton(ns("goRes"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>How to cite:</b>"), actionButton(ns("goCite"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   actionButton(ns("pca_summary"), "Summary"),
                   circle = FALSE,
                   status = "warning",
                   icon = icon("info"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see info!")
                 ))#,
                 #style = "overflow-y: auto; height: 480px"
               ),
               box(
                 title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE,collapsed = FALSE, width = 12,
                 selectInput(ns('group_info'), label = 'Color Variable', choices = ""),
                 materialSwitch(ns('use_shapes'), label = "Add Shapes (2D-Plot only)", status = "success"),
                 materialSwitch(ns('use_cat'), label = "Color Category (2D-Plot only)", status = "success"),
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
                   selectInput(inputId = ns('pca_image_type'), label = 'File', choices = c("jpeg","tiff","png","svg"), selected = "jpeg"),
                   sliderInput(inputId = ns('pca_image_res'), label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                   sliderInput(inputId = ns('pca_image_width'), label = 'Width', value = 10, min = 1, max = 20, step=0.5),
                   sliderInput(inputId = ns('pca_image_height'), label = 'Height', value = 6, min = 1, max = 20, step = 0.5),
                   downloadButton(ns("download_pca"), "Save Image"),
                   circle = FALSE,
                   status = "danger",
                   icon = icon("floppy-disk"), width = "300px", label = "Save",
                   tooltip = tooltipOptions(title = "Click to see inputs!")
                 ))
               )
        ),
        column(width = 8,
               box(
                 title = "PCA Plots", status = "info", solidHeader = FALSE, width = 12, height = 650, maximizable = T,
                 bs4Dash::tabsetPanel(
                   tabPanel("3D-Plot",withSpinner(plotlyOutput(ns("pca_plot"), height = '560px'))),
                   tabPanel("2D-Plot", withSpinner(plotOutput(ns("pca_plot_ggplot"), height = '560px'))),
                   tabPanel("Scree Plot", withSpinner(plotOutput(ns("scree_plot"), height = '560px')))) # Placeholder for plot outputs
               )
        ),
        column(width = 1)
      )
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
#' @import shinydisconnect
#'
#' @noRd
mod_PCA_server <- function(input, output, session, parent_session){

  ns <- session$ns
  #options(warn = -1) #Uncomment to suppress warnings

  # Help links
  observeEvent(input$goPar, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "PCA_tabset",
                      selected = "PCA_par")
    # expand specific box
    updateBox(id = "PCA_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goRes, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "PCA_tabset",
                      selected = "PCA_results")
    # expand specific box
    updateBox(id = "PCA_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goCite, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "PCA_tabset",
                      selected = "PCA_cite")
    # expand specific box
    updateBox(id = "PCA_box", action = "toggle", session = parent_session)
  })
  
  #Default choices
  trait_options <- reactiveValues(
    missing_data = "NA",
    custom_missing = NULL,
    sample_column = NULL,
    file_type = NULL
  )
  
  # Update dropdown menu choices based on uploaded passport file
  passport_table <- reactive({
    validate(
      need(!is.null(input$passport_file), "Upload passport file to access results in this section."),
    )
    info_df <- read.csv(input$passport_file$datapath, header = TRUE, check.names = FALSE, na.strings=trait_options$missing_data)
    info_df[,1] <- as.character(info_df[,1]) #Makes sure that the sample names are characters instead of numeric
    
    updateSelectInput(session, "group_info", choices = colnames(info_df))
    info_df
  })
  
  #PCA specific category selection
  observeEvent(passport_table(), {
    #updateMaterialSwitch(session, inputId = "use_cat", status = "success")
    
    # Get selected column name
    selected_col <- input$group_info
    
    # Extract unique values from the selected column
    unique_values <- unique(passport_table()[[selected_col]])
    
    #Add category selection
    updateVirtualSelect("cat_color", choices = unique_values, session = session)
    
  })
  
  # PCA specific category selection. Update when group_info is changed.
  observeEvent(input$group_info, {
    req(passport_table()) # Ensure passport_table is available
    selected_col <- input$group_info
    unique_values <- unique(passport_table()[[selected_col]])
    updateVirtualSelect("cat_color", choices = unique_values, session = session)
  })
  
  #UI popup window for input
  observeEvent(input$passport_file, {
    req(input$passport_file)
    #Get the column names of the csv file
    info_df <- read.csv(input$passport_file$datapath, header = TRUE, check.names = FALSE, nrows=2)
    info_df[,1] <- as.character(info_df[,1]) #Makes sure that the sample names are characters instead of numeric
    
    # Read first 5 rows for preview
    preview_data <- tryCatch({
      head(read.csv(input$passport_file$datapath, nrows = 5, na.strings=trait_options$missing_data),5)
    }, error = function(e) {
      NULL
    })
    
    showModal(modalDialog(
      title = "Trait File Options",
      size= "l",
      
      selectInput(
        inputId = ns('missing_data'),
        label = 'Missing Data Value',
        choices = c("NA",".","-99","(blank)","Custom"),
        selected = trait_options$missing_data  # Initialize with stored value
      ),
      conditionalPanel(
        condition = "input.missing_data == 'Custom'", ns = ns,
        div(
          textInput(
            inputId = ns('custom_missing'),
            label = 'Custom Missing Value',
            value = trait_options$custom_missing  # Initialize with stored value
          )
        ),
        div(
          id = ns("custom_missing_warning"),
          style = "color: red;",
          textOutput(ns("custom_missing_msg"))
        )
      ),
      selectInput(
        inputId = ns('sample_column'),
        label = 'Sample ID Column',
        choices = colnames(info_df)
      ),
      
      if (!is.null(preview_data)) {
        div(
          h4(
            "File Preview (First 5 Rows)",
            style = "font-size: 18px; color: darkgrey;" # Smaller and purple
          ),
          div(
            style = "background-color: #f0f0f0; padding: 10px; border: 1px solid #ccc;", # Grey box style
            div(
              style = "max-width: 100%; overflow-x: auto;", # Constrain table width and enable horizontal scrolling
              tableOutput(ns("file_preview"))
            )
          )
        )
      } else {
        div(
          p("Could not load file preview.")
        )
      },
      
      footer = tagList(
        actionButton(ns("save_trait_options"), "Save")
      )
    ))
    
    # Render the preview table
    output$file_preview <- renderTable({
      req(preview_data)
      preview_data
    })
    
  })
  
  output$custom_missing_msg <- renderText({
    if (input$missing_data == "Custom" && nchar(input$custom_missing) == 0) {
      "Please enter a custom missing value."
    } else {
      ""
    }
  })
  
  
  #Close popup window when user "saves options"
  observeEvent(input$save_trait_options, {
    trait_options$missing_data <- input$missing_data
    trait_options$custom_missing <- input$custom_missing
    trait_options$sample_column <- input$sample_column
    #trait_options$file_type
    # Save other inputs as needed
    
    if (input$missing_data == "Custom" && nchar(input$custom_missing) == 0) {
      # Validation failed: display warning and prevent modal closure
      showNotification(
        "Please enter a custom missing value.",
        type = "error",
        duration = NULL # Make it persistent
      )
      return() # Stop further execution and keep the modal open
    }
    
    removeModal()  # Close the modal after saving
  })

  #PCA reactive values
  pca_data <- reactiveValues(
    pc_df_pop = NULL,
    variance_explained = NULL,
    my_palette = NULL
  )

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
      vcf <- read.vcfR(geno, verbose = FALSE)

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

    #Start analysis following a genotype data check
    
    if (ncol(genomat) < 5){
      shinyalert(
        title = "Small Genotype File",
        text = "BIGapp currently requires at least 5 samples to perform a PCA",
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
      req(ncol(genomat) > 5) # Stop the analysis if duplicates are found
    }

    # Passport info
    if (!is.null(input$passport_file$datapath) && input$passport_file$datapath != "") {
      #Importing trait file based on user inputs
      if (trait_options$missing_data == "(blank)") {
        info_df <- read.csv(input$passport_file$datapath, header = TRUE, check.names = FALSE, na.strings="")
      } else if (trait_options$missing_data == "Custom") {
        info_df <- read.csv(input$passport_file$datapath, header = TRUE, check.names = FALSE, na.strings = trait_options$custom_missing)
      } else {
        info_df <- read.csv(input$passport_file$datapath, header = TRUE, check.names = FALSE, na.strings = trait_options$missing_data)
      }
      
      # Make the sample ID column the first column in the dataframe
      sample_col_name <- input$sample_column
      info_df <- info_df[, c(sample_col_name, setdiff(names(info_df), sample_col_name))]

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

    # Check ploidy
    if(input$pca_ploidy != max(genomat, na.rm = T)){
      shinyalert(
        title = "Wrong ploidy",
        text = "The input ploidy does not match the maximum dosage found in the genotype file",
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
      req(input$pca_ploidy == max(genomat, na.rm = T))
    }

    #Plotting
    #First build a relationship matrix using the genotype values
    G.mat.updog <- Gmatrix(t(genomat), method = "VanRaden", ploidy = as.numeric(ploidy), missingValue = "NA")

    #PCA
    prin_comp <- prcomp(G.mat.updog, scale = TRUE)
    eig <- get_eigenvalue(prin_comp)

    ###Simple plots
    # Extract the PC scores
    pc_scores <- prin_comp$x

    # Create a data frame with PC scores
    pc_df <- data.frame(PC1 = pc_scores[, 1], PC2 = pc_scores[, 2],
                        PC3 = pc_scores[, 3], PC4 = pc_scores[, 4],
                        PC5 = pc_scores[, 5])


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
    } else {
      if(input$group_info != "") pca_data$pc_df_pop[[input$group_info]] <- as.factor(pca_data$pc_df_pop[[input$group_info]])
    }

    # Similar plotting logic here

    #input$cat_color <- as.character(unique(pca_data$pc_df_pop[[input$group_info]]))
    cat_colors <- c(input$cat_color, "grey")
    
    #Add a check if too many values are needing shapes
    if (length(my_palette) > 25 & input$use_shapes) {
      shinyalert(
        title = "Too Many Values!",
        text = "There are more than 25 unique values selected for the shape option",
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
      
      #Trying to update the ui switch back to "off", but it is not working
      updateMaterialSwitch(session = parent_session, inputId = "use_shapes", value = NULL)
    }
    
    shapes_list <- c(1:25)
    
    #plot
    plot <- ggplot(pca_data$pc_df_pop, aes(
      x = pca_data$pc_df_pop[[input$pc_X]],
      y = pca_data$pc_df_pop[[input$pc_Y]],
      color = if (!is.null(input$group_info) & input$group_info != "") factor(pca_data$pc_df_pop[[input$group_info]]) else NULL,
      shape = if (!is.null(input$group_info) & input$group_info != "" & input$use_shapes) factor(pca_data$pc_df_pop[[input$group_info]]) else NULL
    )) +
      geom_point(size = 2, alpha = 0.8) +
      
      #add color
      {if (input$use_cat & !is.null(my_palette)) 
        scale_color_manual(values = setNames(c(my_palette, "grey"), cat_colors), na.value = selected_grey)
        else if (!is.null(my_palette)) 
          scale_color_manual(values = my_palette)} +
      
      #add shapes
      {if (input$use_shapes)
        scale_shape_manual(values = shapes_list)} +  # Use up to 25 shapes
      
      #adjust legends
      guides(
        color = guide_legend(override.aes = list(size = 5.5), nrow = 17),
        shape = if (input$use_shapes) guide_legend(override.aes = list(size = 5.5), nrow = 17) else NULL
      ) +
      
      #add labels
      labs(
        x = paste0(input$pc_X, " (", pca_data$variance_explained[as.numeric(substr(input$pc_X, 3, 3))], "%)"),
        y = paste0(input$pc_Y, " (", pca_data$variance_explained[as.numeric(substr(input$pc_Y, 3, 3))], "%)"),
        color = if (!is.null(input$group_info)) input$group_info else NULL,
        shape = if (input$use_shapes & !is.null(input$group_info)) input$group_info else NULL
      ) +
      
      #add minimal theme
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16)
      )
    
    plot
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

    tit = paste0('Total Explained Variance =', sum(pca_data$variance_explained[1:3]))

    #Generate colors
    if(input$group_info!= ""){
      unique_countries <- unique(pca_data$pc_df_pop[[input$group_info]])
      palette <- brewer.pal(length(unique_countries),input$color_choice)
      my_palette <- colorRampPalette(palette)(length(unique_countries))

      fig <- plot_ly(pca_data$pc_df_pop, x = ~PC1, y = ~PC2, z = ~PC3, color = as.factor(pca_data$pc_df_pop[[input$group_info]]),
                     colors = my_palette) %>%
        add_markers(size = 12, text = paste0("Sample:",pca_data$pc_df_pop$Row.names))

    } else {
      fig <- plot_ly(pca_data$pc_df_pop, x = ~PC1, y = ~PC2, z = ~PC3) %>%
        add_markers(size = 12, text = paste0("Sample:",pca_data$pc_df_pop$Row.names))
    }

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

  ##Summary Info
  pca_summary_info <- function() {
    # Handle possible NULL values for inputs
    dosage_file_name <- if (!is.null(input$dosage_file$name)) input$dosage_file$name else "No file selected"
    passport_file_name <- if (!is.null(input$passport_file$name)) input$passport_file$name else "No file selected"
    selected_ploidy <- if (!is.null(input$pca_ploidy)) as.character(input$pca_ploidy) else "Not selected"

    # Print the summary information
    cat(
      "BIGapp PCA Summary\n",
      "\n",
      paste0("Date: ", Sys.Date()), "\n",
      paste("R Version:", R.Version()$version.string), "\n",
      "\n",
      "### Input Files ###\n",
      "\n",
      paste("Input Genotype File:", dosage_file_name), "\n",
      paste("Input Passport File:", passport_file_name), "\n",
      "\n",
      "### User Selected Parameters ###\n",
      "\n",
      paste("Selected Ploidy:", selected_ploidy), "\n",
      "\n",
      "### R Packages Used ###\n",
      "\n",
      paste("BIGapp:", packageVersion("BIGapp")), "\n",
      paste("AGHmatrix:", packageVersion("AGHmatrix")), "\n",
      paste("ggplot2:", packageVersion("ggplot2")), "\n",
      paste("plotly:", packageVersion("plotly")), "\n",
      paste("factoextra:", packageVersion("factoextra")), "\n",
      paste("RColorBrewer:", packageVersion("RColorBrewer")), "\n",
      sep = ""
    )
  }

  # Popup for analysis summary
  observeEvent(input$pca_summary, {
    showModal(modalDialog(
      title = "Summary Information",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close"),
        downloadButton("download_pca_info", "Download")
      ),
      pre(
        paste(capture.output(pca_summary_info()), collapse = "\n")
      )
    ))
  })


  # Download Summary Info
  output$download_pca_info <- downloadHandler(
    filename = function() {
      paste("PCA_summary_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      # Write the summary info to a file
      writeLines(paste(capture.output(pca_summary_info()), collapse = "\n"), file)
    }
  )


  #Download figures for PCA
  output$download_pca <- downloadHandler(
    filename = function() {
      if (input$pca_image_type == "jpeg") {
        paste("pca-", Sys.Date(), ".jpg", sep = "")
      } else if (input$pca_image_type == "png") {
        paste("pca-", Sys.Date(), ".png", sep = "")
      } else if (input$pca_image_type == "tiff") {
        paste("pca-", Sys.Date(), ".tiff", sep = "")
      } else {
        paste("pca-", Sys.Date(), ".svg", sep = "")
      }
    },
    content = function(file) {
      req(input$pca_figure)

      if (input$pca_image_type == "jpeg") {
        jpeg(file, width = as.numeric(input$pca_image_width), height = as.numeric(input$pca_image_height), res = as.numeric(input$pca_image_res), units = "in")
      } else if (input$pca_image_type == "png") {
        png(file, width = as.numeric(input$pca_image_width), height = as.numeric(input$pca_image_height), res = as.numeric(input$pca_image_res), units = "in")
      } else if (input$pca_image_type == "tiff") {
        tiff(file, width = as.numeric(input$pca_image_width), height = as.numeric(input$pca_image_height), res = as.numeric(input$pca_image_res), units = "in")
      } else{
        svg(file, width = as.numeric(input$pca_image_width), height = as.numeric(input$pca_image_height))
      }

      # Plot based on user selection
      if (input$pca_figure == "2D Plot") {
        print(pca_2d())
      } else if (input$pca_figure == "Scree Plot") {
        print(pca_scree())
      }

      dev.off()
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

  output$download_pheno <- downloadHandler(
    filename = function() {
      paste0("BIGapp_passport_Example_file.csv")
    },
    content = function(file) {
      ex <- system.file("iris_passport_file.csv", package = "BIGapp")
      file.copy(ex, file)
    })

}

## To be copied in the UI
# mod_PCA_ui("PCA_1")

## To be copied in the server
# mod_PCA_server("PCA_1")
