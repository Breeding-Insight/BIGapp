# Install and load required packages
required_packages <- c("updog", "ggplot2", "VariantAnnotation", "SNPRelate",
                       "adegenet", "future", "scales", "AGHmatrix", "stats", 
                       "factoextra", "readxl", "ggrepel", "dplyr", "shiny",
                       "shinydashboard","randomcoloR","plotly", "DT","RColorBrewer",
                       "dichromat", "bs4Dash", "shinyWidgets", "GWASpoly","data.table",
                       "matrixcalc","Matrix", "shinyalert")

for(package in required_packages) {
  if(!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

# UI
ui <- dashboardPage(
  skin = "black",
  dashboardHeader(title = tagList(
    tags$img(src = 'BIG_logo_edit.png', height = '40', width = '40'),
    tags$span("BreedingInsights w/ Genomics", style = "font-size: 12px; margin-left: 1px;")
  )
),
  dashboardSidebar(
    skin="light", status = "info",
    sidebarMenu(
      flat = FALSE,
      menuItem("Home", tabName = "welcome", icon = icon("house")),
      menuItem("Dosage Calling", tabName = "dosage_calling", icon = icon("diagram-next"),
              menuSubItem("Updog Dosage Calling", tabName = "updog", icon = icon("list-ol")),
              menuSubItem("SNP Filtering", tabName = "filtering", icon = icon("filter"))),
      menuItem("PCA", tabName = "pca", icon = icon("chart-simple")),
      menuItem("GWAS", tabName = "gwas", icon = icon("think-peaks")),
      menuItem("Genomic Diversity", tabName = "diversity", icon = icon("chart-pie")),
      menuItem("QTL Analysis", tabName = "qtl", icon = icon("chart-area")),
      menuItem("Genomic Prediction", tabName = "prediction", icon = icon("right-left")),
      menuItem("Source Code", tabName = "code", icon = icon("circle-info")),
      menuItem("Help", tabName = "help", icon = icon("circle-question"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "welcome",
        # Add welcome content here
        fluidRow(
          box(
              title = "General Info", width = 4,
              style = "overflow-y: auto; height: 400px"
          ),

          box(
              title = "Logo", width = 8,
              img(src="BreedingInsight_Primary_RGBColor_400px.png"),
              style = "overflow-y: auto; height: 400px")
          )
        ),
      tabItem(
        tabName = "filtering",
        fluidRow(
          column(width = 3,
              box(width = 12,
                title = "Quality Filtering", status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                fileInput("updog_rdata","Choose Updog Rdata File", accept = ".rda"),
                textInput("output_name", "Output File Name"),
                #sliderInput("hist_bins","Histogram Bins", min = 1, max = 1200, value = c(500), step = 1),
                sliderInput("size_depth","Minimum Read Depth", min = 0, max = 300, value = 10, step = 1),
                #numericInput("size_depth","Minimum Read Depth", min = 0, max = 300, value = 10, step = 1),
                sliderInput("Bias","Bias (Updog filter)", min = 0, max = 10, value = c(0.5,2), step = 0.1),
                #numericInput("Bias_min","Bias minimum (Updog filter)", min = 0, max = 10, value = 0.5, step = 0.1),
                #numericInput("Bias_max","Bias maximum (Updog filter)", min = 0, max = 10, value = 2, step = 0.1),
                numericInput("OD_filter","OD (Updog filter)", min = 0, value = 0.5),
                numericInput("Prop_mis","Prop_mis (Updog filter)", min = 0, max=1, value = 0.05, step = 0.05),
                numericInput("maxpostprob_filter","maxpostprob (Updog filter)", min = 0, value = 0.9, step = 0.1),
                #numericInput("missing_filter","Remove SNPs with >= % missing data", min = 0, max = 1, value = 0.5, step = 0.1),
                #numericInput("missing_filter","Remove Samples with >= % missing data", min = 0, max = 1, value = 0.5, step = 0.1),
                actionButton("start_updog_filter", "Download Filtered Dosage File", icon = icon("download")),
                  div(style="display:inline-block; float:right",dropdownButton(
                    tags$h3("Updog Filter Parameters"),
                    #selectInput(inputId = 'xcol', label = 'X Variable', choices = names(iris)),
                    #selectInput(inputId = 'ycol', label = 'Y Variable', choices = names(iris), selected = names(iris)[[2]]),
                    #sliderInput(inputId = 'clusters', label = 'Cluster count', value = 3, min = 1, max = 9),
                    "Add description of each filter",
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
              tabPanel("Bias Histogram", icon = icon("image"), plotOutput("bias_hist", height = '550px')),
              tabPanel("OD Histogram", icon = icon("image"), plotOutput("od_hist", height = '550px')),
              tabPanel("MaxPostProb Histogram", icon = icon("image"), plotOutput("maxpostprob_hist", height = '550px')),
              tabPanel("ReadDepth Histogram", icon = icon("image"), plotOutput("depth_hist", height = '550px'))
              #tabPanel("SNP Distribution Plot", icon = icon("image"), plotOutput("snp_dist", height = '550px')),
              #tabPanel("SNP % Missing Histogram", icon = icon("image"), plotOutput("missing_snp_hist", height = '550px')),
              #tabPanel("Sample % Missing Histogram", icon = icon("image"), plotOutput("missing_sample_hist", height = '550px')),
              #tabPanel("Summary Statistics", icon = icon("sliders"), tableOutput("dosages"))
              #plotOutput("coverage"), # Placeholder for plot outputs
              )
            ),
          column(width = 3,
            valueBoxOutput("snps"),
            valueBox(1111,"SNPs Retained", icon = icon("dna"), width = NULL, color = "info"),
            valueBox("78%","SNPs Removed", icon = icon("filter-circle-xmark"), width = NULL, color = "info"), #https://rstudio.github.io/shinydashboard/structure.html#tabbox
            box(title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE,
                sliderInput("hist_bins","Histogram Bins", min = 1, max = 1200, value = c(50), step = 1), width = NULL,
                div(style="display:inline-block; float:left",dropdownButton(
                      tags$h3("Save Image"),
                      selectInput(inputId = 'filter_hist', label = 'Figure', choices = c("Bias Histogram", 
                                                                                         "OD Histogram", 
                                                                                         "MaxPostProb Histogram")),
                      selectInput(inputId = 'image_type', label = 'File Type', choices = c("jpeg","pdf","tiff","png"), selected = "jpeg"),
                      sliderInput(inputId = 'image_res', label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                      sliderInput(inputId = 'image_width', label = 'Width', value = 3, min = 1, max = 10, step=0.5),
                      sliderInput(inputId = 'image_height', label = 'Height', value = 3, min = 1, max = 10, step = 0.5),
                      downloadButton("download_filter_hist", "Save"),
                      circle = FALSE,
                      status = "danger", 
                      icon = icon("floppy-disk"), width = "300px",
                      tooltip = tooltipOptions(title = "Click to see inputs!")
                      ))
            )          
          )
        )
      ),
      tabItem(
        tabName = "updog",
        fluidPage(
          fluidRow(
              box(
              title = "Inputs", status = "info", solidHeader = TRUE, collapsible = FALSE, collapsed = FALSE,
              fileInput("madc_file", "Choose MADC File", accept = c(".csv")),
              #checkboxInput("off-targets","Include off-target loci?"),
              #fileInput("sample_file", "Optional: Choose Sample List (disabled)", accept = c(".csv")),
              textInput("output_name", "Output File Name"),
              numericInput("ploidy", "Species Ploidy", min = 1, value = 2),
              selectInput("updog_model", "Updog Model", choices = c("norm","hw","bb","s1","s1pp","f1","f1pp","flex","uniform"), selected = "norm"),
              numericInput("cores", "Number of CPU Cores", min = 1, max = (future::availableCores() - 1), value = 1),
              actionButton("run_analysis", "Run Analysis"),
              div(style="display:inline-block; float:right",dropdownButton(

                    tags$h3("Updog Population Models"),
                    "Model: What form should the prior (genotype distribution) take?\n
                    The following information is from the Updog manual:\n
                    Possible values of the genotype distribution (values of model) are: \n
                    `norm` A distribution whose genotype frequencies are proportional to the density value of a normal
                    with some mean and some standard deviation. Unlike the `bb` and `hw` options, this will
                    allow for distributions both more and less dispersed than a binomial. This seems to be the
                    most robust to violations in modeling assumptions, and so is the default. This prior class was
                    developed in Gerard and Ferrão (2020).
                    `hw` A binomial distribution that results from assuming that the population is in Hardy-Weinberg
                    equilibrium (HWE). This actually does pretty well even when there are minor to moderate
                    deviations from HWE. Though it does not perform as well as the ‘norm‘ option when there
                    are severe deviations from HWE.
                    `bb` A beta-binomial distribution. This is an overdispersed version of `hw` and can be derived
                    from a special case of the Balding-Nichols model.
                    `s1` This prior assumes the individuals are all full-siblings resulting from one generation of selfing. I.e. there is only one parent. This model assumes a particular type of meiotic behavior:
                    polysomic inheritance with bivalent, non-preferential pairing.
                    `f1` This prior assumes the individuals are all full-siblings resulting from one generation of a
                    bi-parental cross. This model assumes a particular type of meiotic behavior: polysomic inheritance with bivalent, non-preferential pairing.
                    `f1pp` This prior allows for double reduction and preferential pairing in an F1 population of tretraploids.
                    `s1pp` This prior allows for double reduction and preferential pairing in an S1 population of tretraploids.
                    `flex` Generically any categorical distribution. Theoretically, this works well if you have a lot of
                    individuals. In practice, it seems to be much less robust to violations in modeling assumptions.
                    `uniform` A discrete uniform distribution. This should never be used in practice.",
                    circle = FALSE,
                    status = "warning", 
                    icon = icon("info"), width = "300px",
                    tooltip = tooltipOptions(title = "Click to see info!")
                ))
              ),
              valueBoxOutput("MADCsnps")
              #valueBox("Help","Updog Manual", icon = icon("globe"), color = "warning")
          ),
         
         fluidRow(
              box(title = "Status", width = 3, collapsible = TRUE, status = "info",
                progressBar(id = "pb_madc", value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
              )
         
         ) 
        
        )  
      
      ),
      tabItem(
        tabName = "pca",
        # Add PCA content here
        fluidRow(
          column(width = 3,
            box(
              title = "Inputs", width = 12, solidHeader = TRUE, status = "info",
              fileInput("dosage_file", "Choose Genotypes File (Samples as columns)", accept = c(".csv")),
              fileInput("passport_file", "Choose Passport File (Sample IDs in first column)", accept = c(".csv")),
              #textInput("output_name", "Output File Name (disabled)"),
              #Dropdown will update after pasport upload
              numericInput("ploidy", "Species Ploidy", min = 1, value = 2),
              actionButton("pca_start", "Run Analysis"),
              #div(style="display:inline-block; float:right",dropdownButton(
              #     tags$h3("PCA info"),
              #      "Model",
              #      circle = FALSE,
              #      status = "warning", 
              #      icon = icon("info"), width = "300px",
              #      tooltip = tooltipOptions(title = "Click to see info!")
              #  )),
              style = "overflow-y: auto; height: 420px"
            ),

            box(
                title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE,collapsed = FALSE, width = 12,
                "Change the PCA plot parameters", br(),
                selectInput('group_info', label = 'Color Variable (eg, Taxon)', choices = NULL),
                materialSwitch('use_cat', label = "Color Specific Category Within Variable?", status = "success"),
                conditionalPanel(condition = "input.use_cat",
                  virtualSelectInput(
                  inputId = "cat_color",
                  label = "Select Category To Color:",
                  choices = NULL,
                  showValueAsTags = TRUE,
                  search = TRUE,
                  multiple = TRUE
                ),
                  selectInput("grey_choice", "Select Grey", choices = c("Light Grey", "Grey", "Dark Grey", "Black"), selected = "Grey")
                ),
                selectInput("color_choice", "Color Palette", choices = c("YlOrRd","YlOrBr","YlGnBu","YlGn",
                                                                                  "Reds","RdPu","Purples","PuRd","PuBuGn","PuBu",
                                                                                  "OrRd","Oranges","Greys","Greens","GnBu","BuPu",
                                                                                  "BuGn","Blues","Set3","Set2","Set1","Pastel2",
                                                                                  "Pastel1","Paired","Dark2","Accent","Spectral",
                                                                                  "RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn",
                                                                                  "PiYG","BrBG"), selected = "Paired"),
                selectInput("pc_X", "X-Axis (2D-Plot only)", choices = c("PC1","PC2","PC3","PC4","PC5"), selected = "PC1"),
                selectInput("pc_Y", "Y-Axis (2D-Plot only)", choices = c("PC1","PC2","PC3","PC4","PC5"), selected = "PC2"),
                div(style="display:inline-block; float:right",dropdownButton(
                      tags$h3("Save Image"),
                      selectInput(inputId = 'image_type', label = 'File', choices = c("jpeg","pdf","tiff","png"), selected = "jpeg"),
                      sliderInput(inputId = 'image_res', label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                      sliderInput(inputId = 'image_width', label = 'Width', value = 3, min = 1, max = 10, step=0.5),
                      sliderInput(inputId = 'image_height', label = 'Height', value = 3, min = 1, max = 10, step = 0.5),
                      downloadButton("download_pca", "Save"),
                      circle = FALSE,
                      status = "danger", 
                      icon = icon("floppy-disk"), width = "300px",
                      tooltip = tooltipOptions(title = "Click to see inputs!")
                      ))
              )
          ),
          column(width = 8,
            box(title = "Passport Data", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "info", collapsed = FALSE,
                  DTOutput('passport_table'),
                  style = "overflow-y: auto; height: 420px"
              ),

            box(
                title = "PCA Plots", status = "info", solidHeader = FALSE, width = 12, height = 550,
                tabsetPanel(
                  tabPanel("3D-Plot",plotlyOutput("pca_plot", height = '460px')),
                  tabPanel("2D-Plot", plotOutput("pca_plot_ggplot", height = '460px')),
                  tabPanel("Scree Plot", plotOutput("scree_plot", height = '460px'))) # Placeholder for plot outputs
              )


          ),
          column(width = 1, downloadButton("downloasdas","Info", color = "warning"))


        )
      ),
      tabItem(
        tabName = "gwas",
        # Add GWAS content here
        fluidRow(
          column(width = 3,
            box(title="Inputs", width = 12, collapsible = TRUE, collapsed = FALSE, status = "info", solidHeader = TRUE,
              fileInput("gwas_file", "Choose Genotypes File", accept = ".csv"),
              fileInput("phenotype_file", "Choose Phenotype File", accept = ".csv"),
              #textInput("output_name", "Output File Name"),
              numericInput("gwas_ploidy", "Species Ploidy", min = 1, value = 2),
              selectInput('trait_info', label = 'Select Trait (eg, Color):', choices = NULL),
              virtualSelectInput(
                inputId = "fixed_info",
                label = "Select Fixed Effects (eg, Group):",
                choices = NULL,
                showValueAsTags = TRUE,
                search = TRUE,
                multiple = TRUE
              ),
              sliderInput("cores", "Number of CPU Cores", min = 1, max = (future::availableCores() - 1), value = 1, step = 1),
              actionButton("gwas_start", "Run Analysis"),
              #downloadButton("download_pca", "Download All Files"),
              #plotOutput("pca_plot"), # Placeholder for plot outputs
              #checkboxGroupInput("files_to_download", "Select files to download:",
                             #choices = c("PC1vPC2 plot", "PC2vPC3 plot"), selected = c("table1", "table2"))
              div(style="display:inline-block; float:right",dropdownButton(
                    tags$h3("GWAS Parameters"),
                    #selectInput(inputId = 'xcol', label = 'X Variable', choices = names(iris)),
                    #selectInput(inputId = 'ycol', label = 'Y Variable', choices = names(iris), selected = names(iris)[[2]]),
                    #sliderInput(inputId = 'clusters', label = 'Cluster count', value = 3, min = 1, max = 9),
                    "Add description of each filter",
                    circle = FALSE,
                    status = "warning", 
                    icon = icon("info"), width = "300px",
                    tooltip = tooltipOptions(title = "Click to see info!")
                ))#,
              #style = "overflow-y: auto; height: 550px"

              ),
            box(title = "Status", width = 12, collapsible = TRUE, status = "info",

                progressBar(id = "pb_gwas", value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
              )
            ),

          column(width = 6,
            box(
              title = "Plots", status = "info", solidHeader = FALSE, width = 12, height = 600,
              tabsetPanel(
                tabPanel("BIC Plot", plotOutput("bic_plot", height = "500px")),
                tabPanel("Manhattan Plot", plotOutput("manhattan_plot", height = "500px")),
                tabPanel("QQ Plot", plotOutput("qq_plot", height = "500px")),
                tabPanel("BIC Table", DTOutput("bic_table"),style = "overflow-y: auto; height: 500px"),
                tabPanel("Outlier SNPs", DTOutput('gwas_stats'),style = "overflow-y: auto; height: 500px")

                )
              
              )

            ),
          
          column(width = 3,
            valueBox(0,"Outlier SNPs", icon = icon("dna"), width = NULL, color = "info"),
            valueBox("0","QTLs Detected", icon = icon("dna"), width = NULL, color = "info"), #https://rstudio.github.io/shinydashboard/structure.html#tabbox
            box(title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
                #sliderInput("hist_bins","Histogram Bins", min = 1, max = 1200, value = c(50), step = 1), width = NULL,
                selectInput("model_select", label = "Model Selection", choices = c("all","1-dom","2-dom","additive","general","diplo-general","diplo-additive")),
                div(style="display:inline-block; float:left",dropdownButton(
                      tags$h3("Save Image"),
                      selectInput(inputId = 'gwas_figures', label = 'Figure', choices = c("BIC Plot", 
                                                                                         "Manhattan Plot", 
                                                                                         "QTL Statistics",
                                                                                         "QQ Plot")),
                      selectInput(inputId = 'image_type', label = 'File Type', choices = c("jpeg","pdf","tiff","png"), selected = "jpeg"),
                      sliderInput(inputId = 'image_res', label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                      sliderInput(inputId = 'image_width', label = 'Width', value = 3, min = 1, max = 10, step=0.5),
                      sliderInput(inputId = 'image_height', label = 'Height', value = 3, min = 1, max = 10, step = 0.5),
                      downloadButton("download_gwas_table", "Save"),
                      circle = FALSE,
                      status = "danger", 
                      icon = icon("floppy-disk"), width = "300px",
                      tooltip = tooltipOptions(title = "Click to see inputs!")
                      ))
            )          
          
          )
          
        )
      
      ),
      tabItem(
        tabName = "diversity",
        # Add GWAS content here
        fluidRow(
          column(width = 3,
            box(title="Inputs", width = 12, collapsible = TRUE, collapsed = FALSE, status = "info", solidHeader = TRUE,
              fileInput("diversity_file", "Choose Genotypes File"),
              #fileInput("phenotype_file", "Choose Phenotype File"),
              textInput("output_name", "Output File Name"),
              numericInput("ploidy", "Species Ploidy", min = 1, value = 2),
              #numericInput("cores", "Number of CPU Cores", min = 1, max = (future::availableCores() - 1), value = 1),
              actionButton("diversity_start", "Run Analysis"),
              #downloadButton("download_pca", "Download All Files"),
              #plotOutput("pca_plot"), # Placeholder for plot outputs
              #checkboxGroupInput("files_to_download", "Select files to download:",
                             #choices = c("PC1vPC2 plot", "PC2vPC3 plot"), selected = c("table1", "table2"))
              div(style="display:inline-block; float:right",dropdownButton(
                    tags$h3("Diversity Parameters"),
                    #selectInput(inputId = 'xcol', label = 'X Variable', choices = names(iris)),
                    #selectInput(inputId = 'ycol', label = 'Y Variable', choices = names(iris), selected = names(iris)[[2]]),
                    #sliderInput(inputId = 'clusters', label = 'Cluster count', value = 3, min = 1, max = 9),
                    "Add description of each filter",
                    circle = FALSE,
                    status = "warning", 
                    icon = icon("info"), width = "300px",
                    tooltip = tooltipOptions(title = "Click to see info!")
                ))#,
              #style = "overflow-y: auto; height: 550px"

              )
            ),

          column(width = 6,
            box(
              title = "Plots", status = "info", solidHeader = FALSE, width = 12, height = 550,
              tabsetPanel(
                tabPanel("Dosage Plot"),
                tabPanel("AF Plot"),
                tabPanel("MAF Plot"),
                tabPanel("PIC Plot"),
                tabPanel("Diversity Statistics")

                )
              
              )

            ),
          column(width = 3,
            #valueBoxOutput("snps"),
            valueBox(0.365,"Mean Heterozygosity", icon = icon("dna"), width = NULL, color = "info"),
            valueBox(0.245,"Mean Minor-Allele-Frequency", icon = icon("dna"), width = NULL, color = "info"), #https://rstudio.github.io/shinydashboard/structure.html#tabbox
            valueBox(0.301,"Mean PIC", icon = icon("dna"), width = NULL, color = "info"),
            box(title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE,
                sliderInput("hist_bins","Histogram Bins", min = 1, max = 1200, value = c(50), step = 1), width = NULL,
                div(style="display:inline-block; float:left",dropdownButton(
                      tags$h3("Save Image"),
                      selectInput(inputId = 'filter_hist', label = 'Figure', choices = c("Bias Histogram", 
                                                                                         "OD Histogram", 
                                                                                         "MaxPostProb Histogram")),
                      selectInput(inputId = 'image_type', label = 'File Type', choices = c("jpeg","pdf","tiff","png"), selected = "jpeg"),
                      sliderInput(inputId = 'image_res', label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                      sliderInput(inputId = 'image_width', label = 'Width', value = 3, min = 1, max = 10, step=0.5),
                      sliderInput(inputId = 'image_height', label = 'Height', value = 3, min = 1, max = 10, step = 0.5),
                      downloadButton("download_diversity", "Save"),
                      circle = FALSE,
                      status = "danger", 
                      icon = icon("floppy-disk"), width = "300px",
                      tooltip = tooltipOptions(title = "Click to see inputs!")
                      ))
            )          
          )
          
          )
      
      ),
      tabItem(
        tabName = "qtl",
        fluidPage(
          # Add QTL content here
          fileInput("madc_file", "Choose MADC File"),
          textInput("output_name", "Output File Name"),
          numericInput("ploidy", "Ploidy", min = 1, value = 2),
          numericInput("cores", "Number of CPU Cores", min = 1, max = (future::availableCores() - 1), value = 1),
          actionButton("run_analysis", "Run Analysis"),
          downloadButton("download_qtl", "Download All Files"),
          checkboxGroupInput("files_to_download", "Select files to download:",
                             choices = c("table1", "table2"), selected = c("table1", "table2"))
        )
      ),
      tabItem(
        tabName = "prediction",
        fluidPage(
          # Add genomic prediction content here
          fileInput("madc_file", "Choose MADC File"),
          textInput("output_name", "Output File Name"),
          numericInput("ploidy", "Ploidy", min = 1, value = 2),
          numericInput("cores", "Number of CPU Cores", min = 1, max = (future::availableCores() - 1), value = 1),
          actionButton("run_analysis", "Run Analysis"),
          downloadButton("download_prediction", "Download All Files"),
          checkboxGroupInput("files_to_download", "Select files to download:",
                             choices = c("table1", "table2"), selected = c("table1", "table2"))
        )
      ),
      tabItem(
        tabName = "code",
        fluidPage(
          # Add source code here
        )
      ),
      tabItem(
        tabName = "help",
        fluidPage(
          # Add Help content here
        )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Define reactive values to store generated figures and tables
  figures <- reactiveValues()
  tables <- reactiveValues()

  output$resultText <- renderText({
    # Example text output, adapt based on your script's processing
    paste("Processing file:", input$madc_file$name, "\n",
          "with output name:", input$output_name, "\n",
          "and ploidy:", input$ploidy, "\n",
          "on", ifelse(input$cores == 0, "auto-detected number of cores", paste(input$cores, "cores")))
  })

  #Add server configurations
  options(shiny.maxRequestSize = 100000 * 1024^2)  # Set maximum upload size to 100GB

  #shiny.maxRequestSize = 10000 * 1024^2; # 10 GB <- This is for a future limit when using BI's server remotely
  
  
  #PCA reactive values
  pca_data <- reactiveValues(
    pc_df_pop = NULL,
    variance_explained = NULL,
    my_palette = NULL

  )
  #SNP reactive values
  snp_counts <- reactiveValues(
    snp_count = 0,
    snp_percent = 0

  )

  snp_number <- reactiveVal(0)

  #SNP counts value box
  output$MADCsnps <- renderValueBox({
    valueBox(snp_number(), "Markers in MADC File", icon = icon("dna"), color = "info")
    })

##This is for performing Updog Dosage Calling
  observeEvent(input$run_analysis, {
    # Get inputs
    madc_file <- input$madc_file$datapath
    sample_file <- input$sample_file$datapath
    output_name <- input$output_name
    ploidy <- input$ploidy
    cores <- input$cores
    model_select <- input$updog_model
    

    #Status
    updateProgressBar(session = session, id = "pb_madc", value = 0, title = "Formatting Input Files")

    # Perform analysis
    # Modify the get_counts function to accept the MADC file path as an argument
    get_counts <- function(madc_file, output_name) {
      # This function takes the MADC file as input and generates a Ref and Alt counts dataframe as output
      # Note: This assumes that the first 7 rows are not useful here like in the Strawberry DSt23-8501_MADC file
      
      # Read the madc file
      madc_df <- read.csv(madc_file, sep = ',', skip = 7)
      
      # Retain only the Ref and Alt haplotypes
      filtered_df <- madc_df[grep("\\|Ref$|\\|Alt$", madc_df$AlleleID), ]
      
      # Save the csv file for review and use in R
      df_name <- paste0(output_name,'_MADC_alt_ref_counts.csv')
      
      #write.csv(filtered_df, file = df_name, row.names = FALSE)
      return(filtered_df)
    }

    
    # Call the get_counts function with the specified MADC file path and output file path
    result_df <- get_counts(madc_file, output_name)

    #Number of SNPs
    snp_number <- (nrow(result_df) / 2)
    
    output$table1 <- renderTable({
    # Generate table
      result_df
    })
    
    #Get the alt, ref, and size matrix for use in Updog
    #Add functionality here to stop the script if indentical() is False
    get_matrices <- function(result_df) {
      #This function takes the dataframe of ref and alt counts for each sample, and converts them to ref, alt, and size(total count) matrices for Updog
      
      update_df <- result_df
      
      # Filter rows where 'AlleleID' ends with 'Ref'
      ref_df <- subset(update_df, grepl("Ref$", AlleleID))
      
      # Filter rows where 'AlleleID' ends with 'Alt'
      alt_df <- subset(update_df, grepl("Alt$", AlleleID))
      
      #Ensure that each has the same SNPs and that they are in the same order
      identical(alt_df$CloneID,ref_df$CloneID)
      
      ###Convert the ref and alt counts into matrices with the CloneID as the index
      #Set SNP names as index
      row.names(ref_df) <- ref_df$CloneID
      row.names(alt_df) <- alt_df$CloneID
      
      #Remove unwanted columns and convert to matrix
      ref_matrix <- as.matrix(ref_df[, -c(1:16)])
      alt_matrix <- as.matrix(alt_df[, -c(1:16)])
      
      #Make the size matrix by combining the two matrices
      size_matrix <- (ref_matrix + alt_matrix)
      
      #Count the number of cells with 0 count to estimate missing data
      # Count the number of cells with the value 0
      count_zeros <- sum(size_matrix == 0)
      
      # Print the result
      ratio_missing_data <- count_zeros / length(size_matrix)
      cat("Ratio of missing data =", ratio_missing_data, "\n")
      
      # Return the ref and alt matrices as a list
      matrices_list <- list(ref_matrix = ref_matrix, size_matrix = size_matrix)
      return(matrices_list)
    }
    

    #Status
    updateProgressBar(session = session, id = "pb_madc", value = 40, title = "Dosage Calling in Progress")

    #Call the get_matrices function
    matrices <- get_matrices(result_df)
    
    #Run Updog 
    #I initially used the "norm" model
    #I am also taking the ploidy from the max value in the
    print('Performing Updog dosage calling') 
    mout <- multidog(refmat = matrices$ref_matrix, 
                     sizemat = matrices$size_matrix, 
                     ploidy = as.numeric(ploidy),  
                     model = model_select,
                     nc = cores)
    
    #Get genotype matrix of dosage calls
    genomat <- format_multidog(mout, varname = "geno")
    #Save the matrix as a csv file
    updog_file <- paste0(output_name,'_MADC_alt_ref_counts_unfiltered_dose_from_updog_norm_genotype_matrix.csv')
    write.csv(genomat,file=updog_file)

    #Filter dosage calls (I think this is the updog recommended)
    #mout_cleaned <- filter_snp(mout, prop_mis < 0.2 & bias > 0.5 & bias < 2 & od > 0.05) #Recommended filtering by updog
    
    #Save the filtered dosage matrix
    #genomat_cleaned <- format_multidog(mout_cleaned, varname = "geno")
    #head(genomat_cleaned)
    
    #cleaned_name <- paste0(output_name,'_MADC_alt_ref_counts_filtered_prop_miss_0.2_bias_0.5-2_updog_norm_model_dosage_genotype_matrix.csv')
    #Save the matrix as a csv file
    #write.csv(genomat_cleaned,file= cleaned_name)

    #Save rda file for filtering
    save(mout, result_df, file = paste0(output_name,"_MADC_unfiltered_dose_from_updog.rda"))
    
    #Reactive item
    #output$table2 <- renderTable({
    # Generate table
    #  genomat_cleaned
    #}) 
    
    # Display analysis result
    #output$analysis_result <- renderText({
    #  "result" #Can add a variable to print text or figures
    #})
  
  # Update reactive values with generated figures and tables
  #figures$plot1 <- heatmap(G.mat.updog, labCol = NA)# Your plot object
  #tables$table1 <- result_df# Your table object
  #tables$table2 <- genomat_cleaned
  
  #Status
  updateProgressBar(session = session, id = "pb_madc", value = 100, title = "Finished")
  
  })
  
  #Updog filtering
  observeEvent(input$start_updog_filter, {
    req(input$Prop_mis, input$Bias, input$OD_filter, input$size_depth)

    #Variables
    Prop_mis <- input$Prop_mis
    Bias_min <- input$Bias[1]
    Bias_max <- input$Bias[2]
    OD_filter <- input$OD_filter
    size_depth <- input$size_depth
    output_name <- input$output_name


    #Input file
    madc_path <- input$updog_rdata$datapath
    load(madc_path)

    #Number of SNPs
    length(mout$snpdf$bias)

    #Filter dosage calls (I think this is the updog recommended)
    #mout_cleaned <- filter_snp(mout, prop_mis < input$Prop_mis & bias > input$Bias_min & bias < input$Bias_max & od > input$OD_filter) #Recommended filtering by updog
    print("Filtering")
    #Needing to paste the updog filter_snp() method directly to work with shiny
    filter_snp_custom <- function(x, expr) {
      assertthat::assert_that(is.multidog(x))
      cond <- eval(expr = substitute(expr), envir = x$snpdf)
      x$snpdf <- x$snpdf[cond, , drop = FALSE]
      goodsnps <- x$snpdf$snp
      x$inddf <- x$inddf[x$inddf$snp %in% goodsnps, , drop = FALSE]
      return(x)
    }

    mout_cleaned <- filter_snp_custom(mout, prop_mis < Prop_mis & bias > Bias_min & bias < Bias_max & od < OD_filter) #Recommended filtering by updog
    print("Filtering complete")
    # Filter using the dynamically constructed expression
    #mout_cleaned <- filter_snp(mout, filter_expr)
    
    #Notify user if all SNPs were removed and stop filtering
    if (length(mout_cleaned$snpdf$snp) == 0) {
      # If condition is met, show notification toast
        shinyalert(
            title = "Oops",
            text = "No SNPs were found after cleaning the data.\n Adjust filtering parameters.",
            size = "xs",
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            type = "info",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE
          )
                         
      
      # Stop the observeEvent gracefully
      return()
    }
    

    # Replace values in "geno" column with NA where "size" is less than input value
    mout_cleaned$inddf$geno[mout_cleaned$inddf$size < size_depth] <- NA

    #Save the filtered dosage matrix
    genomat_cleaned <- format_multidog(mout_cleaned, varname = "geno")
    print("Generated genotype matrix")
    cleaned_name <- paste0(output_name,'_MADC_alt_ref_counts_filtered_updog_dosage_genotype_matrix.csv')
    #Save the matrix as a csv file
    write.csv(genomat_cleaned,file= cleaned_name)
    print("Write csv complete")
  
  })

  observeEvent(input$updog_rdata, {
    #req(filter_hist$hist_bin_value)

    #Variables


    #Input file
    madc_path <- input$updog_rdata$datapath
    load(madc_path)

    #Number of SNPs
    length(mout$snpdf$bias)

    ###Bias

    #Histogram
    output$bias_hist <- renderPlot({
      hist(mout$snpdf$bias, 
          main = "Unfiltered SNP bias histogram",
          xlab = "bias",
          ylab = "SNPs",
          col = "lightblue",
          border = "black",
          xlim = c(0,5),
          breaks = as.numeric(input$hist_bins))
      axis(1, at = seq(0, 5, by = .2), labels = rep("", length(seq(0, 5, by = 0.2))))  # Add ticks
      abline(v = mean(mout$snpdf$bias), col = "red", lty = 2)  # Mean line
      abline(v = median(mout$snpdf$bias), col = "green", lty = 2)  # Median line
      abline(v = 0.5, col = "black", lty = 2)  # proposed lower line
      abline(v = 2, col = "black", lty = 2)  # proposed upper line
    })

    ###OD
    quantile(mout$snpdf$od, 0.95)
    #Histogram
    output$od_hist <- renderPlot({

      hist(mout$snpdf$od, 
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
      abline(v = mean(mout$snpdf$od), col = "red", lty = 2)  # Mean line
      abline(v = median(mout$snpdf$od), col = "green", lty = 2)  # Median line
      abline(v = 0.05, col = "black", lty = 2)  # proposed filter by updog

    })

    ##MAXPOSTPROB

    #Histogram

    output$maxpostprob_hist <- renderPlot({

        #Histogram
        hist(mout$inddf$maxpostprob, 
            main = "SNP max post probabilities histogram",
            xlab = "Max Post Probability",
            ylab = "Genomic Sites",
            col = "lightblue",
            border = "black",
            xlim = c(0,1),
            breaks = as.numeric(input$hist_bins))
        axis(1, at = seq(0, 1, by = .1), labels = rep("", length(seq(0, 1, by = 0.1))))  # Add ticks

        # Add vertical lines
        abline(v = mean(mout$inddf$maxpostprob), col = "red", lty = 2)  # Mean line
        abline(v = median(mout$inddf$maxpostprob), col = "green", lty = 2)  # Median line
        abline(v = quantile(mout$inddf$maxpostprob, 0.95), col = "blue", lty = 2) 

    })

    ##Read Depth (I would prefer that this show the mean depth for SNPs or Samples instead of all loci/sample cells)

    quantile(mout$inddf$size, 0.95)
    #Histogram
    output$depth_hist <- renderPlot({

      hist(mout$inddf$size, 
          main = "Unfiltered SNP overdispersion parameter histogram",
          xlab = "Read Depth per SNP/Sample",
          ylab = "Genomic Sites",
          col = "lightblue",
          border = "black",
          xlim = c(0,1000),
          breaks = as.numeric(input$hist_bins))
      axis(1, at = seq(0, 1000, by = 20), labels = rep("", length(seq(0, 1000, by = 20))))  # Add ticks
      abline(v = 0.05, col = "black", lty = 2)  # proposed filter by updog

      # Add vertical lines
      abline(v = mean(mout$inddf$size), col = "red", lty = 2)  # Mean line
      abline(v = median(mout$inddf$size), col = "green", lty = 2)  # Median line
      #abline(v = 0.05, col = "black", lty = 2)  # proposed filter by updog

    })

    #Filter dosage calls (I think this is the updog recommended)

    #mout_cleaned <- filter_snp(mout, prop_mis < as.numeric(input$Prop_mis) & bias > input$Bias_min & bias < input$Bias_max & od > input$OD_filter) #Recommended filtering by updog

    # Replace values in "geno" column with NA where "size" is less than input value
    #mout_cleaned$inddf$geno[mout_cleaned$inddf$size < as.numeric(input$size_depth)] <- NA

    #Save the filtered dosage matrix
    #genomat_cleaned <- format_multidog(mout_cleaned, varname = "geno")

  })
  
  #PCA dropdown

  data <- reactiveValues(info_df = NULL)

  # Update dropdown menu choices based on uploaded passport file
  observeEvent(input$passport_file, {
    info_df <- read.csv(input$passport_file$datapath, header = TRUE, check.names = FALSE)
    data$info_df <- info_df  # Store info_df in reactiveValues

    updateSelectInput(session, "group_info", choices = colnames(info_df))

    output$passport_table <- renderDT({info_df}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 4)
    )
  })

  #PCA specific category selection
  observeEvent(input$group_info, {
    req(data$info_df)

    #updateMaterialSwitch(session, inputId = "use_cat", status = "success")

    # Get selected column name
    selected_col <- input$group_info
  
    # Extract unique values from the selected column
    unique_values <- unique(data$info_df[[selected_col]])

    #Add category selection
    updateVirtualSelect("cat_color", choices = unique_values, session = session)

    })

  #PCA events
  observeEvent(input$pca_start, {
    # Get inputs
    geno <- input$dosage_file$datapath
    pedigree_df <- input$passport_file$datapath
    g_info <- as.character(input$group_info)
    output_name <- input$output_name
    ploidy <- input$ploidy
    
    PCX <- input$pc_X
    PCY <- input$pc_Y
    ### Perform analysis

    genomat <- read.csv(geno, header = TRUE, row.names = 1, check.names = FALSE)

    #Passport info
    # Sample dataframe with a column of taxon names
    info_df <- read.csv(pedigree_df, header = TRUE, check.names = FALSE)

    # Print the modified dataframe
    row.names(info_df) <- info_df[,1]

    #Plotting
    #First build a relationship matrix using the genotype values
    #First build a relationship matrix using the genotype values
    G.mat.updog <- AGHmatrix::Gmatrix(t(genomat), method = "VanRaden", ploidy = as.numeric(ploidy), missingValue = "NA")
    
    #PCA...maybe add in a DAPC?
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

    #Add the male parent information for each sample
    pc_df_pop <- merge(pc_df, info_df, by.x = "row.names", by.y = "row.names", all.x = TRUE)
    pc_df_pop[[g_info]] <- as.factor(pc_df_pop[[g_info]])

    #Update global variable
    pca_dataframes <- pc_df_pop

    #Output PC file for manoj
    #write.csv(pc_df_pop, file = 'DAl22-7535_lightly_filtered_bias0.5-2_updog_2382_SNPs_MAF_0.05_PCA_metadata.csv')

    # Generate a distinct color palette
    unique_countries <- unique(pc_df_pop[[g_info]])
    #my_palette <- randomcoloR::distinctColorPalette(length(unique_countries))
    palette <- brewer.pal(length(unique_countries),input$color_choice)
    my_palette <- colorRampPalette(palette)(length(unique_countries))


    # Store processed data in reactive values
    pca_data$pc_df_pop <- pc_df_pop
    pca_data$variance_explained <- variance_explained
    pca_data$my_palette <- my_palette

    #End of PCA section
    }
  )

  #2D PCA plotting
  # Render PCA plot reactive to changes in control options
 # output$pca_plot_ggplot <- renderPlot({
 #   req(pca_data$pc_df_pop, pca_data$variance_explained, pca_data$my_palette)
 #   
 #   #Generate colors
 #   unique_countries <- unique(pca_data$pc_df_pop[[input$group_info]])
 #   #my_palette <- randomcoloR::distinctColorPalette(length(unique_countries))
 #   palette <- brewer.pal(length(unique_countries),input$color_choice)
 #   my_palette <- colorRampPalette(palette)(length(unique_countries))

  #  if input$use_cat == TRUE {}
  #
  #  ggplot(pca_data$pc_df_pop, aes(x = pca_data$pc_df_pop[[input$pc_X]], y = pca_data$pc_df_pop[[input$pc_Y]], color = pca_data$pc_df_pop[[input$group_info]])) +
  #    geom_point(size = 2, alpha=0.8) +
  #    scale_color_manual(values = my_palette) +
  #    guides(color = guide_legend(override.aes = list(size = 5.5), nrow = 17)) +
  #    theme_minimal() +
  #    theme(
  #      panel.border = element_rect(color = "black", fill = NA),
  #      legend.text = element_text(size = 14),
  #      axis.title = element_text(size = 14),
  #      axis.text = element_text(size = 12),
  #      legend.title = element_text(size = 16)
  #    ) +
  #    labs(
  #      x = paste0(input$pc_X, "(", pca_data$variance_explained[as.numeric(substr(input$pc_X, 3, 3))], "%)"),
  #      y = paste0(input$pc_Y, "(", pca_data$variance_explained[as.numeric(substr(input$pc_Y, 3, 3))], "%)"),
  #      color = input$group_info)
  #})

  output$pca_plot_ggplot <- renderPlot({
    req(pca_data$pc_df_pop, pca_data$variance_explained, pca_data$my_palette, input$grey_choice)
    
    # Generate colors
    unique_countries <- unique(pca_data$pc_df_pop[[input$group_info]])
    palette <- brewer.pal(length(unique_countries), input$color_choice)
    my_palette <- colorRampPalette(palette)(length(unique_countries))

    # Define a named vector to map input labels to grey values
    label_to_value <- c("Light Grey" = "grey80",
                        "Grey" = "grey",
                        "Dark Grey" = "grey40",
                        "Black" = "black")

    # Get the corresponding value based on the selected grey
    selected_grey <- label_to_value[[input$grey_choice]]

    if (input$use_cat) {
      cat_colors <- c(input$cat_color, "grey")
      ggplot(pca_data$pc_df_pop, aes(x = pca_data$pc_df_pop[[input$pc_X]], 
                                      y = pca_data$pc_df_pop[[input$pc_Y]], 
                                      color = factor(pca_data$pc_df_pop[[input$group_info]]))) +
        geom_point(size = 2, alpha = 0.8) +
        scale_color_manual(values = setNames(c(my_palette, "grey"), cat_colors), na.value = selected_grey) +
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

    } else {
      ggplot(pca_data$pc_df_pop, aes(x = pca_data$pc_df_pop[[input$pc_X]], 
                                      y = pca_data$pc_df_pop[[input$pc_Y]], 
                                      color = pca_data$pc_df_pop[[input$group_info]])) +
        geom_point(size = 2, alpha = 0.8) +
        scale_color_manual(values = my_palette) +
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
    }
})

  #3D PCA plotting
  output$pca_plot <- renderPlotly({
    #Plotly
    req(pca_data$pc_df_pop, pca_data$variance_explained, pca_data$my_palette)

    #Generate colors
    unique_countries <- unique(pca_data$pc_df_pop[[input$group_info]])
    #my_palette <- randomcoloR::distinctColorPalette(length(unique_countries))
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

  output$scree_plot <- renderPlot({
    #PCA scree plot
    req(pca_data$variance_explained)

    var_explained <- pca_data$variance_explained

    # Create a data frame for plotting
    plot_data <- data.frame(PC = 1:10, Variance_Explained = var_explained[1:10])

    # Use ggplot for plotting
    ggplot(plot_data, aes(x = PC, y = Variance_Explained)) + 
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

    })


  gwas_vars <- reactiveValues(
    gwas_df = NULL


    )
  observeEvent(input$phenotype_file, {
    info_df <- read.csv(input$phenotype_file$datapath, header = TRUE, check.names = FALSE, nrow = 0)
    trait_var <- colnames(info_df)
    trait_var <- trait_var[2:length(trait_var)]
    updateSelectInput(session, "trait_info", choices = c("All", trait_var))
    updateVirtualSelect("fixed_info", choices = trait_var, session = session)

    #output$passport_table <- renderDT({info_df}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 4)
    #)
  })

  #GWAS analysis (Shufen Chen and Meng Lin pipelines)
  observeEvent(input$gwas_start, {

    #pheno_file <- read.csv(input$pheno_file$datapath)
    #geno_file <- read.csv(input$genotype_file$datapath)
    cores <- input$cores

    #Status
    updateProgressBar(session = session, id = "pb_gwas", value = 0, title = "Uploading Data")

    #Make subset phenotype file (need to develop alternative method that does not save a new phenotype file each time.)
    #I think I can subset the read.GWAS file pheno and fixed categories (data@pheno[,c("trait")]) and data@fixed = phenotype_file[,c("List of fixed traits")]
    phenotype_file <- read.csv(input$phenotype_file$datapath, header = TRUE, check.names = FALSE)

    ids <- colnames(phenotype_file)[1]
    traits <- input$trait_info
    fixed <- input$fixed_info
    included_var <- c(ids, traits, fixed)
    ploidy <- as.numeric(input$gwas_ploidy)

    phenotype_file <- phenotype_file[,included_var]

    #Save new phenotype file with selected traits and fixed effects
    write.csv(phenotype_file, file = "phenotypes_selected_traits.csv", row.names = FALSE)

    #Remove the phenotype_file from memory
    rm(phenotype_file)

    #Status
    updateProgressBar(session = session, id = "pb_gwas", value = 5, title = "Upload Complete: Now Formatting GWASpoly Data")

    data <- GWASpoly::read.GWASpoly(ploidy= ploidy, pheno.file="phenotypes_selected_traits.csv", geno.file=input$gwas_file$datapath,
                          format="numeric", n.traits=length(traits), delim=",") #only need to change files here


    data.loco <- set.K(data,LOCO=F,n.core= as.numeric(cores))

    ####Pheno, kinship, PCs from results of GWASpoly
    GE<- data@pheno
    names(GE)
    colnames(GE)[1]<-"Genotype"

    ## kinship
    Kin<- data.loco@K$all

    ## PCs
    PC_all<- eigen(data.loco@K$all)$vectors
    rownames(PC_all) <- rownames(data.loco@K$all)
    PCs<-PC_all[,1:10]
    colnames(PCs)<-c(paste("PC",1:10,sep=""))

    ##
    taxa<-intersect(GE$Genotype,intersect(rownames(PCs),rownames(Kin)))
    PCs<-PCs[which(rownames(PCs) %in% taxa),]
    PCs<-PCs[order(rownames(PCs)),]

    GE<-GE[which(GE$Genotype %in% taxa),]
    GE<-GE[order(GE$Genotype),]

    Kin<-Kin[which(rownames(Kin) %in% taxa),which(rownames(Kin) %in% taxa)] # need check the matrix after this step
    Kin<-Kin[order(rownames(Kin)),order(colnames(Kin))]

    which(rownames(Kin)!=rownames(PCs))
    which(rownames(Kin)!=GE$Genotype)

    #### calculate BIC
    #Status
    updateProgressBar(session = session, id = "pb_gwas", value = 20, title = "Formatting Complete: Now Calculating BIC")
    
    source("FUN/MyFun_BIC_Meng.R") #change directory in your case

    PC<-as.matrix(PCs)
    K=as.matrix(Kin)

    posdefmat <- function(mat) {
      if (is.positive.definite(round(mat, 18))) {
        g = mat
      }
      else {
        g <-nearPD(mat)$mat
        warning("The matrix was adjusted for the nearest positive definite matrix")
      }
      return(g)
    }

    kin.adj<-posdefmat(K)
    kin.test<-as.matrix(kin.adj)


    for (i in 2:ncol(GE)){
  
      #model selection
      y=as.numeric(GE[,i])
  
      BICs<-CalcBIC(y=y,PC=PC,K=kin.test)
      BICs$BIC$withK
      BICs$BIC$withoutK
  
      plotBICs<-cbind(rbind.data.frame(BICs$BIC$withK,BICs$BIC$withoutK),rep(c("w/Kinship","no Kinship"),each=nrow(BICs$BIC$withK)))
      colnames(plotBICs)[ncol(plotBICs)]<-"RelationshipMatrix"
      plotBICs$n.PC<-factor(plotBICs$n.PC,levels=c("0","1","2","3","4","5",
                                                  "6","7","8","9","10"))
      #write.table(plotBICs,paste(colnames(GE)[i],"_model_selection_BIC.txt",sep=""),row.names=F,sep="\t",quote=F)
  
      output$bic_table <- renderDT({plotBICs}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5)
        )

      #setwd("~/Desktop/Alfalfa_GWAS/GWAS by year/ModelSel") #change directory in your case
      #pdf(paste(colnames(GE)[i],"_model_selection_BIC.pdf",sep=""),height=4,width=7)
      
      output$bic_plot <- renderPlot({
      p1<-ggplot(plotBICs, aes(x=n.PC, y=BIC,group=RelationshipMatrix)) + 
        geom_point(aes(color=RelationshipMatrix))+
        geom_line(aes(color=RelationshipMatrix))+
        #geom_point(aes(shape=ISTL1_intro,size=2))+
        theme(axis.text.x = element_text(angle =0),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        labs(color="Relationship Matrix",x = "Number of PCs",y="BIC")
      print(p1)

      })
      #dev.off()
      
      #Status
      updateProgressBar(session = session, id = "pb_gwas", value = 40, title = "BIC Complete: Now Performing GWAS")

      ##GWAS based on model selection
      N <- nrow(data@pheno) #Population size
      model <- c("additive","1-dom","2-dom","general","diplo-general","diplo-additive")
  
      BIC_min <- plotBICs[which.min(plotBICs$BIC),]
      if(BIC_min$n.PC == 0){params <- set.params(geno.freq = 1 - 5/N)}else{params <- set.params(geno.freq = 1 - 5/N,n.PC = BIC_min$n.PC)}
      if(BIC_min$RelationshipMatrix=="no Kinship"){
        data.loco.scan <- GWASpoly(data=data,models=model,traits=colnames(data@pheno[i]),params=params,n.core= as.numeric(cores))}else{
          data.loco.scan <- GWASpoly(data=data.loco,models=model,traits=colnames(data@pheno[i]),params=params,n.core= as.numeric(cores))}
                                                                                                                                                                
      data2 <- set.threshold(data.loco.scan,method="M.eff",level=0.05)
      

      #Save manhattan plots to list (only for single trait analysis)
      #if length(traits) == 1
      manhattan_plot_list <- list()

      #Output the manhattan plots
      output$manhattan_plot <- renderPlot({

          print(manhattan_plot_list[[input$model_select]])

        })

      #plot for six models per trait
      #png(file=paste("Manplot_",colnames(data@pheno[i]),".png",sep=""),width=800, height=550) #change directory in your case
      manhattan_plot_list[["all"]] <- manhattan.plot(data2,traits=colnames(data@pheno[i]))+geom_point(size=3)+theme(text = element_text(size = 25),strip.text = element_text(face = "bold"))
      #print(p1)
      #dev.off()     
  
      #get most significant SNPs per QTL file
      qtl <- get.QTL(data=data2,traits=colnames(data@pheno[i]),bp.window=5e6)
      #knitr::kable(qtl)
      qtl_d <- data.frame(qtl)
      
      output$gwas_stats <-  renderDT({qtl_d}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5)
    )

      #write.csv(qtl_d,file = paste("GWAS_by_year_SNP_",colnames(data@pheno[i]),"_SNP",".csv",sep=""))#change directory in your case
      
      #Status
      updateProgressBar(session = session, id = "pb_gwas", value = 80, title = "GWAS Complete: Now Plotting Results")

      #get qqplot
      data_qq <- cbind.data.frame(SNP=data.loco.scan@map$Marker,Chr=data.loco.scan@map$Chrom, Pos=data.loco.scan@map$Position,10^(-data.loco.scan@scores[[colnames(data@pheno[i])]]))
      
      source("FUN/CMplot.r") #Obtained the CMplot code from GitHub and made edits to allow inline plotting for shiny app

      output$qq_plot <- renderPlot({
        CMplot_shiny(data_qq,plot.type="q",col=c(1:8),
                 ylab.pos=2,
                 #threshold= 10^(-4.54042),
                 #signal.cex=1.2,signal.col="red",signal.pch=c(1:8),
                 file.name=colnames(data@pheno[i]),
                 conf.int=FALSE,
                 box=F,multraits=TRUE,file.output=FALSE)

      })
      #plot for each model per trait
      for (j in 1:6) {
        print(j)
    
        data.loco.scan_2 <- GWASpoly(data=data.loco,models=model[j],
                                 traits=colnames(data@pheno[i]),params=params,n.core= as.numeric(cores))
    
        data3 <- set.threshold(data.loco.scan_2,method="M.eff",level=0.05)
        #png(file=paste("GWAS_by_year_Manplot_",colnames(data@pheno[i]),model[j],".png",sep=""),width=800, height=550)#change directory in your case
        manhattan_plot_list[[model[j]]] <- manhattan.plot(data3,traits=colnames(data@pheno[i]))+geom_point(size=3)+theme(text = element_text(size = 25),strip.text = element_text(face = "bold"))
        #print(p2)
        #dev.off()         
      }
  
    }

    #Status
    updateProgressBar(session = session, id = "pb_gwas", value = 100, status = "success", title = "Finished")

  })

  #Genomic Diversity analysis
  observeEvent(input$Diversity_Start, {

    #Input variables
    ploidy <- input$ploidy
    genomat_cleaned <- read.csv(input$genotype_file$datapath)

    #Diversity statistics (might want to output a dataframe of stats for each marker and sample (avg depth, maf, af, het))
    
    #Output heterozygosity statistics (Need to be able to account for different ploidy - maybe use if statements, or ignore the max min over a range of numbers length ploidy?)
    
    #Can I do a range for sum? like sum(genomat_cleaned == c(1:(ploidy-1)))
    if (ploidy == 2) {
      pop_heterozygosity <- sum(sum(genomat_cleaned == 1)) / length(genomat_cleaned)
    } else if (ploidy == 4) {
      pop_heterozygosity <- sum(sum(genomat_cleaned == 1) + sum(genomat_cleaned == 2) + sum(genomat_cleaned == 3)) / length(genomat_cleaned)
    } else if (ploidy == 6) {
      pop_heterozygosity <- sum(sum(genomat_cleaned == 1) + sum(genomat_cleaned == 2) + sum(genomat_cleaned == 3) + sum(genomat_cleaned == 4) + sum(genomat_cleaned == 5)) / length(genomat_cleaned)
    } else if (ploidy == 8) {
      pop_heterozygosity <- sum(sum(genomat_cleaned == 1) + sum(genomat_cleaned == 2) + sum(genomat_cleaned == 3) + sum(genomat_cleaned == 4) + sum(genomat_cleaned == 5) + sum(genomat_cleaned == 6) + sum(genomat_cleaned == 7)) / length(genomat_cleaned)
    } else if (ploidy > 8) {
      print("Sorry, you will have to manually edit R script to estimate heterozygosity for ploidy > 8")
    }
    
    if (exists("pop_heterozygosity")) {
      cat("Population heterozygosity =", pop_heterozygosity, "\n")
    }
    #print('Population heterozygosity =',sum(sum(genomat_cleaned == 1)) / length(genomat_cleaned))
    
    
    ######Get MAF plot
    # Copy the dataframe (assuming you have 'updog_df' loaded)
    df <- genomat_cleaned
    
    # Switch values so that 0 = homozygous reference
    #df[df == 0] <- 2
    #df[df == 2] <- 0
    
    # Calculate the allele frequency values
    result_values <- apply(df, 1, function(row) sum(row) / (as.numeric(ploidy) * ncol(df)))
    
    # Calculate the row-wise sum
    row_sums <- rowSums(df)
    
    # Calculate the length of each row multiplied by ploidy
    row_lengths <- ncol(df) * as.numeric(ploidy)
    
    # Divide the row sums by the row lengths
    result <- row_sums / row_lengths
    
    # Convert the result to allele frequency list
    result_list <- as.vector(result)
    
    # Get minor allele frequency values
    maf_list <- ifelse(result_list < 0.5, result_list, 1 - result_list)
    
    # Create a histogram
    maf_name <- paste0(output_name,'_MADC_Ref_Alt_filtered_bias0.5-2_updog_SNPs_MAF_histogram.pdf')
    pdf(maf_name)
    hist(maf_list, breaks = 25, col = "lightblue", border = "black", xlab = "Minor Allele Frequency (MAF)",
         ylab = "Frequency", main = "MAF Distribution")
    
    dev.off()
    
    
    ######Make dosage ratio plot
    
    # Calculate the value counts and percentages for each DataFrame (assuming you have 'updog_df' and 'geno_df' loaded)
    value_counts1 <- table(unlist(genomat_cleaned))
    total_values1 <- length(unlist(genomat_cleaned))
    percentages1 <- (value_counts1 / total_values1) * 100
    
    
    # Create a bar chart
    bar_width <- 0.35
    dosage_calls <- unique(names(percentages1))
    dosage_counts1 <- ifelse(names(percentages1) %in% dosage_calls, percentages1, 0)
    
    bar_data <- data.frame(DosageCall = dosage_calls,
                           DArT_Updog = dosage_counts1)
    
    # Plot the bar chart using ggplot2
    
    ggplot(bar_data, aes(x = DosageCall)) +
      geom_bar(aes(y = DArT_Updog), stat = "identity", position = "dodge", fill = "blue", width = bar_width, alpha = 0.7) +
      labs(x = "Dosage Call", y = "Percentage", title = "MADC Updog Alt Ref Dosage Call Ratios") +
      theme_minimal() +
      scale_x_discrete(labels = dosage_calls) +
      scale_fill_manual(values = c("blue", "red")) +
      theme_minimal() +
      scale_y_continuous(limits = c(0, max(max(dosage_counts1)) + 10))
    
    # Display the plot
    dosage_name <- paste0(output_name,'_MADC_Ref_Alt_filtered_bias0.5-2_updog_SNPs_dosage_ratios.png')
    ggsave(dosage_name, width = 8, height = 6, dpi = 300)  # Save the plot as an image
   

  })

  output$downloadData <- downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      paste0(input$output_name,'_selected_files.zip')
    },
    content = function(file) {
      # Create a zip file
      zip(file, files = NULL)
      
      # Add selected files to the zip
      filenames <- input$files_to_download
      
      for (filename in filenames) {
        # Write each file to the zip
        write.csv(get(filename), file.path(dirname(file), paste0(input$output_name, "_", filename, ".csv")))
      }
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)