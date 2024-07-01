# Install and load required packages
required_cran_packages <- c("updog", "ggplot2","devtools","GWASpoly","SNPRelate",
                       "adegenet", "future", "scales", "AGHmatrix", "stats", 
                       "factoextra", "readxl", "ggrepel", "dplyr", "shiny",
                       "shinydashboard","randomcoloR","plotly", "DT","RColorBrewer",
                       "dichromat", "bs4Dash", "shinyWidgets","data.table",
                       "matrixcalc","Matrix", "shinyalert","rrBLUP", "tidyverse",
                       "foreach", "doParallel","VariantAnnotation", "vcfR")

required_bio_packages <- c("SNPRelate","VariantAnnotation")

Dev_tools_packages <- c(
  "GWASpoly" = "jendelman/GWASpoly",
  "BIGr" = "Breeding-Insight/BIGr::latest"
  )

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#Bioconductor
for(package in required_bio_packages){
  if(!require(package, character.only = TRUE)) {
    BiocManager::install(package)
    library(package, character.only = TRUE)
  }
}

#CRAN
for(package in required_cran_packages) {
  if(!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

#GitHub
for (package in names(Dev_tools_packages)) {
  if (!require(package, character.only = TRUE)) {
    devtools::install_github(Dev_tools_packages[package], build_vignettes = FALSE)
    library(package, character.only = TRUE)
  }
}

# UI
ui <- dashboardPage(
  skin = "black",
  dashboardHeader(title = tagList(
    tags$img(src = 'BIG_R_logo.png', height = '40', width = '50'),
    tags$span("Breeding Insight Genomics", style = "font-size: 12px; margin-left: 1px;")
  )#,
    #dropdownMenu(type = "notifications",
    #  notificationItem(
    #    text = "5 new users today",
    #    icon("users")
    #  ),
    #  notificationItem(
    #    text = "12 items delivered",
    #    icon("truck"),
    #    status = "success"
    #  ),
    #  notificationItem(
    #    text = "Server load at 86%",
    #    icon = icon("exclamation-triangle"),
    #    status = "warning"
    #  )
   # )
),
  dashboardSidebar(
    skin="light", status = "info",
    sidebarMenu(
      flat = FALSE,
      menuItem("Home", tabName = "welcome", icon = icon("house")),
      menuItem("Dosage Calling", tabName = "dosage_calling", icon = icon("diagram-next"),
              menuSubItem("DArT Report2VCF", tabName = "dosage2vcf", icon = icon("share-from-square")),
              menuSubItem("Updog Dosage Calling", tabName = "updog", icon = icon("list-ol")),
              menuSubItem("VCF Filtering", tabName = "filtering", icon = icon("filter"))),
      menuItem("Population Structure", tabName = "pop_struct", icon = icon("layer-group"),
              menuSubItem("PCA", tabName = "pca", icon = icon("chart-simple")),
              menuSubItem("DAPC", tabName = "dapc", icon = icon("circle-nodes"))),
      menuItem("Genomic Diversity", tabName = "diversity", icon = icon("chart-pie")),
      menuItem("GWAS", tabName = "gwas", icon = icon("think-peaks")),
      #menuItem("QTL Analysis", tabName = "qtl", icon = icon("chart-area")),
      menuItem(
        span("Genomic Prediction", bs4Badge("beta", position = "right", color = "success")),
        tabName = "prediction", 
        icon = icon("right-left")),
      menuItem("Source Code", icon = icon("circle-info"), href = "https://www.github.com/Breeding-Insight/Genomics_Shiny_App"),
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
              "The app is currently under development",
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
                fileInput("updog_rdata","Choose VCF File", accept = c(".vcf",".vcf.gz")),
                textInput("filter_output_name", "Output File Name"),
                numericInput("filter_ploidy","Ploidy", min = 0, value = NULL),
                numericInput("filter_maf","MAF filter", min = 0, max=1, value = 0.05, step = 0.01),
                sliderInput("size_depth","Minimum Read Depth", min = 0, max = 300, value = 10, step = 1),
                numericInput("snp_miss","Remove SNPs with >= % missing data", min = 0, max = 1, value = 0.5, step = 0.1),
                numericInput("sample_miss","Remove Samples with >= % missing data", min = 0, max = 1, value = 0.5, step = 0.1),
                "Updog Filtering Parameters",
                checkboxInput("use_updog", "Use Updog Filtering Parameters?", value = FALSE),
                conditionalPanel(
                  condition = "input.use_updog == true",
                  div(
                    numericInput("OD_filter", "Max OD (Updog filter)", min = 0, value = 0.05),
                    sliderInput("Bias", "Bias (Updog filter)", min = 0, max = 10, value = c(0.5, 2), step = 0.1),
                    numericInput("Prop_mis", "Max Prop_mis (Updog filter)", min = 0, max = 1, value = 0.05, step = 0.05),
                    numericInput("maxpostprob_filter", "Minimum maxpostprob (Updog filter)", min = 0, value = 0.5, step = 0.1)
                  )
                ),
                downloadButton("start_updog_filter", "Download Filtered VCF", icon = icon("download")),
                  div(style="display:inline-block; float:right",dropdownButton(
                    tags$h3("Updog Filter Parameters"),
                    #selectInput(inputId = 'xcol', label = 'X Variable', choices = names(iris)),
                    #selectInput(inputId = 'ycol', label = 'Y Variable', choices = names(iris), selected = names(iris)[[2]]),
                    #sliderInput(inputId = 'clusters', label = 'Cluster count', value = 3, min = 1, max = 9),
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
              tabPanel("Bias Histogram", icon = icon("image"), plotOutput("bias_hist", height = '550px')),
              tabPanel("OD Histogram", icon = icon("image"), plotOutput("od_hist", height = '550px')),
              tabPanel("Prop_mis Histogram", icon = icon("image"), plotOutput("maxpostprob_hist", height = '550px'))
              #tabPanel("ReadDepth Histogram", icon = icon("image"), plotOutput("depth_hist", height = '550px'))
              #tabPanel("SNP Distribution Plot", icon = icon("image"), plotOutput("snp_dist", height = '550px')),
              #tabPanel("SNP % Missing Histogram", icon = icon("image"), plotOutput("missing_snp_hist", height = '550px')),
              #tabPanel("Sample % Missing Histogram", icon = icon("image"), plotOutput("missing_sample_hist", height = '550px')),
              #tabPanel("Summary Statistics", icon = icon("sliders"), tableOutput("dosages"))
              #plotOutput("coverage"), # Placeholder for plot outputs
              )
            ),
          column(width = 3,
            valueBoxOutput("snp_retained_box", width = NULL),
            valueBoxOutput("snp_removed_box", width = NULL),
            #valueBox("0%","SNPs Removed", icon = icon("filter-circle-xmark"), width = NULL, color = "info"), #https://rstudio.github.io/shinydashboard/structure.html#tabbox
            box(title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE,
                sliderInput("hist_bins","Histogram Bins", min = 1, max = 1200, value = c(50), step = 1), width = NULL,
                div(style="display:inline-block; float:left",dropdownButton(
                      tags$h3("Save Image"),
                      selectInput(inputId = 'filter_hist', label = 'Figure', choices = c("Bias Histogram", 
                                                                                         "OD Histogram", 
                                                                                         "Prop_mis Histogram")),
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
            ),
            box(title = "Status", width =12, collapsible = TRUE, status = "info",
                progressBar(id = "pb_filter", value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
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
              fileInput("madc_file", "Choose MADC or VCF File", accept = c(".csv",".vcf",".vcf.gz")),
              #checkboxInput("off-targets","Include off-target loci?"),
              #fileInput("sample_file", "Optional: Choose Sample List (disabled)", accept = c(".csv")),
              textInput("output_name", "Output File Name"),
              selectInput("markers", "Select Markers", choices = c("All Loci (not supported)", "Target Loci Only"), selected = "Target Loci Only"),
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
        tabName = "dosage2vcf",
        fluidPage(
          fluidRow(
              box(
              title = "Inputs", status = "info", solidHeader = TRUE, collapsible = FALSE, collapsed = FALSE,
              fileInput("report_file", "Choose DArT Dose Report File", accept = c(".csv")),
              fileInput("counts_file", "Choose DArT Counts File", accept = c(".csv")),
              #checkboxInput("off-targets","Include off-target loci?"),
              #fileInput("sample_file", "Optional: Choose Sample List (disabled)", accept = c(".csv")),
              textInput("d2v_output_name", "Output File Name"),
              numericInput("dosage2vcf_ploidy", "Species Ploidy", min = 1, value = 2),
              downloadButton("download_d2vcf", "Download VCF File"),
              div(style="display:inline-block; float:right",dropdownButton(

                    tags$h3("DArT File Converstion"),
                    "Converting DArT report files to VCF format. The VCF file will automatically
                    download when complete.",
                    circle = FALSE,
                    status = "warning", 
                    icon = icon("info"), width = "300px",
                    tooltip = tooltipOptions(title = "Click to see info!")
                ))
              ),
              valueBoxOutput("ReportSnps")
              #valueBox("Help","Updog Manual", icon = icon("globe"), color = "warning")
          ),
         
         fluidRow(
              box(title = "Status", width = 3, collapsible = TRUE, status = "info",
                progressBar(id = "dosage2vcf_pb", value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
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
              fileInput("dosage_file", "Choose Genotypes File"),#, accept = c(".csv",".vcf",".vcf.gz")),
              fileInput("passport_file", "Choose Passport File (Sample IDs in first column)", accept = c(".csv")),
              #textInput("output_name", "Output File Name (disabled)"),
              #Dropdown will update after pasport upload
              numericInput("pca_ploidy", "Species Ploidy", min = 1, value = 2),
              actionButton("pca_start", "Run Analysis"),
              div(style="display:inline-block; float:right",dropdownButton(

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
                      selectInput(inputId = 'pca_figure', label = 'Figure', choices = c("2D Plot", "Scree Plot"), selected = "2D Plot"),
                      selectInput(inputId = 'pca_image_type', label = 'File', choices = c("jpeg","tiff","png"), selected = "jpeg"),
                      sliderInput(inputId = 'pca_image_res', label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                      sliderInput(inputId = 'pca_image_width', label = 'Width', value = 3, min = 1, max = 10, step=0.5),
                      sliderInput(inputId = 'pca_image_height', label = 'Height', value = 3, min = 1, max = 10, step = 0.5),
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
          column(width = 1)


        )
      ),
      tabItem(
        tabName = "dapc",
        # Add PCA content here
        fluidRow(
          column(width = 3,
            box(
              title = "Inputs", width = 12, solidHeader = TRUE, status = "info",
              tabsetPanel(
                tabPanel("Step 1:(K)", width = 12,
                  fileInput("dosage_file", "Choose Genotypes File", accept = c(".csv",".vcf",".vcf.gz")),
                  #fileInput("passport_file", "Choose Passport File (Sample IDs in first column)", accept = c(".csv")),
                  #textInput("output_name", "Output File Name (disabled)"),
                  #Dropdown will update after pasport upload
                  numericInput("dapc_kmax", "Maximum K", min = 1, value = 5),
                  numericInput("dapc_ploidy", "Species Ploidy", min = 1, value = 2),
                  actionButton("K_start", "Run Analysis"),
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
                  fileInput("dosage_file", "Choose Genotypes File", accept = c(".csv",".vcf",".vcf.gz")),
                  #fileInput("passport_file", "Choose Passport File (Sample IDs in first column)", accept = c(".csv")),
                  #textInput("output_name", "Output File Name (disabled)"),
                  #Dropdown will update after pasport upload
                  numericInput("dapc_k", "Number of Clusters (K)", min = 1, value = NULL),
                  numericInput("dapc_ploidy", "Species Ploidy", min = 1, value = 2),
                  actionButton("dapc_start", "Run Analysis"),
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


              ),

              style = "overflow-y: auto; height: 420px"
            ),

            box(
                title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE,collapsed = FALSE, width = 12,
                "Change the plot parameters", br(),
                #selectInput('group_info', label = 'Color Variable (eg, Taxon)', choices = NULL),
                selectInput("color_choice", "Color Palette", choices = c("YlOrRd","YlOrBr","YlGnBu","YlGn",
                                                                                  "Reds","RdPu","Purples","PuRd","PuBuGn","PuBu",
                                                                                  "OrRd","Oranges","Greys","Greens","GnBu","BuPu",
                                                                                  "BuGn","Blues","Set3","Set2","Set1","Pastel2",
                                                                                  "Pastel1","Paired","Dark2","Accent","Spectral",
                                                                                  "RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn",
                                                                                  "PiYG","BrBG"), selected = "Paired"),
                materialSwitch('plot_BICX', label = "Label Suggested K?", status = "success", value = TRUE),
                #selectInput("pc_X", "X-Axis (2D-Plot only)", choices = c("PC1","PC2","PC3","PC4","PC5"), selected = "PC1"),
                #selectInput("pc_Y", "Y-Axis (2D-Plot only)", choices = c("PC1","PC2","PC3","PC4","PC5"), selected = "PC2"),
                div(style="display:inline-block; float:right",dropdownButton(
                      tags$h3("Save Output"),
                      selectInput(inputId = 'dapc_figure', label = 'Figure', choices = c("BIC Plot", "DAPC Plot"), selected = "BIC Plot"),
                      selectInput(inputId = 'dapc_image_type', label = 'File', choices = c("jpeg","tiff","png"), selected = "jpeg"),
                      sliderInput(inputId = 'dapc_image_res', label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                      sliderInput(inputId = 'dapc_image_width', label = 'Width', value = 3, min = 1, max = 10, step=0.5),
                      sliderInput(inputId = 'dapc_image_height', label = 'Height', value = 3, min = 1, max = 10, step = 0.5),
                      fluidRow(
                      downloadButton("download_dapc_image", "Save Image"),
                      downloadButton("download_dapc_file", "Save Files")),
                      circle = FALSE,
                      status = "danger", 
                      icon = icon("floppy-disk"), width = "300px",
                      tooltip = tooltipOptions(title = "Click to see inputs!")
                    ))
              ),

              box(title = "Status", width = 12, collapsible = TRUE, status = "info",

                progressBar(id = "pb_dapc", value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
              )
          ),
          column(width = 8,
            box(title = "DAPC Data", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "info", collapsed = FALSE,
                tabsetPanel(
                  tabPanel("BIC Values",DTOutput('BIC_table')),
                  tabPanel("DAPC Values", DTOutput('DAPC_table'))), # Placeholder for plot outputs
                  style = "overflow-y: auto; height: 420px"
            ),

            box(title = "DAPC Plots", status = "info", solidHeader = FALSE, width = 12, height = 550,
                tabsetPanel(
                  tabPanel("BIC Plot",plotOutput("BIC_plot", height = '460px')),
                  tabPanel("DAPC Plot", plotOutput("DAPC_plot", height = '460px')),
                  tabPanel("STRUCTURE Plot", "Not yet supported"))
                  #tabPanel("STRUCTURE Plot", plotOutput("STRUCTURE_plot", height = '460px')))
              )


          ),
          column(width = 1)


        )
      ),
      tabItem(
        tabName = "gwas",
        # Add GWAS content here
        fluidRow(
          column(width = 3,
            box(title="Inputs", width = 12, collapsible = TRUE, collapsed = FALSE, status = "info", solidHeader = TRUE,
              fileInput("gwas_file", "Choose Genotypes File", accept = c(".csv",".vcf",".gz")),
              fileInput("phenotype_file", "Choose Phenotype File", accept = ".csv"),
              #textInput("output_name", "Output File Name"),
              numericInput("gwas_ploidy", "Species Ploidy", min = 1, value = 2),
              selectInput('gwas_threshold', label='Significance Threshold Method', choices = c("M.eff","Bonferroni","FDR","permute"), selected="M.eff"),
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
              fileInput("diversity_file", "Choose Genotypes File", accept = c(".csv",".vcf",".vcf.gz")),
              #fileInput("pop_file", "Choose Passport File"),
              #textInput("output_name", "Output File Name"),
              numericInput("diversity_ploidy", "Species Ploidy", min = 1, value = 2),
              selectInput("zero_value", "What are the Dosage Calls?", choices = c("Reference Allele Counts", "Alternate Allele Counts"), selected = NULL),
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

              ),
              box(title = "Plot Controls", width=12, status = "warning", solidHeader = TRUE, collapsible = TRUE,
                sliderInput("hist_bins","Histogram Bins", min = 1, max = 200, value = c(20), step = 1),
                div(style="display:inline-block; float:left",dropdownButton(
                      tags$h3("Save Image"),
                      selectInput(inputId = 'div_figure', label = 'Figure', choices = c("Dosage Plot",
                                                                                        "AF Histogram", 
                                                                                         "MAF Histogram", 
                                                                                         "OHet Histogram")),
                      selectInput(inputId = 'div_image_type', label = 'File Type', choices = c("jpeg","pdf","tiff","png"), selected = "jpeg"),
                      sliderInput(inputId = 'div_image_res', label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                      sliderInput(inputId = 'div_image_width', label = 'Width', value = 3, min = 1, max = 10, step=0.5),
                      sliderInput(inputId = 'div_image_height', label = 'Height', value = 3, min = 1, max = 10, step = 0.5),
                      fluidRow(
                      downloadButton("download_div_figure", "Save Image"),
                      downloadButton("download_div_file", "Save Files")),
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
              tabsetPanel(
                tabPanel("Dosage Plot", plotOutput('dosage_plot'),style = "overflow-y: auto; height: 500px"),
                tabPanel("AF Plot", plotOutput('af_plot'),style = "overflow-y: auto; height: 500px"),
                tabPanel("MAF Plot", plotOutput('maf_plot'),style = "overflow-y: auto; height: 500px"),
                tabPanel("OHet Plot", plotOutput('het_plot'),style = "overflow-y: auto; height: 500px"),
                tabPanel("Sample Table",DTOutput('sample_table'),style = "overflow-y: auto; height: 470px"),
                tabPanel("SNP Table",DTOutput('snp_table'),style = "overflow-y: auto; height: 470px")


                )
              
              )

            ),
          column(width = 3,
            #valueBoxOutput("snps"),
            valueBoxOutput("mean_het_box", width = NULL),
            valueBoxOutput("mean_maf_box", width = NULL),
            box(title = "Status", width = 12, collapsible = TRUE, status = "info",
              progressBar(id = "pb_diversity", value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
            )
            #valueBoxOutput("mean_pic_box", width = NULL),
            #valueBox(0,"Mean Heterozygosity", icon = icon("dna"), width = NULL, color = "info"),
            #valueBox(0,"Mean Minor-Allele-Frequency", icon = icon("dna"), width = NULL, color = "info"), #https://rstudio.github.io/shinydashboard/structure.html#tabbox
            #valueBox(0,"Mean PIC", icon = icon("dna"), width = NULL, color = "info")
          
          )
          
          )
      
      ),
      tabItem(
        tabName = "prediction",
        # Add GWAS content here
        fluidRow(
          column(width = 3,
            box(title="Inputs", width = 12, collapsible = TRUE, collapsed = FALSE, status = "info", solidHeader = TRUE,
              fileInput("pred_file", "Choose Genotypes File", accept = ".csv"),
              fileInput("trait_file", "Choose Passport File", accept = ".csv"),
              #textInput("output_name", "Output File Name"),
              numericInput("pred_ploidy", "Species Ploidy", min = 1, value = 2),
              numericInput("pred_cv", "Cross-Validations", min = 1, value = 10),
              numericInput("pred_train", "Training Subset Size (%)", min = 1, max = 99, value = 60),
              #selectInput('pred_trait_info', label = 'Select Trait (eg, Color):', choices = NULL),
              virtualSelectInput(
                inputId = "pred_trait_info",
                label = "Select Trait (eg, Color):",
                choices = NULL,
                showValueAsTags = TRUE,
                search = TRUE,
                multiple = TRUE
              ),
              virtualSelectInput(
                inputId = "pred_fixed_info",
                label = "Select Fixed Effects (eg, Group):",
                choices = NULL,
                showValueAsTags = TRUE,
                search = TRUE,
                multiple = TRUE
              ),
              sliderInput("pred_cores", "Number of CPU Cores", min = 1, max = (future::availableCores() - 1), value = 1, step = 1),
              actionButton("prediction_start", "Run Analysis"),
              #downloadButton("download_pca", "Download All Files"),
              #plotOutput("pca_plot"), # Placeholder for plot outputs
              #checkboxGroupInput("files_to_download", "Select files to download:",
                             #choices = c("PC1vPC2 plot", "PC2vPC3 plot"), selected = c("table1", "table2"))
              div(style="display:inline-block; float:right",dropdownButton(
                    tags$h3("GP Parameters"),
                    #selectInput(inputId = 'xcol', label = 'X Variable', choices = names(iris)),
                    #selectInput(inputId = 'ycol', label = 'Y Variable', choices = names(iris), selected = names(iris)[[2]]),
                    #sliderInput(inputId = 'clusters', label = 'Cluster count', value = 3, min = 1, max = 9),
                    "GP uses the rrBLUP package: It can impute missing data, maybe adapt to different ploidy, perform cross validations, define training size, run multiple traits, and accept multiple fixed effects.",
                    circle = FALSE,
                    status = "warning", 
                    icon = icon("info"), width = "300px",
                    tooltip = tooltipOptions(title = "Click to see info!")
                ))#,
              #style = "overflow-y: auto; height: 550px"
 
              ),
            box(title = "Status", width = 12, collapsible = TRUE, status = "info",
 
                progressBar(id = "pb_prediction", value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
              )
            ),

          column(width = 6,
            box(
              title = "Plots", status = "info", solidHeader = FALSE, width = 12, height = 600,
              tabsetPanel(
                tabPanel("Violin Plot", plotOutput("pred_violin_plot", height = "500px")),
                tabPanel("Box Plot", plotOutput("pred_box_plot", height = "500px")),
                tabPanel("Heatmap Plot", plotOutput("pred_heat_plot", height = "500px")),
                tabPanel("Results Table", DTOutput("pred_table"),style = "overflow-y: auto; height: 500px")

                )
              
              )

            ),
          
          column(width = 3,
            valueBox(0,"Samples in Genotype File", icon = icon("dna"), width = NULL, color = "info"),
            valueBox(0,"Samples with Phenotype Information", icon = icon("dna"), width = NULL, color = "info"),
            #valueBox("0","QTLs Detected", icon = icon("dna"), width = NULL, color = "info"), #https://rstudio.github.io/shinydashboard/structure.html#tabbox
            box(title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
                #sliderInput("hist_bins","Histogram Bins", min = 1, max = 1200, value = c(50), step = 1), width = NULL,
                selectInput("pred_trait_select", label = "Trait Selection", choices = c("all","1-dom","2-dom","additive","general","diplo-general","diplo-additive")),
                div(style="display:inline-block; float:left",dropdownButton(
                      tags$h3("Save Image"),
                      selectInput(inputId = 'pred_figures', label = 'Figure', choices = c("Scatter Plot", 
                                                                                         "Box Plot", 
                                                                                         "Heatmap Plot",
                                                                                         "Blank")),
                      selectInput(inputId = 'pred_image_type', label = 'File Type', choices = c("jpeg","pdf","tiff","png"), selected = "jpeg"),
                      sliderInput(inputId = 'pred_image_res', label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                      sliderInput(inputId = 'pred_image_width', label = 'Width', value = 3, min = 1, max = 10, step=0.5),
                      sliderInput(inputId = 'pred_image_height', label = 'Height', value = 3, min = 1, max = 10, step = 0.5),
                      downloadButton("download_pred_table", "Save"),
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
        tabName = "code",
        fluidPage(
          # Add source code here
        )
      ),
      tabItem(
        tabName = "help",
        fluidPage(
          column(width=12),
          column(width=12,
            box(title="Dosage Calling", width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
              tabsetPanel(
                tabPanel("DArT Report2VCF",
                  "**Draft**This tab is designed to convert the DArT Dose Report and Counts files to a VCF file. **DArT Website**"
                ),
                tabPanel("Updog Dosage Calling",
                  "**Draft**This tab is designed to handle the process of dosage calling in genomic data. Dosage calling is essential for determining the number of copies of a particular allele at each genomic location. The app likely includes functionalities to upload raw genomic data, apply various filtering criteria, and generate plots to visualize the distribution of dosages. Users can examine histograms for SNP max post probabilities and read depths, which help in assessing the quality and accuracy of the dosage calls.**Updog**"
                ),
                tabPanel("SNP Filtering",
                  "Filtering the genotypes"
                ))
            ),
            box(title="Population Structure", width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
              tabsetPanel(
                tabPanel("PCA",
                  "**Draft**This tab focuses on analyzing the population structure using Discriminant Analysis of Principal Components (DAPC) and Principal Component Analysis (PCA). These methods are used to identify and visualize genetic diversity and structure within the population. The app provides options to perform PCA to reduce the dimensionality of the genomic data and visualize principal components. DAPC is used to find clusters within the data and visualize these clusters, helping users understand the genetic relationships and structure in their dataset."
                ),
                tabPanel("DAPC"))

            ),
            box(title="Genomic Diversity", width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
              "**Draft**This tab is dedicated to analyzing genomic diversity within the population. It calculates various diversity metrics such as heterozygosity and minor allele frequency (MAF). The app includes functionalities to visualize these metrics through histograms and other plots. Users can download the calculated diversity metrics as CSV files. This tab helps in understanding the genetic variability and distribution of alleles within the population."
            ),
            box(title="GWAS", width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
              "**Draft**This tab is for conducting Genome-Wide Association Studies (GWAS) to identify associations between genetic variants and traits of interest. Users can input phenotypic data and specify parameters for the GWAS analysis. The app performs statistical tests to identify significant associations between SNPs and traits, and visualizes the results using Manhattan plots and Q-Q plots. This tab helps in identifying potential genetic markers linked to specific traits.**List R packages utilized"
            ),
            box(title="Genomic Prediction/Selection", width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
              "**Draft**This tab provides functionalities for genomic prediction, which involves predicting phenotypic traits based on genomic data. Users can input phenotypic and genotypic data, and specify parameters such as the number of cross-validation folds, training percentage, and fixed effects. The app performs genomic prediction using methods such as rrBLUP, and displays the results including cross-validation performance metrics. Users can download the prediction results for further analysis.**List R packages utilized"
            ),
            box(title="How to Cite", width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
              "**Draft**Instructions for citing the app and packages used in analyses"
            ),
          ),
          column(width=2)
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
  
  #Plots for downloading
  all_plots <- reactiveValues(
    pca_2d = NULL,
    pca_3d = NULL,
    pca_scree = NULL

    )
  
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

  report_snp_number <- reactiveVal(0)

  #SNP counts value box
  output$ReportSnps <- renderValueBox({
    valueBox(snp_number(), "Markers in VCF File", icon = icon("dna"), color = "info")
    })

##This is for the DArT files conversion to VCF
  output$download_d2vcf <- downloadHandler(
    filename = function() {
      paste0(input$d2v_output_name, ".vcf.gz")
    },
    content = function(file) {
      # Ensure the files are uploaded
      req(input$report_file, input$counts_file, input$d2v_output_name, input$dosage2vcf_ploidy)
      
      # Get the uploaded file paths
      dosage_file <- input$report_file$datapath
      counts_file <- input$counts_file$datapath
      ploidy <- input$dosage2vcf_ploidy
      
      # Use a temporary file path without appending .vcf
      temp_base <- tempfile()
      
      #Status
      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 50, title = "Converting DArT files to VCF")

      # Convert to VCF using the BIGr package
      cat("Running BIGr::dosage2vcf...\n")
      BIGr::dosage2vcf(
        dart.report = dosage_file,
        dart.counts = counts_file,
        output.file = temp_base,
        ploidy = as.numeric(ploidy)
      )
      
      # The output file should be temp_base.vcf
      output_name <- paste0(temp_base, ".vcf")
      
      # Check if the VCF file was created
      if (file.exists(output_name)) {
        cat("VCF file created successfully.\n")
        
        # Compress the VCF file using gzip
        gzip_file <- paste0(output_name, ".gz")
        gz <- gzfile(gzip_file, "w")
        writeLines(readLines(output_name), gz)
        close(gz)
        
        # Check if the gzip file was created
        if (file.exists(gzip_file)) {
          cat("Gzip file created successfully.\n")
          
          # Move the compressed file to the path specified by 'file'
          file.copy(gzip_file, file)
          
          # Delete the temporary files
          unlink(gzip_file)
          unlink(output_name)
          
          cat("Temporary files deleted successfully.\n")
        } else {
          stop("Error: Failed to create the gzip file.")
        }
      } else {
        stop("Error: Failed to create the VCF file.")
      }

      #Status
      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 100, title = "Complete! - Downloading VCF")
    }
  )

##This is for performing Updog Dosage Calling
  observeEvent(input$run_analysis, {
    # Get inputs
    madc_file <- input$madc_file$datapath
    sample_file <- input$sample_file$datapath
    output_name <- input$output_name
    ploidy <- input$ploidy
    cores <- input$cores
    model_select <- input$updog_model
    marker_set <- (markers == "Target Loci Only")
    

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
      filtered_df <- madc_df[!grepl("\\|AltMatch|\\|RefMatch", madc_df$AlleleID), ]

      #Remove extra text after Ref and Alt (_001 or _002)
      filtered_df$AlleleID <- sub("\\|Ref.*", "|Ref", filtered_df$AlleleID)
      filtered_df$AlleleID <- sub("\\|Alt.*", "|Alt", filtered_df$AlleleID)
      
      # Save the csv file for review and use in R
      #df_name <- paste0(output_name,'_MADC_alt_ref_counts.csv')
      
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
      
      #remove alt or ref rows that do not have a counterpart in the other dataframe
      if (nrow(ref_df) > nrow(alt_df)) {
        ref_df <- ref_df[ref_df$CloneID %in% alt_df$CloneID,]
      } else if (nrow(ref_df) < nrow(alt_df)) {
        alt_df <- alt_df[alt_df$CloneID %in% ref_df$CloneID,]
      } else {
        alt_df <- alt_df[alt_df$CloneID %in% ref_df$CloneID,]
      }

      #Ensure that each has the same SNPs and that they are in the same order
      identical(alt_df$CloneID,ref_df$CloneID)
      
      ###Convert the ref and alt counts into matrices with the CloneID as the index
      #Set SNP names as index
      row.names(ref_df) <- ref_df$CloneID
      row.names(alt_df) <- alt_df$CloneID
      
      #Remove unwanted columns and convert to matrix
      #Probably best to just remove the column names that aren't wanted instead of the first 16 columns.
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
  
  #vcf
  filtering_files <- reactiveValues(
      raw_vcf_df = NULL

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
      snp_miss <- input$snp_miss
      sample_miss <- input$sample_miss
      ploidy <- as.numeric(input$filter_ploidy)
      maf_filter <- input$filter_maf

      
      temp_file <- tempfile(fileext = ".vcf.gz")

      #Status
      updateProgressBar(session = session, id = "dosage2vcf_pb", value = 50, title = "Converting DArT files to VCF")

      # Convert to VCF using the BIGr package
      cat("Running BIGr::dosage2vcf...\n")
      updateProgressBar(session = session, id = "pb_filter", value = 10, title = "Processing VCF file")

      #Input file
      vcf <- vcfR::read.vcfR(input$updog_rdata$datapath)
      #Starting SNPs
      starting_snps <- nrow(vcf)
      #export INFO dataframe
      filtering_files$raw_vcf_df <- data.frame(vcf@fix)

      #Pb
      updateProgressBar(session = session, id = "pb_filter", value = 40, title = "Filtering VCF file")

      #Filtering
      vcf <- BIGr::filterVCF(vcf.file = vcf,
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
      
      #Pb
      updateProgressBar(session = session, id = "pb_filter", value = 70, title = "Exporting Filtered VCF")

      #Writing file
      vcfR::write.vcf(vcf, file = temp_file)

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

  ##Updog file stats
  #Consider Extracting the GT info or UD info if present as a datafrfame,
  #Obtaining the info in the INFO column as it's own dataframe with a column for each value
  #Then remove the VCF file and use the remaining dataframes for producing the figures
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

      ##Read Depth (I would prefer that this show the mean depth for SNPs or Samples instead of all loci/sample cells)

      quantile(as.numeric(new_df$DP), 0.95)
      #Histogram
      output$depth_hist <- renderPlot({

        hist(as.numeric(new_df$DP), 
          main = "Unfiltered SNP Total Read Depth Across All Samples",
          xlab = "Total Read Depth per SNP",
          ylab = "Genomic Sites",
          col = "lightblue",
          border = "black",
          xlim = c(0,1000),
          breaks = as.numeric(input$hist_bins))
        axis(1, at = seq(0, 1000, by = 20), labels = rep("", length(seq(0, 1000, by = 20))))  # Add ticks
        abline(v = 0.05, col = "black", lty = 2)  # proposed filter by updog

        # Add vertical lines
        abline(v = mean(as.numeric(new_df$DP)), col = "red", lty = 2)  # Mean line
        abline(v = median(as.numeric(new_df$DP)), col = "green", lty = 2)  # Median line
        #abline(v = 0.05, col = "black", lty = 2)  # proposed filter by updog

      })
  
  })
  
  #PCA dropdown

  data <- reactiveValues(info_df = NULL)

  # Update dropdown menu choices based on uploaded passport file
  observeEvent(input$passport_file, {
    info_df <- read.csv(input$passport_file$datapath, header = TRUE, check.names = FALSE)
    info_df[,1] <- as.character(info_df[,1]) #Makes sure that the sample names are characters instead of numeric
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
    req(input$pca_ploidy)
    # Get inputs
    geno <- input$dosage_file$datapath
    pedigree_df <- input$passport_file$datapath
    g_info <- as.character(input$group_info)
    output_name <- input$output_name
    ploidy <- input$pca_ploidy
    
    PCX <- input$pc_X
    PCY <- input$pc_Y

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
      genomat <- extract.gt(vcf, element = "UD")
      class(genomat) <- "numeric"
      rm(vcf) #Remove vcf
    }else{
      #Extract GT and convert to numeric calls
      genomat <- extract.gt(vcf, element = "GT")
      genomat <- apply(genomat, 2, convert_to_dosage)
      rm(vcf) #Remove VCF
    }

    #Add support for genotype matrix
    #} else {
      #Import genotype matrix
   #  genomat <- read.csv(geno, header = TRUE, row.names = 1, check.names = FALSE)
   # } 


    #Start analysis

    #Passport info
    # Sample dataframe with a column of taxon names
    info_df <- read.csv(pedigree_df, header = TRUE, check.names = FALSE)

    # Print the modified dataframe
    row.names(info_df) <- info_df[,1]

    #Plotting
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

  ##2D PCA plotting
  observe({
    req(pca_data$pc_df_pop, pca_data$variance_explained, pca_data$my_palette, input$grey_choice)
    
    # Generate colors
    unique_countries <- unique(pca_data$pc_df_pop[[input$group_info]])
    palette <- brewer.pal(length(unique_countries), input$color_choice)
    my_palette <- colorRampPalette(palette)(length(unique_countries))

    # Define a named vector to map input labels to grey values
    label_to_value <- c("Light Grey" = "grey80",
                        "Grey" = "grey60",
                        "Dark Grey" = "grey40",
                        "Black" = "black")

    # Get the corresponding value based on the selected grey
    selected_grey <- label_to_value[[input$grey_choice]]

    # Similar plotting logic here
    if (input$use_cat) {
        # cat plotting logic
      cat_colors <- c(input$cat_color, "grey")
      plot <- ggplot(pca_data$pc_df_pop, aes(x = pca_data$pc_df_pop[[input$pc_X]], 
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
        # non-cat plotting logic
        plot <- ggplot(pca_data$pc_df_pop, aes(x = pca_data$pc_df_pop[[input$pc_X]], 
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

    all_plots$pca_2d <- plot  # Assign the plot to your reactiveValues
  })

  #Plot the 2d plot
  output$pca_plot_ggplot <- renderPlot({
    req(all_plots$pca_2d)  # Ensure the plot is ready
    all_plots$pca_2d
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

  observe({
    #PCA scree plot
    req(pca_data$variance_explained)

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

      all_plots$pca_scree <- plot

    })

  #Scree plot
  output$scree_plot <- renderPlot({
    req(all_plots$pca_scree)  # Ensure the plot is ready
    all_plots$pca_scree
  })

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

  ##GWAS items
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

    # Create a temporary file for the selected phenotype data
    temp_pheno_file <- tempfile(fileext = ".csv")

    #Save new phenotype file with selected traits and fixed effects
    #write.csv(phenotype_file, file = "phenotypes_selected_traits.csv", row.names = FALSE)
    write.csv(phenotype_file, file = temp_pheno_file, row.names = FALSE)

    #Remove the phenotype_file from memory
    rm(phenotype_file)

    #Status
    updateProgressBar(session = session, id = "pb_gwas", value = 5, title = "Upload Complete: Now Formatting GWASpoly Data")

    #Geno file path
    file_path <- input$gwas_file$datapath

    #Geno.file conversion if needed
    if (grepl("\\.csv$", file_path)) {
        data <- GWASpoly::read.GWASpoly(ploidy= ploidy, pheno.file= temp_pheno_file, geno.file=input$gwas_file$datapath,
                          format="numeric", n.traits=length(traits), delim=",") #only need to change files here

    } else if (grepl("\\.vcf$", file_path) || grepl("\\.vcf\\.gz$", file_path)) {
      # Create a temporary file for the selected phenotype data
      temp_geno_file <- tempfile(fileext = ".csv")

      #Function to convert GT to dosage calls (add to BIGr)
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
      
      #Convert VCF file if submitted
      vcf <- vcfR::read.vcfR(input$gwas_file$datapath)

      #Extract GT
      geno_mat <- extract.gt(vcf, element = "GT")
      geno_mat <- apply(geno_mat, 2, convert_to_dosage)
      info <- data.frame(vcf@fix)
      gpoly_df <- cbind(info[,c("ID","CHROM","POS")], geno_mat)
      write.csv(gpoly_df, file = temp_geno_file, row.names = FALSE)

      data <- GWASpoly::read.GWASpoly(ploidy= ploidy, pheno.file= temp_pheno_file, geno.file=temp_geno_file,
                          format="numeric", n.traits=length(traits), delim=",")
    } else {

        # If condition is met, show notification toast
        shinyalert(
            title = "Oops",
            text = "No valid genotype file detected",
            size = "xs",
            closeOnEsc = TRUE,
            closeOnClickOutside = FALSE,
            html = TRUE,
            type = "info",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE,
          )

      #Stop the analysis
      return()
    }

    data.loco <- set.K(data,LOCO=F,n.core= as.numeric(cores))

    #Delete temp pheno file
    unlink(temp_pheno_file)

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
      plotBICs_kinship <- subset(plotBICs,plotBICs$RelationshipMatrix =="w/Kinship")
      #write.table(plotBICs_kinship,paste(colnames(GE)[i],"_model_selection_BIC.txt",sep=""),row.names=F,sep="\t",quote=F)
  
      output$bic_table <- renderDT({plotBICs_kinship}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5)
        )

      #setwd("~/Desktop/Alfalfa_GWAS/GWAS by year/ModelSel") #change directory in your case
      #pdf(paste(colnames(GE)[i],"_model_selection_BIC.pdf",sep=""),height=4,width=7)
      
      output$bic_plot <- renderPlot({
        p1<-ggplot(plotBICs_kinship, aes(x=n.PC, y=BIC,group=RelationshipMatrix)) + 
          geom_line(color="grey")+
          geom_point(shape=21, color="black", fill="#d95f0e", size=3)+
          #geom_point(aes(shape=ISTL1_intro,size=2))+
          theme(text=element_text(size=15),axis.text.x = element_text(angle =0),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))+
          labs(x = "Number of PCs",y="BIC")
        print(p1)

      })
      #dev.off()
      
      #Status
      updateProgressBar(session = session, id = "pb_gwas", value = 40, title = "BIC Complete: Now Performing GWAS")

      ##GWAS based on model selection
      N <- nrow(data@pheno) #Population size
      model <- c("additive","1-dom","2-dom","general","diplo-general","diplo-additive")
  
      BIC_min <- plotBICs_kinship[which.min(plotBICs_kinship$BIC),]
      if(BIC_min$n.PC == 0){params <- set.params(geno.freq = 1 - 5/N)}else{params <- set.params(geno.freq = 1 - 5/N,n.PC = as.numeric(levels(BIC_min$n.PC))[BIC_min$n.PC])}
      data.loco.scan <- GWASpoly(data=data.loco,models=model,traits=colnames(data@pheno[i]),params=params,n.core=9)                                                                                                                                                              
      #Consider adding options for different thresholds
      data2 <- set.threshold(data.loco.scan,method=input$gwas_threshold,level=0.05)
      

      #Save manhattan plots to list (only for single trait analysis)
      #if length(traits) == 1
      manhattan_plot_list <- list()

      #plot for six models per trait
      #png(file=paste("Manplot_",colnames(data@pheno[i]),".png",sep=""),width=800, height=550) #change directory in your case
      manhattan_plot_list[["all"]] <- manhattan.plot(data2,traits=colnames(data@pheno[i]))+geom_point(size=3)+theme(text = element_text(size = 25),strip.text = element_text(face = "bold"))
      #print(p1)
      #dev.off()  

      #Output the manhattan plots
      output$manhattan_plot <- renderPlot({

          print(manhattan_plot_list[[input$model_select]])

        })
   
  
      #get most significant SNPs per QTL file
      qtl <- get.QTL(data=data2,traits=colnames(data@pheno[i]),bp.window=5e6)
      #knitr::kable(qtl)
      qtl_d <- data.frame(qtl)
      
      output$gwas_stats <-  renderDT({qtl_d}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

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
  


  #######Genomic Diversity analysis

  #Genomic Diversity output files
  diversity_items <- reactiveValues(
    diversity_df = NULL,
    dosage_df = NULL,
    box_plot = NULL,
    het_df = NULL,
    maf_df = NULL


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

    req(input$diversity_file, input$diversity_ploidy, input$zero_value)

    #Input variables (need to add support for VCF file)
    ploidy <- as.numeric(input$diversity_ploidy)
    geno <- input$diversity_file$datapath
    #geno_mat <- read.csv(input$diversity_file$datapath, header = TRUE, check.names = FALSE, row.names = 1)
    #pheno <- read.csv(input$pop_file$datapath, header = TRUE, check.names = FALSE)

    #Status
    updateProgressBar(session = session, id = "pb_diversity", value = 20, title = "Importing VCF")

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
    
    #Status
    updateProgressBar(session = session, id = "pb_diversity", value = 40, title = "Converting to Numeric")

    #Get the genotype values if the updog dosage calls are present
    if ("UD" %in% info_ids) {
      geno_mat <- extract.gt(vcf, element = "UD")
      class(geno_mat) <- "numeric"
      rm(vcf) #Remove vcf
    }else{
      #Extract GT and convert to numeric calls
      geno_mat <- extract.gt(vcf, element = "GT")
      geno_mat <- apply(geno_mat, 2, convert_to_dosage)
      rm(vcf) #Remove VCF
    }

    #} else {
      #Import genotype matrix
     # geno_mat <- read.csv(geno, header = TRUE, row.names = 1, check.names = FALSE)
    #}

    print(class(geno_mat))
    #Convert genotypes to alternate counts if they are the reference allele counts
    #Importantly, the dosage plot is based on the input format NOT the converted genotypes
    is_reference <- (input$zero_value == "Reference Allele Counts")
    convert_genotype_counts <- function(df, ploidy, is_reference = TRUE) {
      if (is_reference) {
        # Convert from reference to alternate alleles
        return(abs(df - ploidy))
      } else {
        # Data already represents alternate alleles
        return(df)
      }
    }


    print("Genotype file successfully imported")
    ######Get MAF plot (Need to remember that the VCF genotypes are likely set as 0 = homozygous reference, where the dosage report is 0 = homozygous alternate) 
    
    #Updated MAF function
    calculateMAF <- function(df, ploidy) {
      if (is.matrix(df)) {
        df <- as.data.frame(df)
      }
      
      #Convert the elements to numeric if they are characters
      df[] <- lapply(df, function(x) if(is.character(x)) as.numeric(as.character(x)) else x)

      allele_frequencies <- apply(df, 1, function(row) {
        non_na_count <- sum(!is.na(row))
        allele_sum <- sum(row, na.rm = TRUE)
        #print(paste("Non-NA count:", non_na_count, "Allele sum:", allele_sum))
        if (non_na_count > 0) {
          allele_sum / (ploidy * non_na_count)
        } else {
          NA
        }
      })

      maf <- ifelse(allele_frequencies <= 0.5, allele_frequencies, 1 - allele_frequencies)
  
      df$AF <- allele_frequencies
      df$MAF <- maf
  
      maf_df <- df[,c("AF", "MAF"), drop = FALSE]

      #Make the row names (SNP ID) the first column
      maf_df <- maf_df %>%
        tibble::rownames_to_column(var = "SNP_ID")
  
      return(maf_df)
    }
   
    # Function to calculate percentages for each genotype in each sample
    calculate_percentages <- function(matrix_data, ploidy) {
      apply(matrix_data, 2, function(col) {
        counts <- table(col)
        prop <- prop.table(counts) * 100
        #max_val <- max(as.numeric(names(counts)))  # Find the maximum value in the column
        prop[as.character(0:ploidy)]  # Adjust the range based on the max value (consider entering the ploidy value explicitly for max_val)
      })
    } 
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

    #Heterozygosity function
    calculate_heterozygosity <- function(genotype_matrix, ploidy = 2) {
      # Determine the heterozygous values based on ploidy
      heterozygous_values <- seq(1, ploidy - 1)
  
      # Create a logical matrix where TRUE represents heterozygous loci
      is_heterozygous <- sapply(genotype_matrix, function(x) x %in% heterozygous_values)
  
      # Count the number of heterozygous loci per sample, ignoring NAs
      heterozygosity_counts <- colSums(is_heterozygous, na.rm = TRUE)
  
      # Calculate the total number of non-NA loci per sample
      total_non_na_loci <- colSums(!is.na(genotype_matrix))
  
      # Compute the proportion of heterozygous loci
      heterozygosity_proportion <- heterozygosity_counts / total_non_na_loci
  
      # Create a dataframe with Sample ID and Observed Heterozygosity
      result_df <- data.frame(
        SampleID = colnames(genotype_matrix),
        ObservedHeterozygosity = heterozygosity_proportion,
        row.names = NULL,
        check.names = FALSE
      )
  
      return(result_df)
    }

    #Convert the genotype calls prior to het,af, and maf calculation
    geno_mat <- data.frame(convert_genotype_counts(df = geno_mat, ploidy = ploidy, is_reference),
                            check.names = FALSE)

    # Calculating heterozygosity for a tetraploid organism
    diversity_items$het_df <- calculate_heterozygosity(geno_mat, ploidy = ploidy)

    print("Heterozygosity success")
    diversity_items$maf_df <- calculateMAF(geno_mat, ploidy = ploidy)

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

  observe({

    req(diversity_items$dosage_df)

    #Plotting
    #pdf("alfalfa_11_GBS_and_Realignment_doubletons_filtered_dosages.pdf")
    box <- ggplot(diversity_items$dosage_df, aes(x=Dosage, y=Percentage, fill=Data)) +
    #geom_point(aes(color = Data), position = position_dodge(width = 0.8), width = 0.2, alpha = 0.5) +  # Add jittered points
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.9) +
    labs(x = "\nDosage", y = "Percentage\n", title = "Genotype Distribution by Sample") +
    #scale_fill_manual(values = c("GBS loci" = "tan3", "DArTag Realignment loci" = "beige")) +
    theme_bw() +
      theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14)
    )
    #dev.off()
    
    diversity_items$box_plot <- box

  })

  output$dosage_plot <- renderPlot({
    
    req(diversity_items$box_plot)
    diversity_items$box_plot

  })

  #Het plot
  output$het_plot <- renderPlot({
    
    req(diversity_items$het_df, input$hist_bins)
    
    #Plot
    #pdf("meng_filtered_alfalfa_11_GBS_DArTag_no_doubletons_sample_heterozygosity.pdf")
    hist(diversity_items$het_df$ObservedHeterozygosity, breaks = as.numeric(input$hist_bins), col = "tan3", border = "black", xlim= c(0,1),
      xlab = "Observed Heterozygosity",
      ylab = "Number of Samples",
      main = "Sample Observed Heterozygosity")

    axis(1, at = seq(0, 1, by = 0.1), labels = TRUE)


    #legend("topright", legend = c("GBS", "DArTag Realignment"), 
    #     fill = c("tan3", "beige"), 
    #     border = c("black", "black"))
    #dev.off()

  })

  #AF Plot
  output$af_plot <- renderPlot({
    
    req(diversity_items$maf_df, input$hist_bins)
    
    #Plot
    hist(diversity_items$maf_df$AF, breaks = as.numeric(input$hist_bins), col = "grey", border = "black", xlab = "Alternate Allele Frequency",
      ylab = "Frequency", main = "Alternate Allele Frequency Distribution")

  })

  #MAF plot
  output$maf_plot <- renderPlot({
    
    req(diversity_items$maf_df, input$hist_bins)
    
    #Plot
    hist(diversity_items$maf_df$MAF, breaks = as.numeric(input$hist_bins), col = "grey", border = "black", xlab = "Minor Allele Frequency (MAF)",
      ylab = "Frequency", main = "Minor Allele Frequency Distribution")

  })

  observe({

    req(diversity_items$het_df)

    output$sample_table <- renderDT({diversity_items$het_df}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))
  
  })

  observe({

    req(diversity_items$maf_df)

    output$snp_table <- renderDT({diversity_items$maf_df}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

    #Plot

  
  })

  ####Genomic Prediction
  #This tab involved 3 observeEvents
    #1) to get the traits listed in the phenotype file
    #2) to input and validate the input files
    #3) to perform the genomic prediction

   #1) Get traits  
  observeEvent(input$trait_file, {
    info_df <- read.csv(input$trait_file$datapath, header = TRUE, check.names = FALSE, nrow = 0)
    trait_var <- colnames(info_df)
    trait_var <- trait_var[2:length(trait_var)]
    #updateSelectInput(session, "pred_trait_info", choices = c("All", trait_var))
    updateVirtualSelect("pred_fixed_info", choices = trait_var, session = session)
    updateVirtualSelect("pred_trait_info", choices = trait_var, session = session)

    #output$passport_table <- renderDT({info_df}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 4)
    #)
  })

  #2) Error check for prediction and save input files
  continue_prediction <- reactiveVal(NULL)
  pred_inputs <- reactiveValues(
    pheno_input = NULL,
    geno_input = NULL
    )
  pred_outputs <- reactiveValues(
    corr_output = NULL,
    box_plot = NULL,
    violin_plot = NULL
    )

  observeEvent(input$prediction_start, {
    #req(pred_inputs$pheno_input, pred_inputs$geno_input)

    #Variables
    ploidy <- as.numeric(input$pred_ploidy)
    geno_path <- input$pred_file$datapath
    pheno <- read.csv(input$trait_file$datapath, header = TRUE, check.names = FALSE)
    row.names(pheno) <- pheno[,1]
    traits <- input$pred_trait_info
    CVs <- as.numeric(input$pred_cv)
    train_perc <- as.numeric(input$pred_train)


  #Make sure at least one trait was input
  if (length(traits) == 0) {

        # If condition is met, show notification toast
        shinyalert(
            title = "Oops",
            text = "No traits were selected",
            size = "xs",
            closeOnEsc = TRUE,
            closeOnClickOutside = FALSE,
            html = TRUE,
            type = "info",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE,
          )
                         
      
        # Stop the observeEvent gracefully
        return()

  }


  #Getting genotype matrix
  geno <- read.csv(geno_path, header = TRUE, row.names = 1, check.names = FALSE)

  #Check that the ploidy entered is correct
  if (ploidy != max(geno, na.rm = TRUE)) {
      # If condition is met, show notification toast
        shinyalert(
            title = "Ploidy Mismatch",
            text = paste0("The maximum value in the genotype file (",max(geno, na.rm = TRUE),") does not equal the ploidy entered"),
            size = "xs",
            closeOnEsc = FALSE,
            closeOnClickOutside = FALSE,
            html = TRUE,
            type = "warning",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            #closeOnConfirm = TRUE,
            #closeOnCancel = TRUE,
            imageUrl = "",
            animation = TRUE
          )
                         
      
        # Stop the observeEvent gracefully
        #return()
      }


  # Function to convert genotype matrix according to ploidy
  convert_genotype <- function(genotype_matrix, ploidy) {
  normalized_matrix <- 2 * (genotype_matrix / ploidy) - 1
  return(normalized_matrix)
  }

  #tranforming genotypes
  geno_adj_init <- convert_genotype(geno, as.numeric(ploidy))

  #Make sure the trait file and genotype file are in the same order
  # Column names for geno (assuming these are the individual IDs)
  colnames_geno <- colnames(geno)
  # Assuming the first column in Pheno contains the matching IDs
  ids_pheno <- pheno[, 1]
  # Find common identifiers
  common_ids <- intersect(colnames_geno, ids_pheno)
  
  #Throw an error if there are less matching samples in the phenotype file than the genotype file     
      if (length(common_ids) == 0) {

        # If condition is met, show notification toast
        shinyalert(
            title = "Oops",
            text = "All samples were missing from the phenotype file",
            size = "xs",
            closeOnEsc = TRUE,
            closeOnClickOutside = FALSE,
            html = TRUE,
            type = "info",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE,
          )
                         
      
        # Stop the observeEvent gracefully
        return()

      } else if (length(common_ids) < length(colnames_geno)) {
      # If condition is met, show notification toast
        shinyalert(
            title = "Data Mismatch",
            text = paste0((length(colnames_geno)-length(common_ids))," samples were removed for not having trait information"),
            size = "xs",
            closeOnEsc = FALSE,
            closeOnClickOutside = FALSE,
            html = TRUE,
            type = "warning",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            #closeOnConfirm = TRUE,
            #closeOnCancel = TRUE,
            imageUrl = "",
            animation = TRUE
          )
                         
      
        # Stop the observeEvent gracefully
        #return()
      } 
      



  #Final check before performing analyses
      shinyalert(
            title = "Ready?",
            text = "Inputs have been checked",
            size = "xs",
            closeOnEsc = FALSE,
            closeOnClickOutside = FALSE,
            html = TRUE,
            type = "info",
            showConfirmButton = TRUE,
            confirmButtonText = "Proceed",
            confirmButtonCol = "#004192",
            showCancelButton = TRUE,
            #closeOnConfirm = TRUE,
            #closeOnCancel = TRUE,
            imageUrl = "",
            animation = TRUE,
            callbackR = function(value) {
              if (isTRUE(value)) {
                # Proceed with adjusted data
                continue_prediction(TRUE)
              } else {
                # Stop or change the process
                continue_prediction(FALSE)
              }
            }
          )
                         

  # Subset and reorder geno and pheno to ensure they only contain and are ordered by common IDs
  geno_adj <- geno_adj_init[, common_ids]  # Assuming that the columns can be directly indexed by IDs
  pheno <- pheno[match(common_ids, ids_pheno), ]

  #Save to reactive values
  pred_inputs$pheno_input <- pheno
  #pred_inputs$geno_adj_input <- geno_adj
  pred_inputs$geno_input <- geno_adj

  })

  #3) Analysis only proceeds once continue_prediction is converted to TRUE
  observe({

  req(continue_prediction(),pred_inputs$pheno_input, pred_inputs$geno_input)

  # Stop analysis if cancel was selected
  if (isFALSE(continue_prediction())) {
      return()
  }

  #Variables
  ploidy <- as.numeric(input$pred_ploidy)
  geno_adj <- pred_inputs$geno_input
  pheno <- pred_inputs$pheno_input
  traits <- input$pred_trait_info
  CVs <- as.numeric(input$pred_cv)
  train_perc <- as.numeric(input$pred_train)
  fixed_traits <- input$pred_fixed_info
  cores <- input$pred_cores


  ##Need to add ability for the use of parallelism for the for cross-validation
  ##Example at this tutorial also: https://www.youtube.com/watch?v=ARWjdQU6ays

  # Function to perform genomic prediction
  ##Make sure this is correct (I think I need to be generating a relationship matrix A.mat() to account for missing data, but I am not sure how that works with this)
  genomic_prediction <- function(geno, Pheno, traits, fixed_effects = NULL, k = 5, percentage = 60, cores = 1) {
  
  # Define variables
  traits <- traits
  cycles <- as.numeric(k)
  total_population <- ncol(geno)
  train_size <- floor(percentage / 100 * total_population)
  fixed_traits <- fixed_effects
  cores <- as.numeric(cores)
  print(cores)
  
  # Establish results matrix
  results <- matrix(nrow = cycles, ncol = length(traits))
  colnames(results) <- traits  # Set the column names to be the traits

  # Initialize a list to store GEBVs for all traits and cycles
  GEBVs <- list()

  #Cross validation number
  pb_value = 10

  #Remove the fixed traits from the Pheno file
  if (length(fixed_traits) == 0) {
    Pheno <- Pheno
  } else {
    #Pheno <- subset(Pheno, select = -fixed_traits)
    Fixed <- subset(Pheno, select = fixed_traits)
  }

  #Make kinship matrix of all individuals?
  #Kin_mat <- A.mat(t(geno), n.core = 1) ##Need to explore whether or not to use a Kinship matrix and if it makes a significant improvement to accuracy
  #If wanting to use Kkinship matrix, will then need to see how to implement it here

  #For now, I am just imputing the missing sites using mean, but EM is more accurate, but slower (can use multiple cores).
  impute = (A.mat(t(geno), max.missing=0.5,impute.method="mean",return.imputed=TRUE))
  geno <- impute$imputed

  # For loop
  for (r in 1:cycles) {

    #Status bar length
    #pb_value = pb_value + floor(70 / as.numeric(cycles))
    pb_value = pb_value + (70 / as.numeric(cycles))

    #Status
    updateProgressBar(session = session, id = "pb_prediction", value = as.numeric(pb_value), status = "info", title = paste0("Performing Cross-Validation:", r, "of", cycles))
    
    train <- as.matrix(sample(1:nrow(geno), train_size))
    test <- setdiff(1:nrow(geno), train)

    #Subset datasets
    if (length(fixed_traits) == 0) {
      Fixed_train = NULL
    } else{
      Fixed_train <- Fixed[train, ]
    }
    Pheno_train <- Pheno[train, ] # Subset the phenotype df to only retain the relevant samples from the training set
    m_train <- geno[train, ]
    Pheno_test <- Pheno[test, ]
    #Fixed_test <- Fixed[test, ] #Where would the Fixed_test be used?
    m_valid <- geno[test, ]

    # Initialize a matrix to store GEBVs for this cycle
    GEBVs_cycle <- matrix(nrow = train_size, ncol = length(traits))
    colnames(GEBVs_cycle) <- traits
    rownames(GEBVs_cycle) <- paste("Cycle", r, "Ind", train, sep="_")

    #Evaluate each trait using the same train and testing samples for each
    for (trait_idx in 1:length(traits)) {
      trait <- Pheno_train[, traits[trait_idx]] # Get the trait of interest
      trait_answer <- mixed.solve(y= trait, Z = m_train, K = NULL, X = Fixed_train, SE = FALSE, return.Hinv = FALSE)
      TRT <- trait_answer$u
      e <- as.matrix(TRT)
      pred_trait_test <- m_valid %*% e
      pred_trait <- pred_trait_test[, 1] + c(trait_answer$beta) # Make sure this still works when using multiple traits
      trait_test <- Pheno_test[, traits[trait_idx]]
      results[r, trait_idx] <- cor(pred_trait, trait_test, use = "complete")

      # Extract GEBVs
      # Check if Fixed_train is not NULL and include beta if it is
      if (!is.null(Fixed_train) && !is.null(trait_answer$beta)) {
        # Calculate GEBVs including fixed effects
        GEBVs_cycle[, trait_idx] <- m_train %*% trait_answer$u + Fixed_train %*% matrix(trait_answer$beta, nrow = length(trait_answer$beta), ncol = 1)
      } else {
        # Calculate GEBVs without fixed effects
        GEBVs_cycle[, trait_idx] <- m_train %*% trait_answer$u
      }

      }

      # Store GEBVs for this cycle
      GEBVs[[r]] <- GEBVs_cycle

    }

    # Combine all GEBVs into a single DataFrame
    GEBVs_df <- do.call(rbind, GEBVs)
  
    results <- as.data.frame(results)
    return(list(GEBVs = GEBVs_df, PredictionAccuracy = results))
  } 

  # Example call to the function
  #This is slow when using 3k markers and 1.2k samples...will need to parallelize if using this script...
  results <- genomic_prediction(geno_adj, pheno, traits = traits, fixed_effects = fixed_traits, k= CVs, percentage= train_perc, cores = cores)

  print(results$PredictionAccuracy)
  #With fixed effects (need to inforporate the ability for fixed effects into the prediction?)
  #results <- genomic_prediction(geno_matrix, phenotype_df, c("height", "weight"), "~ age + sex")

  #Save to reactive value
  pred_outputs$corr_output <- results$PredictionAccuracy

  #TESTING!!!
  write.csv(results$GEBVs, "GEBVs_test.csv")

  #Status
  updateProgressBar(session = session, id = "pb_prediction", value = 90, status = "info", title = "Generating Results")

  ##Figures and Tables

  #Status
  updateProgressBar(session = session, id = "pb_prediction", value = 100, status = "success", title = "Finished")

  #End the event
  continue_prediction(NULL)
  })

#####Analysis Outputs

  observe({
    req(pred_outputs$corr_output)

    df <- pred_outputs$corr_output

    #Probably want to add the ability for the user to select which trait(s) to display here

    #Convert to long format for ggplot
    df_long <- pivot_longer(
    df,
    cols = colnames(df),  # Exclude the Cycle column from transformation
    names_to = "Trait",  # New column for trait names
    values_to = "Correlation"  # New column for correlation values
    )

    #This can be adapted if we start comparing more than one GP model
    #Also consider a violin plot to show each cor value
    #plot <- ggplot(df_long, aes(x = factor(Trait), y = Correlation, fill = "red"), fill = "red") +
    plot <- ggplot(df_long, aes(x = "rrBLUP", y = Correlation, fill = "red"), fill = "red") +
    #geom_boxplot(position = position_dodge(width = 0.8), color = "black", width = 0.7, outlier.size = 0.2) +
    geom_boxplot() +
    facet_wrap(~ Trait, nrow = 1) +  # Facet by trait, allowing different y-scales
    labs(title = "Predicton Accuracy by Trait",
       x = " ",
       y = "Pearson Correlation") +
    #theme_minimal() +                      # Using a minimal theme
    theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), 
        strip.text.x = element_text(face = "bold"))

    plot_violin <- ggplot(df_long, aes(x = "rrBLUP", y = Correlation, fill = "red")) +
    geom_violin(trim = TRUE) +  # Add violin plot
    geom_point(position = position_jitter(width = 0.1), color = "black", size = 1.5) +  # Add jittered points
    facet_wrap(~ Trait, nrow = 1) +  # Facet by trait, allowing different y-scales
    labs(title = "Prediction Accuracy by Trait",
         x = " ",  # x-label is blank because it's not relevant per facet
         y = "Pearson Correlation") +
    theme(legend.position = "none",
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), 
          strip.text.x = element_text(face = "bold"))

  pred_outputs$box_plot <- plot
  pred_outputs$violin_plot <- plot_violin

  })

  #Output the genomic prediction correlation box plots
  output$pred_box_plot <- renderPlot({
    req(pred_outputs$box_plot)
    pred_outputs$box_plot
    })

  #Output the genomic prediction correlation box plots
  output$pred_violin_plot <- renderPlot({
    req(pred_outputs$violin_plot)
    pred_outputs$violin_plot
    })

##Saving analysis outputs

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
      #req(all_plots$pca_2d, all_plots$pca3d, all_plots$scree, input$pca_image_type, input$pca_image_res, input$pca_image_width, input$pca_image_height) #Get the plots
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
        print(all_plots$pca_2d)
      } else if (input$pca_figure == "Scree Plot") {
        print(all_plots$pca_scree)
      #} else if (input$pca_figure == "3D Plot") {
        #print(all_plots$pca3d)  # Assuming you might want a 3D plot as well
      } else {
        plot(x = 1:10, y = 1:10, main = "Fallback Simple Test Plot")  # Fallback simple test plot
      }

      dev.off()
    }

  )

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
      #req(all_plots$pca_2d, all_plots$pca3d, all_plots$scree, input$pca_image_type, input$pca_image_res, input$pca_image_width, input$pca_image_height) #Get the plots
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
          req(diversity_items$box_plot)
          print(diversity_items$box_plot)
      
      } else if (input$div_figure == "AF Histogram") {
          req(diversity_items$maf_df, input$hist_bins)
    
          #Plot
          hist(diversity_items$maf_df$AF, breaks = as.numeric(input$hist_bins), col = "grey", border = "black", xlab = "Alternate Allele Frequency",
            ylab = "Frequency", main = "Alternate Allele Frequency Distribution")

      } else if (input$div_figure == "MAF Histogram") {
        req(diversity_items$maf_df, input$hist_bins)
    
        #Plot
        hist(diversity_items$maf_df$MAF, breaks = as.numeric(input$hist_bins), col = "grey", border = "black", xlab = "Minor Allele Frequency (MAF)",
          ylab = "Frequency", main = "Minor Allele Frequency Distribution")

      } else if (input$div_figure == "OHet Histogram") {
          req(diversity_items$het_df, input$hist_bins)
 
          hist(diversity_items$het_df$ObservedHeterozygosity, breaks = as.numeric(input$hist_bins), col = "tan3", border = "black", xlim= c(0,1),
            xlab = "Observed Heterozygosity",
            ylab = "Number of Samples",
            main = "Sample Observed Heterozygosity")

          axis(1, at = seq(0, 1, by = 0.1), labels = TRUE)

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