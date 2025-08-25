#' gwas UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import GWASpoly
#' @importFrom future availableCores
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyWidgets virtualSelectInput
#' @import shinydisconnect
#'
mod_gwas_ui <- function(id){
  ns <- NS(id)
  tagList(
    # Add GWAS content here
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
             box(title="Inputs", width = 12, collapsible = TRUE, collapsed = FALSE, status = "info", solidHeader = TRUE,
                 "* Required",
                 fileInput(ns("gwas_file"), "Choose VCF File*", accept = c(".csv",".vcf",".gz")),
                 fileInput(ns("phenotype_file"), "Choose Trait File*", accept = ".csv"),
                 numericInput(ns("gwas_ploidy"), "Species Ploidy*", min = 1, value = NULL),
                 numericInput(ns("bp_window_before"), "Base pair window (Mb)*", min = 0, value = 2),
                 selectInput(ns('gwas_threshold'), label='Significance Threshold*', choices = c("M.eff","Bonferroni","FDR","permute"), selected="M.eff"),
                 selectInput(ns('trait_info'), label = 'Select Trait*', choices = NULL),
                 virtualSelectInput(
                   inputId = ns("fixed_info"),
                   label = "Select Fixed Effects:",
                   choices = NULL,
                   showValueAsTags = TRUE,
                   search = TRUE,
                   multiple = TRUE
                 ),
                 sliderInput(ns("cores"), "Number of CPU Cores*", min = 1, max = (availableCores() - 1), value = 1, step = 1),
                 actionButton(ns("gwas_start"), "Run Analysis"),
                 div(style="display:inline-block; float:right",dropdownButton(
                   HTML("<b>Input files</b>"),
                   p(downloadButton(ns('download_vcf'),""), "VCF Example File"),
                   p(downloadButton(ns('download_pheno'),""), "Trait Example File"), hr(),
                   p(HTML("<b>Parameters description:</b>"), actionButton(ns("goGWASpar"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>Results description:</b>"), actionButton(ns("goGWASgraph"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>How to cite:</b>"), actionButton(ns("goGWAScite"), icon("arrow-up-right-from-square", verify_fa = FALSE) )), hr(),
                   p(HTML("<b>GWASpoly tutorial:</b>"), actionButton(ns("goGWASpoly"), icon("arrow-up-right-from-square", verify_fa = FALSE), onclick ="window.open('https://jendelman.github.io/GWASpoly/GWASpoly.html', '_blank')" )),hr(),
                   actionButton(ns("gwas_summary"), "Summary"),
                   circle = FALSE,
                   status = "warning",
                   icon = icon("info"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see info!")
                 ))
             )
      ),
      column(width = 6,
             box(
               title = "Plots", status = "info", solidHeader = FALSE, width = 12,
               bs4Dash::tabsetPanel(
                 tabPanel("BIC Plot", withSpinner(plotOutput(ns("bic_plot"), height = "500px"))),
                 tabPanel("Manhattan Plot", withSpinner(plotOutput(ns("manhattan_plot"), height = "500px"))),
                 tabPanel("QQ Plot", withSpinner(plotOutput(ns("qq_plot"), height = "500px"))),
                 tabPanel("BIC Table", withSpinner(DTOutput(ns("bic_table"))),style = "overflow-y: auto; height: 500px"),
                 tabPanel("QTL - significant markers", withSpinner(DTOutput(ns("all_qtl"))),style = "overflow-y: auto; height: 500px"),
                 tabPanel("Filter QTL by LD window",
                          br(),
                          box(
                            title = "LD plot",solidHeader = FALSE, width = 12,
                            plotlyOutput(ns("LD_plot"), height = "500px"), br(),
                            sliderInput(ns("bp_window_after"), label = "Adjust base pair window here to filter QTLs", min = 0,
                                        max = 100, value = 5, step = 1)
                          ),
                          box(
                            title = "Filtered QTL", solidHeader = FALSE, width = 12,
                            withSpinner(DTOutput(ns('gwas_stats')))
                          )
                 ),
                 tabPanel("Multiple QTL model",
                          br(),
                          pickerInput(
                            inputId = ns("sele_models"),
                            label = "Select model",
                            choices = "will be updated",
                            options = list(
                              `actions-box` = TRUE),
                            multiple = FALSE
                          ), hr(),
                          pickerInput(
                            inputId = ns("sele_qtl"),
                            label = "Select QTL",
                            choices = "will be updated",
                            options = list(
                              `actions-box` = TRUE),
                            multiple = TRUE
                          ), hr(),
                          withSpinner(DTOutput(ns('gwas_fitqtl'))),
                          style = "overflow-y: auto; height: 500px")
               )
             )
      ),
      column(width = 3,
             valueBoxOutput(ns("qtls_detected"), width = NULL),
             box(title = "Status", width = 12, collapsible = TRUE, status = "info",
                 progressBar(id = ns("pb_gwas"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
             ),
             box(title = "Plot Controls", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
                 selectInput(ns('model_select'), label = 'Model Selection', choices = NULL),
                 div(style="display:inline-block; float:left",dropdownButton(
                   tags$h3("Save Image"),
                   selectInput(inputId = ns('gwas_figures'), label = 'Figure', choices = c("BIC Plot",
                                                                                           "Manhattan Plot",
                                                                                           "QQ Plot",
                                                                                           "LD Plot")),
                   selectInput(inputId = ns('gwas_image_type'), label = 'File Type', choices = c("jpeg","tiff","png","svg"), selected = "jpeg"),
                   sliderInput(inputId = ns('gwas_image_res'), label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                   sliderInput(inputId = ns('gwas_image_width'), label = 'Width', value = 9, min = 1, max = 20, step=0.5),
                   sliderInput(inputId = ns('gwas_image_height'), label = 'Height', value = 5, min = 1, max = 20, step = 0.5),
                   fluidRow(
                     downloadButton(ns("download_gwas_figure"), "Save Image"),
                     downloadButton(ns("download_gwas_file"), "Save Files"),
                     downloadButton(ns("download_viewpoly"), "VIEWpoly input")),
                   circle = FALSE,
                   status = "danger",
                   icon = icon("floppy-disk"), width = "300px", label = "Save",
                   tooltip = tooltipOptions(title = "Click to see inputs!")
                 ))
             )
      )
    )
  )
}

#' gwas Server Functions
#'
#' @importFrom DT renderDT
#' @importFrom vcfR read.vcfR
#' @importFrom Matrix nearPD
#' @importFrom stats BIC as.formula lm logLik median model.matrix na.omit prcomp qbeta quantile runif sd setNames
#' @importFrom bs4Dash updatebs4TabItems updateBox
#' @importFrom shiny updateTabsetPanel
#' @importFrom plotly ggplotly
#' @import dplyr
#' @noRd
mod_gwas_server <- function(input, output, session, parent_session){

  ns <- session$ns

  # Help links
  observeEvent(input$goGWASpar, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "GWAS_tabset",
                      selected = "GWAS_par")
    # expand specific box
    updateBox(id = "GWAS_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goGWASgraph, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "GWAS_tabset",
                      selected = "GWAS_results")
    # expand specific box
    updateBox(id = "GWAS_box", action = "toggle", session = parent_session)
  })

  observeEvent(input$goGWAScite, {
    # change to help tab
    updatebs4TabItems(session = parent_session, inputId = "MainMenu",
                      selected = "help")

    # select specific tab
    updateTabsetPanel(session = parent_session, inputId = "GWAS_tabset",
                      selected = "GWAS_cite")
    # expand specific box
    updateBox(id = "GWAS_box", action = "toggle", session = parent_session)
  })
  
  #Default choices
  trait_options <- reactiveValues(
    missing_data = "NA",
    custom_missing = NULL,
    sample_column = NULL,
    file_type = NULL
  )
  
  #UI popup window for input
  observeEvent(input$phenotype_file, {
    req(input$phenotype_file)
    #Get the column names of the csv file
    info_df <- read.csv(input$phenotype_file$datapath, header = TRUE, check.names = FALSE, nrows=2)
    info_df[,1] <- as.character(info_df[,1]) #Makes sure that the sample names are characters instead of numeric
    
    # Read first 5 rows for preview
    preview_data <- tryCatch({
      head(read.csv(input$phenotype_file$datapath, nrows = 5, na.strings=trait_options$missing_data),5)
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
  

  #Call some plots to NULL so that the spinners do not show before analysis
  output$bic_plot <- renderDT(NULL)
  output$manhattan_plot <- renderDT(NULL)
  output$qq_plot <- renderDT(NULL)
  output$bic_table <- renderDT(NULL)
  output$gwas_stats <- renderDT(NULL)

  ##GWAS items
  gwas_data <- reactiveValues(
    data2 = NULL,
    phenos = NULL
  )

  gwas_vars <- reactiveValues(
    gwas_df = NULL,
    gwas_df_filt = NULL,
    fit_qtl = NULL,
    manhattan_plots = NULL,
    LD_plot = NULL,
    bp_window = NULL,
    qq_plots = NULL,
    bic_df = NULL,
    BIC_ggplot = NULL
  )

  output$qtls_detected <- renderValueBox({
    valueBox(
      value = 0,
      subtitle = "QTLs Detected",
      icon = icon("dna"),
      color = "info"
    )
  })

  observeEvent(input$phenotype_file, {
    info_df <- read.csv(input$phenotype_file$datapath, header = TRUE, check.names = FALSE, nrow = 0)
    trait_var <- colnames(info_df)
    trait_var <- trait_var[2:length(trait_var)]
    updateSelectInput(session, "trait_info", choices = c(trait_var))
    updateVirtualSelect("fixed_info", choices = trait_var, session = session)
  })

  #GWAS analysis (Shufen Chen and Meng Lin pipelines)
  observeEvent(input$gwas_start, {
    toggleClass(id = "gwas_ploidy", class = "borderred", condition = (is.na(input$gwas_ploidy) | is.null(input$gwas_ploidy)))
    toggleClass(id = "trait_info", class = "borderred", condition = (all(is.na(input$trait_info)) | all(is.null(input$trait_info))))

    if (is.null(input$phenotype_file$datapath) | is.null(input$gwas_file$datapath)) {
      shinyalert(
        title = "Missing input!",
        text = "Upload VCF and phenotype files",
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
    req(input$phenotype_file$datapath, input$gwas_file$datapath, input$gwas_ploidy, input$trait_info)

    cores <- input$cores
    #Status
    updateProgressBar(session = session, id = "pb_gwas", value = 0, title = "Uploading Data")

    #Make subset phenotype file (need to develop alternative method that does not save a new phenotype file each time.)
    #I think I can subset the read.GWAS file pheno and fixed categories (data@pheno[,c("trait")]) and data@fixed = phenotype_file[,c("List of fixed traits")]
    # assign the missing data and sample ID column based on user selections
    if (trait_options$missing_data == "(blank)") {
      phenotype_file <- read.csv(input$phenotype_file$datapath, header = TRUE, check.names = FALSE, na.strings="")
    } else if (trait_options$missing_data == "Custom") {
      phenotype_file <- read.csv(input$phenotype_file$datapath, header = TRUE, check.names = FALSE, na.strings = trait_options$custom_missing)
    } else {
      phenotype_file <- read.csv(input$phenotype_file$datapath, header = TRUE, check.names = FALSE, na.strings = trait_options$missing_data)
    }
    
    # Make the sample ID column the first column in the dataframe
    sample_col_name <- input$sample_column
    phenotype_file <- phenotype_file[, c(sample_col_name, setdiff(names(phenotype_file), sample_col_name))]

    # Remove empty lines
    rm.empty <- which(apply(phenotype_file, 1, function(x) all(is.na(x) | x == "")))
    if(length(rm.empty) > 0){
      warning(paste("Removing", length(rm.empty),"empty lines"))
      phenotype_file <- phenotype_file[-rm.empty,]
    }

    ids <- colnames(phenotype_file)[1]
    traits <- input$trait_info
    fixed <- input$fixed_info
    included_var <- c(ids, traits, fixed)
    ploidy <- as.numeric(input$gwas_ploidy)

    # Check if traits are numerical
    n_traits <- as.matrix(phenotype_file[,traits])
    n_traits <- apply(n_traits, 2, function(x) all(is.na(as.numeric(x))))

    if(any(n_traits)){
      shinyalert(
        title = "Input not supported",
        text = paste("All selected traits must be numerical. Categorial traits found:",if(length(n_traits) > 1) names(which(n_traits)) else input$trait_info),
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

    validate(
      need(!any(n_traits), "The selected traits must be numerical.")
    )

    phenotype_file <- phenotype_file[,included_var]

    # Create a temporary file for the selected phenotype data
    temp_pheno_file <- tempfile(fileext = ".csv")

    #Save new phenotype file with selected traits and fixed effects
    write.csv(phenotype_file, file = temp_pheno_file, row.names = FALSE)

    #Status
    updateProgressBar(session = session, id = "pb_gwas", value = 5, title = "Upload Complete: Now Formatting GWASpoly Data")

    #Geno file path
    file_path <- input$gwas_file$datapath

    #Geno.file conversion if needed
    if (grepl("\\.csv$", file_path)) {
      #TODO: Add check for matches of sample names in genotype and phenotype data

      data <- read.GWASpoly(ploidy= ploidy, pheno.file= temp_pheno_file, geno.file=input$gwas_file$datapath,
                            format="numeric", n.traits=length(traits), delim=",") #only need to change files here

    } else if (grepl("\\.vcf$", file_path) || grepl("\\.gz$", file_path)) {
      # Create a temporary file for the selected phenotype data
      temp_geno_file <- tempfile(fileext = ".csv")

      #Convert VCF file if submitted
      #### VCF sanity check
      checks <- vcf_sanity_check(input$gwas_file$datapath, max_markers = 10000)
      
      error_if_false <- c(
        "VCF_header", "VCF_columns", "unique_FORMAT", "GT",
        "samples", "chrom_info", "pos_info", "VCF_compressed"
      )
      
      error_if_true <- c(
        "multiallelics", "phased_GT",  "mixed_ploidies",
        "duplicated_samples", "duplicated_markers"
      )
      
      warning_if_false <- c("ref_alt","max_markers")
      
      checks_result <- vcf_sanity_messages(checks, 
                                           error_if_false, 
                                           error_if_true, 
                                           warning_if_false = warning_if_false, 
                                           warning_if_true = NULL,
                                           input_ploidy = ploidy)
      
      if(checks_result) return() # Stop the analysis if checks fail
      #########
      
      vcf <- read.vcfR(input$gwas_file$datapath, verbose = FALSE)

      #Extract GT
      geno_mat <- extract.gt(vcf, element = "GT")
      geno_mat <- apply(geno_mat, 2, convert_to_dosage)
      class(geno_mat) <- "numeric"
      info <- data.frame(vcf@fix)
      gpoly_df <- cbind(info[,c("ID","CHROM","POS")], geno_mat)

      #Fix ID column if empty
      if (is.na(gpoly_df[1,1])) {
        gpoly_df$ID <- row.names(gpoly_df)
      }

      if(!any(colnames(gpoly_df) %in% phenotype_file[,1])) {
        shinyalert(
          title = "Samples ID do not match",
          text = paste("Check if passport/phenotype files have same sample ID as the VCF/genotype file."),
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
      validate(
        need(any(colnames(gpoly_df) %in% phenotype_file[,1]), "The selected traits must be numerical.")
      )

      write.csv(gpoly_df, file = temp_geno_file, row.names = FALSE)

      data <- read.GWASpoly(ploidy= ploidy, pheno.file= temp_pheno_file, geno.file=temp_geno_file,
                            format="numeric", n.traits=length(traits), delim=",")
      rm(geno_mat)
      rm(gpoly_df)
      rm(vcf)

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

    gwas_vars$LD_plot <- LD.plot(data)
    lim.d <- max(gwas_vars$LD_plot$data$d)

    if(input$bp_window_before > lim.d)
      shinyalert(
        title = "Adjust base pair window",
        text = paste0("Base pair window larger than maximum distance (",lim.d,"). Reduce window size."),
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

    validate(
      need(input$bp_window_before <= lim.d, paste0("Base pair window larger than maximum distance (",lim.d,"). Reduce window size."))
    )

    gwas_vars$bp_window <- input$bp_window_before
    updateSliderInput(session = session, inputId = "bp_window_after", min = 0, max = round(lim.d,2), value = gwas_vars$bp_window, step = round(lim.d/150,4))

    data.loco <- set.K(data,LOCO=F,n.core= as.numeric(cores))

    #Delete temp pheno file
    unlink(temp_pheno_file)

    ####Pheno, kinship, PCs from results of GWASpoly
    GE<- data@pheno
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

    #### calculate BIC
    #Status
    updateProgressBar(session = session, id = "pb_gwas", value = 20, title = "Formatting Complete: Now Calculating BIC")

    PC<-as.matrix(PCs)
    K=as.matrix(Kin)

    kin.adj<-posdefmat(K)
    kin.test<-as.matrix(kin.adj)

    phenos <- vector()
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
      output$bic_table <- renderDT({plotBICs_kinship}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5)
      )

      p1<-ggplot(plotBICs_kinship, aes(x=n.PC, y=BIC,group=RelationshipMatrix)) +
        geom_line(color="grey")+
        geom_point(shape=21, color="black", fill="#d95f0e", size=3)+
        theme(text=element_text(size=15),axis.text.x = element_text(angle =0),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        labs(x = "Number of PCs",y="BIC")

      #Save BIC plot
      gwas_vars$BIC_ggplot <- p1

      #Save BIC plot info
      gwas_vars$bic_df <- plotBICs_kinship

      #Status
      updateProgressBar(session = session, id = "pb_gwas", value = 40, title = "BIC Complete: Now Performing GWAS")

      ##GWAS based on model selection
      N <- nrow(data@pheno) #Population size
      #Select models depending on ploidy
      if (ploidy > 2) {
        model <- c("additive","1-dom","2-dom","general","diplo-general","diplo-additive")
        updateSelectInput(session, "model_select", choices = c("all", model))
      }else{
        model <- c("additive", "1-dom", "general")
        updateSelectInput(session, "model_select", choices = c("all", model))
      }

      BIC_min <- plotBICs_kinship[which.min(plotBICs_kinship$BIC),]
      if(BIC_min$n.PC == 0){params <- set.params(geno.freq = 1 - 5/N)}else{params <- set.params(geno.freq = 1 - 5/N,n.PC = as.numeric(levels(BIC_min$n.PC))[BIC_min$n.PC])}
      data.loco.scan <- GWASpoly(data=data.loco,models=model,traits=colnames(data@pheno[i]),params=params,n.core=as.numeric(cores))
      #Consider adding options for different thresholds
      data2 <- set.threshold(data.loco.scan,method=input$gwas_threshold,level=0.05)

      #Save manhattan plots to list (only for single trait analysis)
      #if length(traits) == 1
      manhattan_plot_list <- list()

      #plot for six models per trait
      manhattan_plot_list[["all"]] <- manhattan.plot(data2,traits=colnames(data@pheno[i]), models = model)+
        geom_point(size=3)+
        theme(text = element_text(size = 25),strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

      #Status
      updateProgressBar(session = session, id = "pb_gwas", value = 80, title = "GWAS Complete: Now Plotting Results")

      #get qqplot
      data_qq <- cbind.data.frame(SNP=data.loco.scan@map$Marker,Chr=data.loco.scan@map$Chrom, Pos=data.loco.scan@map$Position,10^(-data.loco.scan@scores[[colnames(data@pheno[i])]]))

      #Save qq_plot info
      gwas_vars$qq_plots <- data_qq

      #plot for each model per trait
      for (j in 1:length(model)) {
        data.loco.scan_2 <- GWASpoly(data=data.loco,models=model[j],
                                     traits=colnames(data@pheno[i]),params=params,n.core= as.numeric(cores))

        data3 <- set.threshold(data.loco.scan_2,method=input$gwas_threshold,level=0.05)
        manhattan_plot_list[[model[j]]] <- manhattan.plot(data3,traits=colnames(data@pheno[i]))+geom_point(size=3)+
          theme(text = element_text(size = 25),
                strip.text = element_text(face = "bold"),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      }

      #Save manhattan plots
      gwas_vars$manhattan_plots <- manhattan_plot_list
      phenos[i] <- colnames(data@pheno[i])
    }

    gwas_data$data2 <- data2
    gwas_data$phenos <- phenos[-which(is.na(phenos))]

    qtl <- get.QTL(data=gwas_data$data2,traits=gwas_data$phenos,bp.window=0)
    gwas_vars$gwas_df <- data.frame(qtl)

    #Status
    updateProgressBar(session = session, id = "pb_gwas", value = 100, status = "success", title = "Finished")
  })

  #Checking if any QTLs were detected and returning a user notice if not
  observe({
    req(gwas_vars$gwas_df)
    if(dim(gwas_vars$gwas_df)[1] == 0) {
      shinyalert(
        title = "No QTL Detected",
        text = "No QTL detected for this trait.",
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

    #Gracefully abort
    return()
  })

  #Updating value boxes
  output$qtls_detected <- renderValueBox({
    valueBox(
      value = length(unique(gwas_vars$gwas_df_filt$Position)),
      subtitle = "QTLs Detected",
      icon = icon("dna"),
      color = "info"
    )
  })

  # Tables
  output$all_qtl <-  renderDT({
    #get most significant SNPs per QTL file
    validate(
      need(dim(gwas_vars$gwas_df)[1] > 0, "No QTL detected.")
    )
    gwas_vars$gwas_df
  }, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))


  output$gwas_stats <-  renderDT({
    #get most significant SNPs per QTL file
    lim.d <- max(gwas_vars$LD_plot$data$d)

    validate(
      need(gwas_vars$bp_window <= lim.d, paste0("Base pair window larger than maximum distance (",lim.d,"). Reduce window size."))
    )
    if(is.null(input$bp_window_after)) {
      line <- gwas_vars$bp_window
    } else line <- input$bp_window_after

    qtl <- get.QTL(data=gwas_data$data2,traits=gwas_data$phenos,bp.window=line*1000000)
    gwas_vars$gwas_df_filt <- data.frame(qtl)

    validate(
      need(dim(gwas_vars$gwas_df_filt)[1] > 0, "No QTL detected.")
    )
    gwas_vars$gwas_df_filt
  }, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))


  observe({
    req(gwas_vars$gwas_df_filt, nrow(gwas_vars$gwas_df_filt) > 0)
    updatePickerInput(session = session, inputId = "sele_models", choices = unique(gwas_vars$gwas_df_filt$Model), selected = unique(gwas_vars$gwas_df_filt$Model)[1])
  })

  observe({
    req(gwas_vars$gwas_df_filt, nrow(gwas_vars$gwas_df_filt) > 0)

    df <- gwas_vars$gwas_df_filt %>% filter(Model %in% input$sele_models)
    updatePickerInput(session = session, inputId = "sele_qtl", choices = unique(paste0(df$Marker, "_", df$Model)),
                      selected = unique(paste0(df$Marker, "_", df$Model)))
  })

  output$gwas_fitqtl <-  renderDT({
    validate(
      need(dim(gwas_vars$gwas_df_filt)[1] > 0, "No QTL detected.")
    )

    df <- gwas_vars$gwas_df_filt[which(paste0(gwas_vars$gwas_df_filt$Marker, "_", gwas_vars$gwas_df_filt$Model) %in% input$sele_qtl),]

    rm.qtl <- which(df$Model %in% c("diplo-general", "diplo-additive"))
    if(length(rm.qtl) > 0){
      shinyalert(
        title = "Oops",
        text = "QTL detected by the models diplo-general and diplo-additive are not supported in the fit.QTL current version",
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
    }

    if(length(df$Model) >0){
      rm.qtl <- which(df$Model %in% c("diplo-general", "diplo-additive"))
      if(length(rm.qtl) > 0){
        warning("QTL detected by the models diplo-general and diplo-additive are not supported in the fit.QTL current version")
        qtl <- df[-rm.qtl,]
      } else qtl <- df

      validate(
        need(dim(qtl)[1] > 0, "No QTL evaluated")
      )

      fit.ans_temp <- fit.QTL(data=gwas_data$data2,
                              trait=input$trait_info,
                              qtl=qtl[,c("Marker","Model")])
      gwas_vars$fit_qtl <- fit.ans_temp
    } else gwas_vars$fit_qtl <- NULL

    gwas_vars$fit_qtl
  }, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

  # Plots
  #Output the manhattan plots
  output$manhattan_plot <- renderPlot({
    validate(
      need(!is.null(gwas_vars$manhattan_plots), "Upload the input files, set the parameters and click 'run analysis' to access results in this session.")
    )
    print(gwas_vars$manhattan_plots[[input$model_select]])

  })

  output$qq_plot <- renderPlot({
    validate(
      need(!is.null(gwas_vars$qq_plots), "Upload the input files, set the parameters and click 'run analysis' to access results in this session.")
    )
    CMplot_shiny(gwas_vars$qq_plots,plot.type="q",col=c(1:8),
                 ylab.pos=2,
                 file.name=input$trait_info,
                 conf.int=FALSE,
                 box=F,
                 multraits=TRUE,
                 file.output=FALSE)
  })

  #Display BIC figure
  output$bic_plot <- renderPlot({
    validate(
      need(!is.null(gwas_vars$BIC_ggplot), "Upload the input files, set the parameters and click 'run analysis' to access results in this session.")
    )

    print(gwas_vars$BIC_ggplot)
  })

  output$LD_plot <- renderPlotly({

    validate(
      need(!is.null(gwas_vars$LD_plot), "Upload the input files, set the parameters and click 'run analysis' to access results in this session.")
    )

    lim.d <- max(gwas_vars$LD_plot$data$d)

    validate(
      need(gwas_vars$bp_window <= lim.d, paste0("Base pair window larger than maximum distance (",lim.d,"). Reduce window size."))
    )
    if(is.null(input$bp_window_after)) {
      line <- gwas_vars$bp_window
    } else line <- input$bp_window_after

    p <- gwas_vars$LD_plot + geom_vline(aes(xintercept=line, color = "bp window"),linetype="dashed") +
      theme(legend.title=element_blank(), legend.position=c(1,1),legend.justification = c(1,1), text = element_text(size = 15)) +
      labs(y = "R-squared")

    updateNumericInput(session = session, inputId = "bp_window_before",value = line)

    ggplotly(p) %>%
      layout(
        legend = list(
          title = list(text = ''),  # Explicitly remove legend title in plotly
          x = 1,  # X position (right)
          y = 1,  # Y position (top)
          xanchor = 'right',  # Anchor the legend at the right
          yanchor = 'top'  # Anchor the legend at the top
        )
      )
  })

  
  
  output$download_viewpoly <- downloadHandler(
    filename = function() {
      paste0("BIGapp_Viewpoly_GWASpoly.RData")
    },
    content = function(file) {
      temp <- gwas_data$data2
      save(temp, file = file)
    }
  )
  
  #Download files for GWAS
  output$download_gwas_file <- downloadHandler(
    filename = function() {
      paste0("GWAS-results-", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Temporary files list
      temp_dir <- tempdir()
      temp_files <- c()

      if (!is.null(gwas_vars$gwas_df)) {
        # Create a temporary file for assignments
        gwas_file <- file.path(temp_dir, paste0("QTL-Significant_Markers-statistics-", Sys.Date(), ".csv"))
        write.csv(gwas_vars$gwas_df, gwas_file, row.names = FALSE)
        temp_files <- c(temp_files, gwas_file)
      }

      if (!is.null(gwas_vars$gwas_df_filt)) {
        # Create a temporary file for assignments
        gwas_file <- file.path(temp_dir, paste0("QTL-LD-filtered-statistics-", Sys.Date(), ".csv"))
        write.csv(gwas_vars$gwas_df_filt, gwas_file, row.names = FALSE)
        temp_files <- c(temp_files, gwas_file)
      }

      if (!is.null(gwas_vars$fit_qtl)) {
        # Create a temporary file for assignments
        gwas_file <- file.path(temp_dir, paste0("Multiple-QTL-model-statistics-", Sys.Date(), ".csv"))
        write.csv(gwas_vars$fit_qtl, gwas_file, row.names = FALSE)
        temp_files <- c(temp_files, gwas_file)
      }

      if (!is.null(gwas_vars$bic_df)) {
        # Create a temporary file for BIC data frame
        bic_file <- file.path(temp_dir, paste0("GWAS-BIC-statistics-", Sys.Date(), ".csv"))
        write.csv(gwas_vars$bic_df, bic_file, row.names = FALSE)
        temp_files <- c(temp_files, bic_file)
      }

      # Zip files only if there's something to zip
      if (length(temp_files) > 0) {
        zip(file, files = temp_files, extras = "-j") # Using -j to junk paths
      }

      # Optionally clean up
      file.remove(temp_files)
    }
  )


  #Download Figures for GWAS Tab (Need to convert figures to ggplot)
  output$download_gwas_figure <- downloadHandler(

    filename = function() {
      if (input$gwas_image_type == "jpeg") {
        paste("GWAS-", Sys.Date(), ".jpg", sep="")
      } else if (input$gwas_image_type == "png") {
        paste("GWAS-", Sys.Date(), ".png", sep="")
      } else if (input$gwas_image_type == "tiff") {
        paste("GWAS-", Sys.Date(), ".tiff", sep="")
      } else {
        paste("GWAS-", Sys.Date(), ".svg", sep="")
      }
    },
    content = function(file) {
      req(input$gwas_figures)

      if (input$gwas_image_type == "jpeg") {
        jpeg(file, width = as.numeric(input$gwas_image_width), height = as.numeric(input$gwas_image_height), res= as.numeric(input$gwas_image_res), units = "in")
      } else if (input$gwas_image_type == "png") {
        png(file, width = as.numeric(input$gwas_image_width), height = as.numeric(input$gwas_image_height), res= as.numeric(input$gwas_image_res), units = "in")
      } else if (input$gwas_image_type == "tiff") {
        tiff(file, width = as.numeric(input$gwas_image_width), height = as.numeric(input$gwas_image_height), res= as.numeric(input$gwas_image_res), units = "in")
      } else {
        svg(file, width = as.numeric(input$gwas_image_width), height = as.numeric(input$gwas_image_height))
      }

      # Conditional plotting based on input selection
      if (input$gwas_figures == "BIC Plot") {
        req(gwas_vars$BIC_ggplot)
        print(gwas_vars$BIC_ggplot)

      } else if (input$gwas_figures == "LD Plot") {
        req(gwas_vars$LD_plot)
        #Plot
        print(gwas_vars$LD_plot)

      } else if (input$gwas_figures == "Manhattan Plot") {
        req(gwas_vars$manhattan_plots, input$model_select)
        #Plot
        print(gwas_vars$manhattan_plots[[input$model_select]])

      } else if (input$gwas_figures == "QQ Plot") {
        req(gwas_vars$qq_plots)

        source("R/CMplot.r") #Obtained the CMplot code from GitHub and made edits to allow inline plotting for shiny app

        #Plot
        CMplot_shiny(gwas_vars$qq_plots,plot.type="q",col=c(1:8),
                     ylab.pos=2,
                     conf.int=FALSE,
                     box=F,multraits=TRUE,file.output=FALSE)

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

  ##Summary Info
  gwas_summary_info <- function() {
    #Handle possible NULL values for inputs
    dosage_file_name <- if (!is.null(input$gwas_file$name)) input$gwas_file$name else "No file selected"
    passport_file_name <- if (!is.null(input$phenotype_file$name)) input$phenotype_file$name else "No file selected"
    selected_ploidy <- if (!is.null(input$gwas_ploidy)) as.character(input$gwas_ploidy) else "Not selected"

    #Print the summary information
    cat(
      "BIGapp GWAS Summary\n",
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
      paste("Significance Threshold Method:", input$gwas_threshold), "\n",
      paste("Selected Trait:", input$trait_info), "\n",
      paste("Selected Fixed Effects:", input$fixed_info), "\n",
      "\n",
      "### R Packages Used ###\n",
      "\n",
      paste("BIGapp:", packageVersion("BIGapp")), "\n",
      paste("AGHmatrix:", packageVersion("AGHmatrix")), "\n",
      paste("ggplot2:", packageVersion("ggplot2")), "\n",
      paste("GWASpoly:", packageVersion("GWASpoly")), "\n",
      paste("vcfR:", packageVersion("vcfR")), "\n",
      paste("Matrix:", packageVersion("Matrix")), "\n",
      sep = ""
    )
  }

  # Popup for analysis summary
  observeEvent(input$gwas_summary, {
    showModal(modalDialog(
      title = "Summary Information",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close"),
        downloadButton("download_gwas_info", "Download")
      ),
      pre(
        paste(capture.output(gwas_summary_info()), collapse = "\n")
      )
    ))
  })


  # Download Summary Info
  output$download_gwas_info <- downloadHandler(
    filename = function() {
      paste("gwas_summary_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      # Write the summary info to a file
      writeLines(paste(capture.output(gwas_summary_info()), collapse = "\n"), file)
    }
  )
}

## To be copied in the UI
# mod_gwas_ui("gwas_1")

## To be copied in the server
# mod_gwas_server("gwas_1")
