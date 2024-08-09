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
#'
mod_gwas_ui <- function(id){
  ns <- NS(id)
  tagList(
    # Add GWAS content here
    fluidRow(
      column(width = 3,
             box(title="Inputs", width = 12, collapsible = TRUE, collapsed = FALSE, status = "info", solidHeader = TRUE,
                 fileInput(ns("gwas_file"), "Choose Genotypes File", accept = c(".csv",".vcf",".gz")),
                 fileInput(ns("phenotype_file"), "Choose Phenotype File", accept = ".csv"),
                 numericInput(ns("gwas_ploidy"), "Species Ploidy", min = 1, value = NULL),
                 selectInput(ns('gwas_threshold'), label='Significance Threshold Method', choices = c("M.eff","Bonferroni","FDR","permute"), selected="M.eff"),
                 selectInput(ns('trait_info'), label = 'Select Trait (eg, Color):', choices = NULL),
                 virtualSelectInput(
                   inputId = ns("fixed_info"),
                   label = "Select Fixed Effects (optional):",
                   choices = NULL,
                   showValueAsTags = TRUE,
                   search = TRUE,
                   multiple = TRUE
                 ),
                 sliderInput(ns("cores"), "Number of CPU Cores", min = 1, max = (availableCores() - 1), value = 1, step = 1),
                 actionButton(ns("gwas_start"), "Run Analysis"),
                 div(style="display:inline-block; float:right",dropdownButton(
                   tags$h3("GWAS Parameters"),
                   "Add description of each filter",
                   circle = FALSE,
                   status = "warning",
                   icon = icon("info"), width = "300px",
                   tooltip = tooltipOptions(title = "Click to see info!")
                 ))
             )
      ),
      column(width = 6,
             box(
               title = "Plots", status = "info", solidHeader = FALSE, width = 12, height = 600,
               bs4Dash::tabsetPanel(
                 tabPanel("BIC Plot", withSpinner(plotOutput(ns("bic_plot"), height = "500px"))),
                 tabPanel("Manhattan Plot", withSpinner(plotOutput(ns("manhattan_plot"), height = "500px"))),
                 tabPanel("QQ Plot", withSpinner(plotOutput(ns("qq_plot"), height = "500px"))),
                 tabPanel("BIC Table", withSpinner(DTOutput(ns("bic_table"))),style = "overflow-y: auto; height: 500px"),
                 tabPanel("QTL - significant markers",
                          withSpinner(DTOutput(ns('gwas_stats'))),style = "overflow-y: auto; height: 500px")
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
                                                                                           "QQ Plot")),
                   selectInput(inputId = ns('gwas_image_type'), label = 'File Type', choices = c("jpeg","tiff","png"), selected = "jpeg"),
                   sliderInput(inputId = ns('gwas_image_res'), label = 'Resolution', value = 300, min = 50, max = 1000, step=50),
                   sliderInput(inputId = ns('gwas_image_width'), label = 'Width', value = 9, min = 1, max = 20, step=0.5),
                   sliderInput(inputId = ns('gwas_image_height'), label = 'Height', value = 5, min = 1, max = 20, step = 0.5),
                   fluidRow(
                     downloadButton(ns("download_gwas_figure"), "Save Image"),
                     downloadButton(ns("download_gwas_file"), "Save Files")),
                   circle = FALSE,
                   status = "danger",
                   icon = icon("floppy-disk"), width = "300px",
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
#' @importFrom matrixcalc is.positive.definite
#' @importFrom Matrix nearPD
#' @importFrom stats BIC as.formula lm logLik median model.matrix na.omit prcomp qbeta quantile runif sd setNames
#' @noRd
mod_gwas_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    #Call some plots to NULL so that the spinners do not show before analysis
    output$bic_plot <- renderDT(NULL)
    output$manhattan_plot <- renderDT(NULL)
    output$qq_plot <- renderDT(NULL)
    output$bic_table <- renderDT(NULL)
    output$gwas_stats <- renderDT(NULL)

    ##GWAS items
    gwas_vars <- reactiveValues(
      gwas_df = NULL,
      manhattan_plots = NULL,
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
      updateSelectInput(session, "trait_info", choices = c("All", trait_var))
      updateVirtualSelect("fixed_info", choices = trait_var, session = session)
    })

    #GWAS analysis (Shufen Chen and Meng Lin pipelines)
    observeEvent(input$gwas_start, {
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
      write.csv(phenotype_file, file = temp_pheno_file, row.names = FALSE)

      #Remove the phenotype_file from memory
      rm(phenotype_file)

      #Status
      updateProgressBar(session = session, id = "pb_gwas", value = 5, title = "Upload Complete: Now Formatting GWASpoly Data")

      #Geno file path
      file_path <- input$gwas_file$datapath

      #Geno.file conversion if needed
      if (grepl("\\.csv$", file_path)) {
        data <- read.GWASpoly(ploidy= ploidy, pheno.file= temp_pheno_file, geno.file=input$gwas_file$datapath,
                              format="numeric", n.traits=length(traits), delim=",") #only need to change files here

      } else if (grepl("\\.vcf$", file_path) || grepl("\\.gz$", file_path)) {
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
        vcf <- read.vcfR(input$gwas_file$datapath)

        #Extract GT
        geno_mat <- extract.gt(vcf, element = "GT")
        geno_mat <- apply(geno_mat, 2, convert_to_dosage)
        class(geno_mat) <- "numeric"
        info <- data.frame(vcf@fix)
        gpoly_df <- cbind(info[,c("ID","CHROM","POS")], geno_mat)
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

      source("R/MyFun_BIC_Meng.R") #change directory in your case

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

        #Display BIC figure
        output$bic_plot <- renderPlot({
          print(p1)
        })
        #dev.off()

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
          model <- c("additive", "1-dom")
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
        manhattan_plot_list[["all"]] <- manhattan.plot(data2,traits=colnames(data@pheno[i]), models = model)+geom_point(size=3)+theme(text = element_text(size = 25),strip.text = element_text(face = "bold"))

        #Output the manhattan plots
        output$manhattan_plot <- renderPlot({

          print(manhattan_plot_list[[input$model_select]])

        })


        #get most significant SNPs per QTL file
        qtl <- get.QTL(data=data2,traits=colnames(data@pheno[i]),bp.window=5e6)
        qtl_d <- data.frame(qtl)

        #Save QTL info
        gwas_vars$gwas_df <- qtl_d

        output$gwas_stats <-  renderDT({qtl_d}, options = list(scrollX = TRUE,autoWidth = FALSE, pageLength = 5))

        #Updating value boxes
        output$qtls_detected <- renderValueBox({
          valueBox(
            value = length(unique(qtl_d$Position)),
            subtitle = "QTLs Detected",
            icon = icon("dna"),
            color = "info"
          )
        })

        #Status
        updateProgressBar(session = session, id = "pb_gwas", value = 80, title = "GWAS Complete: Now Plotting Results")

        #get qqplot
        data_qq <- cbind.data.frame(SNP=data.loco.scan@map$Marker,Chr=data.loco.scan@map$Chrom, Pos=data.loco.scan@map$Position,10^(-data.loco.scan@scores[[colnames(data@pheno[i])]]))

        source("R/CMplot.r") #Obtained the CMplot code from GitHub and made edits to allow inline plotting for shiny app

        #Save qq_plot info
        gwas_vars$qq_plots <- data_qq

        output$qq_plot <- renderPlot({
          CMplot_shiny(data_qq,plot.type="q",col=c(1:8),
                       ylab.pos=2,
                       file.name=colnames(data@pheno[i]),
                       conf.int=FALSE,
                       box=F,multraits=TRUE,file.output=FALSE)
        })

        #plot for each model per trait
        for (j in 1:length(model)) {
          print(j)

          data.loco.scan_2 <- GWASpoly(data=data.loco,models=model[j],
                                       traits=colnames(data@pheno[i]),params=params,n.core= as.numeric(cores))

          data3 <- set.threshold(data.loco.scan_2,method="M.eff",level=0.05)
          manhattan_plot_list[[model[j]]] <- manhattan.plot(data3,traits=colnames(data@pheno[i]))+geom_point(size=3)+theme(text = element_text(size = 25),strip.text = element_text(face = "bold"))
        }

        #Save manhattan plots
        gwas_vars$manhattan_plots <- manhattan_plot_list

      }

      #Status
      updateProgressBar(session = session, id = "pb_gwas", value = 100, status = "success", title = "Finished")

    })

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
          gwas_file <- file.path(temp_dir, paste0("QTL-statistics-", Sys.Date(), ".csv"))
          write.csv(gwas_vars$gwas_df, gwas_file, row.names = FALSE)
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
        } else {
          paste("GWAS-", Sys.Date(), ".tiff", sep="")
        }
      },
      content = function(file) {
        req(input$gwas_figures)

        if (input$gwas_image_type == "jpeg") {
          jpeg(file, width = as.numeric(input$gwas_image_width), height = as.numeric(input$gwas_image_height), res= as.numeric(input$gwas_image_res), units = "in")
        } else if (input$gwas_image_type == "png") {
          png(file, width = as.numeric(input$gwas_image_width), height = as.numeric(input$gwas_image_height), res= as.numeric(input$gwas_image_res), units = "in")
        } else {
          tiff(file, width = as.numeric(input$gwas_image_width), height = as.numeric(input$gwas_image_height), res= as.numeric(input$gwas_image_res), units = "in")
        }

        # Conditional plotting based on input selection
        if (input$gwas_figures == "BIC Plot") {
          req(gwas_vars$BIC_ggplot)
          print(gwas_vars$BIC_ggplot)

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
  })
}

## To be copied in the UI
# mod_gwas_ui("gwas_1")

## To be copied in the server
# mod_gwas_server("gwas_1")
