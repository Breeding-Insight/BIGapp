#' DosageCall UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shinyjs enable disable useShinyjs
#' @importFrom shiny NS tagList
#' @importFrom future availableCores
#' @importFrom bs4Dash renderValueBox
#'
#'
mod_DosageCall_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        box(
          title = "Inputs", status = "info", solidHeader = TRUE, collapsible = FALSE, collapsed = FALSE,
          fileInput(ns("madc_file"), "Choose MADC or VCF File", accept = c(".csv",".vcf",".gz")),
          textInput(ns("output_name"), "Output File Name"),
          numericInput(ns("ploidy"), "Species Ploidy", min = 1, value = NULL),
          selectInput(ns("updog_model"), "Updog Model", choices = c("norm","hw","bb","s1","s1pp","f1","f1pp","flex","uniform"), selected = "norm"),
          numericInput(ns("cores"), "Number of CPU Cores", min = 1, max = (future::availableCores() - 1), value = 1),
          actionButton(ns("run_analysis"), "Run Analysis"),
          useShinyjs(),
          downloadButton(ns('download_updog_vcf'), "Download VCF File", class = "butt"),

          div(style="display:inline-block; float:right",dropdownButton(

            tags$h3("Updog Population Models"),
            "Model: What form should the prior (genotype distribution) take?\n
                    The following information is from the Updog manual:\n
                    Possible values of the genotype distribution (values of model) are: \n
                    `norm` A distribution whose genotype frequencies are proportional to the density value of a normal
                    with some mean and some standard deviation. Unlike the `bb` and `hw` options, this will
                    allow for distributions both more and less dispersed than a binomial. This seems to be the
                    most robust to violations in modeling assumptions, and so is the default. This prior class was
                    developed in Gerard and Ferrao (2020).
                    `hw` A binomial distribution that results from assuming that the population is in Hardy-Weinberg
                    equilibrium (HWE). This actually does pretty well even when there are minor to moderate
                    deviations from HWE. Though it does not perform as well as the `norm` option when there
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
        valueBoxOutput(ns("MADCsnps"))
      ),

      fluidRow(
        box(title = "Status", width = 3, collapsible = TRUE, status = "info",
            progressBar(id = ns("pb_madc"), value = 0, status = "info", display_pct = TRUE, striped = TRUE, title = " ")
        )
      )
    )
  )
}

#' DosageCall Server Functions
#'
#' @import vcfR
#' @import updog
#' @importFrom BIGr updog2vcf
#' @importFrom shinyjs enable disable
#'
#' @noRd
mod_DosageCall_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    snp_number <- reactiveVal(0)

    #SNP counts value box
    output$MADCsnps <- renderValueBox({
      valueBox(snp_number(), "Markers in uploaded file", icon = icon("dna"), color = "info")
    })

    disable("download_updog_vcf")

    ##This is for performing Updog Dosage Calling
    updog_out <- eventReactive(input$run_analysis,{

      # Missing input with red border and alerts
      toggleClass(id = "ploidy", class = "borderred", condition = (is.na(input$ploidy) | is.null(input$ploidy)))
      toggleClass(id = "output_name", class = "borderred", condition = (is.na(input$output_name) | is.null(input$output_name) | input$output_name == ""))

      if (is.null(input$madc_file$datapath)) {
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
      req(input$madc_file$datapath, input$output_name, input$ploidy)

      # Get inputs
      madc_file <- input$madc_file$datapath
      output_name <- input$output_name
      ploidy <- input$ploidy
      cores <- input$cores
      model_select <- input$updog_model

      # Status
      updateProgressBar(session = session, id = "pb_madc", value = 0, title = "Formatting Input Files")
      #Import genotype info if genotype matrix format
      if (grepl("\\.csv$", madc_file)) {
        # Call the get_counts function with the specified MADC file path and output file path
        #Status
        result_df <- get_counts(madc_file, output_name)

        #Call the get_matrices function
        matrices <- get_matrices(result_df)

        #Number of SNPs
        snp_number <- (nrow(result_df) / 2)

        #SNP counts value box
        output$MADCsnps <- renderValueBox({
          valueBox(snp_number, "Markers in MADC File", icon = icon("dna"), color = "info")
        })

      } else {

        #Initialize matrices list
        matrices <- list()

        #Import genotype information if in VCF format
        vcf <- read.vcfR(madc_file, verbose = FALSE)

        #Get items in FORMAT column
        info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT
        chrom <- vcf@fix[,1]
        pos <- vcf@fix[,2]

        info_ids <- extract_info_ids(info[1])

        if (("DP" %in% info_ids) && (("RA" %in% info_ids) | ("AD" %in% info_ids))) {
          #Extract DP and RA and convert to matrices
          matrices$size_matrix <- extract.gt(vcf, element = "DP")
          if("RA" %in% info_ids){
            matrices$ref_matrix <- extract.gt(vcf, element = "RA")
          } else {
            ad_matrix <- extract.gt(vcf, element = "AD")
            matrices$ref_matrix <- matrix(sapply(strsplit(ad_matrix, ","), "[[", 1), nrow = nrow(matrices$size_matrix))
            colnames(matrices$ref_matrix) <- colnames(matrices$size_matrix)
            rownames(matrices$ref_matrix) <- rownames(matrices$size_matrix)
          }

          class(matrices$size_matrix) <- "numeric"
          class(matrices$ref_matrix) <- "numeric"
          rownames(matrices$size_matrix) <- rownames(matrices$ref_matrix) <- paste0(chrom, "_", pos)

          rm(vcf) #Remove VCF

          snp_number <- (nrow(matrices$size_matrix))

          #SNP counts value box
          output$MADCsnps <- renderValueBox({
            valueBox(snp_number, "Markers in VCF File", icon = icon("dna"), color = "info")
          })

        }else{
          ##Add user warning about read depth and allele read depth not found
          stop(safeError("Error: DP and RA/AD FORMAT flags not found in VCF file"))
        }
      }

      #Run Updog
      #I initially used the "norm" model
      #I am also taking the ploidy from the max value in the
      updateProgressBar(session = session, id = "pb_madc", value = 40, title = "Dosage Calling in Progress")
      print('Performing Updog dosage calling')
      mout <- multidog(refmat = matrices$ref_matrix,
                       sizemat = matrices$size_matrix,
                       ploidy = as.numeric(ploidy),
                       model = model_select,
                       nc = cores)
      #Status
      updateProgressBar(session = session, id = "pb_madc", value = 100, title = "Finished")
      mout
    })

    # Only make available the download button when analysis is finished
    observe({
      if (!is.null(updog_out())) {
        Sys.sleep(1)
        # enable the download button
        enable("download_updog_vcf")
      } else {
        disable("download_updog_vcf")
      }
    })

    output$download_updog_vcf <- downloadHandler(
      filename = function() {
        paste0(input$output_name, ".vcf.gz")
      },
      content = function(file) {
        #Save Updog output as VCF file
        temp <- tempfile()
        updog2vcf(
          multidog.object = updog_out(),
          output.file = temp,
          updog_version = packageVersion("updog"),
          compress = TRUE
        )

        # Move the file to the path specified by 'file'
        file.copy(paste0(temp, ".vcf.gz"), file)

        # Delete the temporary file
        unlink(temp)
      })
  })
}

## To be copied in the UI
# mod_DosageCall_ui("DosageCall_1")

## To be copied in the server
# mod_DosageCall_server("DosageCall_1")
