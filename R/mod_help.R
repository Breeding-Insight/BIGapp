#' help UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_help_ui <- function(id){
  ns <- NS(id)
  tagList(
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
}

#' help Server Functions
#'
#' @noRd
mod_help_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_help_ui("help_1")

## To be copied in the server
# mod_help_server("help_1")
