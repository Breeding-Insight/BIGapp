#' help UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList includeMarkdown
mod_help_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      column(width=12),
      column(width=12,
             box(title="DArT Report2VCF", id = "DArT_Report2VCF_box",width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
                 "**Draft**This tab is designed to convert the DArT Dose Report and Counts files to a VCF file. **DArT Website**",
                 br(), br(),
                 bs4Dash::tabsetPanel(id = "DArT_Report2VCF_tabset",
                                      tabPanel("Parameters description", value = "DArT_Report2VCF_par",
                                               includeMarkdown(system.file("help_files/DArT_Report2VCF_par.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("Results description", value = "DArT_Report2VCF_results",
                                               includeMarkdown(system.file("help_files/DArT_Report2VCF_res.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("How to cite", value = "DArT_Report2VCF_cite",
                                               includeMarkdown(system.file("help_files/DArT_Report2VCF_cite.Rmd", package = "BIGapp"))
                                      ))
             ),
             box(title="Updog Dosage Calling", id = "Updog_Dosage_Calling_box",width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
                 "This tab is designed to handle the process of dosage calling in genomic data. Dosage calling is essential for determining the number of copies of a particular allele at each genomic location.",
                 br(), br(),
                 bs4Dash::tabsetPanel(id = "Updog_Dosage_Calling_tabset",
                                      tabPanel("Parameters description", value = "Updog_Dosage_Calling_par", br(),
                                               includeMarkdown(system.file("help_files/Updog_Dosage_Calling_par.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("Results description", value = "Updog_Dosage_Calling_results", br(),
                                               includeMarkdown(system.file("help_files/Updog_Dosage_Calling_res.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("How to cite", value = "Updog_Dosage_Calling_cite", br(),
                                               includeMarkdown(system.file("help_files/Updog_Dosage_Calling_cite.Rmd", package = "BIGapp"))
                                      ))
             ),
             box(title="VCF Filtering", id = "VCF_Filtering_box",width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
                 bs4Dash::tabsetPanel(id = "VCF_Filtering_tabset",
                                      tabPanel("Parameters description", value = "VCF_Filtering_par",
                                               includeMarkdown(system.file("help_files/VCF_Filtering_par.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("Results description", value = "VCF_Filtering_results",
                                               includeMarkdown(system.file("help_files/VCF_Filtering_par.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("How to cite", value = "Updog_Dosage_Calling_cite",
                                               includeMarkdown(system.file("help_files/VCF_Filtering_cite.Rmd", package = "BIGapp"))
                                      ))
             ),
             box(title="PCA", id = "PCA_box",width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
                 bs4Dash::tabsetPanel(id = "PCA_tabset",
                                      tabPanel("Parameters description", value = "PCA_par",
                                               includeMarkdown(system.file("help_files/PCA_par.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("Results description", value = "PCA_results",
                                               includeMarkdown(system.file("help_files/PCA_res.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("How to cite", value = "PCA_cite",
                                               includeMarkdown(system.file("help_files/PCA_cite.Rmd", package = "BIGapp"))
                                      ))
             ),
             box(title="DAPC", id = "DAPC_box",width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
                 bs4Dash::tabsetPanel(id = "DAPC_tabset",
                                      tabPanel("Parameters description", value = "DAPC_par",
                                               includeMarkdown(system.file("help_files/DAPC_par.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("Results description", value = "DAPC_results",
                                               includeMarkdown(system.file("help_files/DAPC_res.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("How to cite", value = "DAPC_cite",
                                               includeMarkdown(system.file("help_files/DAPC_cite.Rmd", package = "BIGapp"))
                                      ))
             ),
             box(title="Genomic Diversity", id = "Genomic_Diversity_box",width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
                 "**Draft**This tab is dedicated to analyzing genomic diversity within the population. It calculates various diversity metrics such as heterozygosity and minor allele frequency (MAF). The app includes functionalities to visualize these metrics through histograms and other plots. Users can download the calculated diversity metrics as CSV files. This tab helps in understanding the genetic variability and distribution of alleles within the population.",
                 br(), br(),
                 bs4Dash::tabsetPanel(id = "Genomic_Diversity_tabset",
                                      tabPanel("Parameters description", value = "Genomic_Diversity_par",
                                               includeMarkdown(system.file("help_files/Genomic_Diversity_par.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("Results description", value = "Genomic_Diversity_results",
                                               includeMarkdown(system.file("help_files/Genomic_Diversity_res.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("How to cite", value = "Genomic_Diversity_cite",
                                               includeMarkdown(system.file("help_files/Genomic_Diversity_cite.Rmd", package = "BIGapp"))
                                      ))
             ),
             box(title="GWAS", id = "GWAS_box",width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
                 "The tab is for conducting Genome-Wide Association Studies (GWAS) to identify associations between genetic variants and traits of interest. Users can input phenotypic data and specify parameters for the GWAS analysis. The app performs statistical tests to identify significant associations between SNPs and traits, and visualizes the results using Manhattan plots and Q-Q plots. The tab helps in identifying potential genetic markers linked to specific traits. GWASpoly package is used to perform the analysis.",
                 br(), br(),
                 bs4Dash::tabsetPanel(id = "GWAS_tabset",
                                      tabPanel("Parameters description", value = "GWAS_par", br(),
                                               includeMarkdown(system.file("help_files/GWAS_par.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("Results description", value = "GWAS_results", br(),
                                               includeMarkdown(system.file("help_files/GWAS_res.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("How to cite", value = "GWAS_cite", br(),
                                               includeMarkdown(system.file("help_files/GWAS_cite.Rmd", package = "BIGapp"))
                                      ))
             ),
             box(title="Predictive Ability", id = "Predictive_Ability_box",width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
                 bs4Dash::tabsetPanel(id = "Predictive_Ability_tabset",
                                      tabPanel("Parameters description", value = "Predictive_Ability_par",
                                               includeMarkdown(system.file("help_files/Predictive_Ability_par.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("Results description", value = "Predictive_Ability_results",
                                               includeMarkdown(system.file("help_files/Predictive_Ability_res.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("How to cite", value = "Predictive_Ability_cite",
                                               includeMarkdown(system.file("help_files/Predictive_Ability_cite.Rmd", package = "BIGapp"))
                                      ))
             ),
             box(title="Genomic Prediction", id = "Genomic_Prediction_box",width = 12, collapsible = TRUE, collapsed = TRUE, status = "info", solidHeader = TRUE,
                 bs4Dash::tabsetPanel(id = "Genomic_Prediction_tabset",
                                      tabPanel("Parameters description", value = "Genomic_Prediction_par",
                                               includeMarkdown(system.file("help_files/Genomic_Prediction_par.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("Results description", value = "Genomic_Prediction_results",
                                               includeMarkdown(system.file("help_files/Genomic_Prediction_res.Rmd", package = "BIGapp"))
                                      ),
                                      tabPanel("How to cite", value = "Genomic_Prediction_cite",
                                               includeMarkdown(system.file("help_files/Genomic_Prediction_cite.Rmd", package = "BIGapp"))
                                      ))
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
mod_help_server <- function(input, output, session, parent_session){

  ns <- session$ns

}

## To be copied in the UI
# mod_help_ui("help_1")

## To be copied in the server
# mod_help_server("help_1")
