#' Home UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shinyjs enable disable useShinyjs
#' @importFrom shiny NS tagList
#' @importFrom bs4Dash renderValueBox valueBox
#'
#'
mod_Home_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        column(width = 4,
               box(
                 title = "Breeding Insight Genomics App", status = "info", solidHeader = FALSE, width = 12, collapsible = FALSE,
                 HTML(
                   "<strong>The app is under development</strong>
              <p>Breeding Insight provides bioinformatic processing support for our external collaborators. This R shiny app provides a web-based user friendly way for users to analyze genomic data without needing to use command-line tools.</p>

              <p><b>Supported Analyses</b></p>
              Initial supported analyses includes the mature genomics/bioinformatics pipelines developed within Breeding Insight:
              <ul>
                <li>Genotype processing</li>
                <li>Summary metrics</li>
                <li>Population Structure</li>
                <li>GWAS</li>
                <li>GS</li>
              </ul>"
                 ),
                 style = "overflow-y: auto; height: 500px"

               )
        ),
        column(width = 4,
               box(
                 title = "About Breeding Insight", status = "success", solidHeader = FALSE, width = 12, collapsible = FALSE,
                 HTML(
                   "We provide scientific consultation and data management software to the specialty crop and animal breeding communities.
            <ul>
              <li>Genomics</li>
              <li>Phenomics</li>
              <li>Data Management</li>
              <li>Software Tools</li>
              <li>Analysis</li>
            </ul>
            Breeding Insight is funded by the U.S. Department of Agriculture (USDA) Agricultural Research Service (ARS) through Cornell University.
            <div style='text-align: center; margin-top: 20px;'>
              <img src='www/BreedingInsight.png' alt='Breeding Insight' style='width: 85px; height: 85px;'>
            </div>"
                 ),
                 style = "overflow-y: auto; height: 500px"
               )
        ),
        column(width = 4,
               a(
                 href = "https://www.breedinginsight.org",  # Replace with your desired URL
                 target = "_blank",  # Optional: opens the link in a new tab
                 valueBox(
                   value = NULL,
                   subtitle = "Learn More About Breeding Insight",
                   icon = icon("link"),
                   color = "purple",
                   gradient = TRUE,
                   width = 11
                 ),
                 style = "text-decoration: none; color: inherit;"  # Optional: removes underline and retains original color
               ),
               a(
                 href = "https://breedinginsight.org/contact-us/",  # Replace with your desired URL
                 target = "_blank",  # Optional: opens the link in a new tab
                 valueBox(
                   value = NULL,
                   subtitle = "Contact Us",
                   icon = icon("envelope"),
                   color = "danger",
                   gradient = TRUE,
                   width = 11
                 ),
                 style = "text-decoration: none; color: inherit;"  # Optional: removes underline and retains original color
               ),
               a(
                 href = "https://scribehow.com/page/BIGapp_Tutorials__FdLsY9ZxQsi6kgT9p-U2Zg",  # Replace with your desired URL
                 target = "_blank",  # Optional: opens the link in a new tab
                 valueBox(
                   value = NULL,
                   subtitle = "BIGapp Tutorials (in-progress)",
                   icon = icon("compass"),
                   color = "warning",
                   gradient = TRUE,
                   width = 11
                 ),
                 style = "text-decoration: none; color: inherit;"  # Optional: removes underline and retains original color
               )
        )
      )
    )
  )
}

#' Home Server Functions
#'
#'
#' @noRd
mod_Home_server <- function(input, output, session, parent_session){

  ns <- session$ns

}

## To be copied in the UI
# mod_Home_ui("Home_1")

## To be copied in the server
# mod_Home_server("Home_1")
