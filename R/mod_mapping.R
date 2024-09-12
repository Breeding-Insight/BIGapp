#' mapping UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_mapping_ui <- function(id){
  ns <- NS(id)
  tagList(

  )
}

#' mapping Server Functions
#'
#' @noRd
mod_mapping_server <- function(input, output, session, parent_session){

  ns <- session$ns

}

## To be copied in the UI
# mod_mapping_ui("mapping_1")

## To be copied in the server
# mod_mapping_server("mapping_1")
