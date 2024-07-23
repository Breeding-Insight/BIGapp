#' GenoCall UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_GenoCall_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' GenoCall Server Functions
#' 
#' @import vcfR
#' @import updog
#' 
#' @noRd 
mod_GenoCall_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    
    
 
  })
}
    
## To be copied in the UI
# mod_GenoCall_ui("GenoCall_1")
    
## To be copied in the server
# mod_GenoCall_server("GenoCall_1")
