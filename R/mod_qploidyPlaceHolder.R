#' qploidy placeholder UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList includeMarkdown
mod_qploidyPlaceHolder_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      column(
        width = 12,
        h3("Qploidy not installed"),
        p("Install the Qploidy package to access this feature."),
        div(
          style = "margin-top: 12px;",
          actionButton(ns("install_qploidy"), "Install Qploidy", icon = icon("download")),
          tags$span(style="display:inline-block; width: 8px;")
        ),
        tags$hr(),
        verbatimTextOutput(ns("install_log"))
      )
    )
  )
}

#' qploidyPlaceHolder Server Functions
#'
#' @noRd
mod_qploidyPlaceHolder_server <- function(input, output, session, parent_session){
  
  ns <- session$ns
  
  output$install_log <- renderText("")
  
  observeEvent(input$install_qploidy, {
    err_msg <- NULL
    ok <- FALSE
    
    showNotification("Installing Qploidy...")
    
    # Prefer remotes (lighter than devtools)
    tryCatch({
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes")
      }
      # Adjust owner/ref if needed
      remotes::install_github("Cristianetaniguti/Qploidy", ref = "add_shiny", upgrade = "never", quiet = TRUE)
      ok <- requireNamespace("Qploidy", quietly = TRUE)
    }, error = function(e){
      err_msg <<- conditionMessage(e)
    })
    
    if (ok) {
      shiny::showNotification("Qploidy installed. Restart the app.", type = "message", duration = 8)
      output$install_log <- renderText("")
    } else {
      shiny::showNotification("Installation failed. See details below.", type = "error", duration = NULL)
      output$install_log <- renderText(if (is.null(err_msg)) "Unknown error (check server permissions/logs)." else err_msg)
    }
  })
}

## To be copied in the UI
# mod_qploidyPlaceHolder_ui("qploidyPlaceHolder_1")

## To be copied in the server
# mod_qploidyPlaceHolder_server("qploidyPlaceHolder_1")
