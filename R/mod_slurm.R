#' slurm UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom DT DTOutput
#' @importFrom shiny NS tagList
#'
mod_slurm_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      fluidRow(
        box(title="Submitted Jobs", width = 12, collapsible = FALSE, status = "info", solidHeader = TRUE,
            DTOutput(ns("job_table")),style = "overflow-y: auto; height: 500px",
            div(style = "position: absolute; bottom: 20px; width: 100%;",
                actionButton(ns("run"), "Cancel Submitted Jobs",
                             style = "color: #fff; background-color: #ff0000; border-color: #ff0000;"))
        )
      )
    )
  )
}

#' slurm Server Functions
#'
#' @importFrom DT renderDT datatable
#' @noRd
mod_slurm_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    ####Command to get information from SLURM as dataframe
    #add button to ui for user to cancel their scheduled jobs...
    get_slurm_jobs <- reactive({
      # Run the SLURM squeue command to get job information
      slurm_output <- system("squeue --Format=UserName,JobName,TimeUsed --noheader", intern = TRUE)

      # Process the output
      job_list <- strsplit(slurm_output, "\\s+")

      # Create a data frame
      job_df <- do.call(rbind, job_list)
      colnames(job_df) <- c("userID", "JobType", "Duration")

      # Convert to data frame
      job_df <- as.data.frame(job_df, stringsAsFactors = FALSE)

      return(job_df)
    })

    #Display job queue to user
    output$job_table <- renderDT({
      #job_df <- get_slurm_jobs() #Use this when I get the above code working on the server
      job_df <- data.frame(userID = c("User1","User2","User3","User4","User5"),
                           JobID = c("000303","000312","000335","000342", "000348"),
                           JobType = c("Updog Dosage Calling","Updog Dosage Calling", "Updog Dosage Calling", "GWAS", "GWAS"),
                           Duration = c("Completed: Email notification sent","06:11:43", "03:31:01", "00:46:00", "Scheduled"))

      # Render the data table
      datatable(job_df, options = list(pageLength = 10))
    })
  })
}

## To be copied in the UI
# mod_slurm_ui("slurm_1")

## To be copied in the server
# mod_slurm_server("slurm_1")
