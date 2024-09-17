#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic

  ##Add server configurations
  options(shiny.maxRequestSize = 10000 * 1024^2)  # Set maximum upload size to 10GB
  #shiny.maxRequestSize = 10000 * 1024^2; # 10 GB <- This is for a future limit when using BI's server remotely

  callModule(mod_DosageCall_server,
             "DosageCall_1",
             parent_session = session)
  callModule(mod_Filtering_server,
             "Filtering_1",
             parent_session = session)
  callModule(mod_dosage2vcf_server,
             "dosage2vcf_1",
             parent_session = session)
  callModule(mod_PCA_server,
             "PCA_1",
             parent_session = session)
  callModule(mod_dapc_server,
             "dapc_1",
             parent_session = session)
  callModule(mod_gwas_server,
             "gwas_1",
             parent_session = session)
  callModule(mod_diversity_server,
             "diversity_1",
             parent_session = session)
  callModule(mod_GS_server,
             "GS_1",
             parent_session = session)
  callModule(mod_GSAcc_server,
            "GSAcc_1",
            parent_session = session)
  callModule(mod_slurm_server,
             "slurm_1",
             parent_session = session)

  # mod_DosageCall_server("DosageCall_1")
  # mod_Filtering_server("Filtering_1")
  # mod_dosage2vcf_server("dosage2vcf_1")
  # mod_PCA_server("PCA_1")
  # mod_dapc_server("dapc_1")
  # mod_gwas_server("gwas_1")
  # mod_diversity_server("diversity_1")
  # mod_GS_server("GS_1")
  # mod_GSAcc_server("GSAcc_1")
  # mod_slurm_server("slurm_1")

  #Session info popup
  observeEvent(input$session_info_button, {
    showModal(modalDialog(
      title = "Session Information",
      size = "l",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close"),
        downloadButton("download_session_info", "Download")
      ),
      pre(
        paste(capture.output(sessionInfo()), collapse = "\n")
      )
    ))
  })

  #Download Session Info
  output$download_session_info <- downloadHandler(
    filename = function() {
      paste("session_info_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      writeLines(paste(capture.output(sessionInfo()), collapse = "\n"), file)
    }
  )

}
