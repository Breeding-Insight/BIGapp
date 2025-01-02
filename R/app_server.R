#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom httr GET content status_code
#' @importFrom curl new_handle curl_fetch_memory
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
  
  #Check for updates from GitHub for BIGapp
  get_latest_github_commit <- function(repo, owner) {
    url <- paste0("https://api.github.com/repos/", owner, "/", repo, "/releases/latest")
    response <- GET(url)
    content <- content(response, "parsed")
    
    if (status_code(response) == 200) {
      tag_name <- content$tag_name
      clean_tag_name <- sub("-.*", "", tag_name)
      clean_tag_name <- sub("v", "", clean_tag_name)
      return(clean_tag_name)
    } else {
      return(NULL)
    }
  }
  
  is_internet_connected <- function() {
    handle <- new_handle()
    success <- tryCatch({
      curl_fetch_memory("https://www.google.com", handle = handle)
      TRUE
    }, error = function(e) {
      FALSE
    })
    return(success)
  }
  
  observeEvent(input$updates_info_button, {
    # Check internet connectivity
    if (!is_internet_connected()) {
      # Display internet connectivity issues message
      showModal(modalDialog(
        title = "No Internet Connection",
        easyClose = TRUE,
        footer = tagList(
          modalButton("Close")
        ),
        "Please check your internet connection and try again."
      ))
      return()
    }
    
    package_name <- "BIGapp"
    repo_name <- "BIGapp" # GitHub repo name
    repo_owner <- "Breeding-Insight" # User or organization name
    
    # Get the installed version
    installed_version <- as.character(packageVersion(package_name))
    
    # Get the latest version from GitHub (can be tag version or latest commit)
    latest_commit <- get_latest_github_commit(repo_name, repo_owner)
    
    # Compare versions and prepare message
    if (latest_commit > installed_version) {
      update_status <- "A new version is available. Please update your package."
      # Prepare styled HTML text for the modal
      message_html <- paste(
        "Installed version:", installed_version, "<br>",
        #"Latest version commit SHA:", latest_commit, "<br>",
        "<span>A new version is available on GitHub!</span><br>",
        "<span style='color: red;'>Please update your package.</span>"
      )
    } else {
      update_status <- "Your package is up-to-date!"
      # Prepare non-styled text for no update needed
      message_html <- paste(
        "Installed version:", installed_version, "<br>",
        #"Latest version commit SHA:", latest_commit, "<br>",
        update_status
      )
    }
    
    # Display message in a Shiny modal
    showModal(modalDialog(
      title = "BIGapp Updates",
      size = "m",
      easyClose = TRUE,
      footer = tagList(
        modalButton("Close")
      ),
      # Use HTML to format the message and include styling
      HTML(message_html)
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
