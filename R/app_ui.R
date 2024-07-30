#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom bs4Dash bs4Badge bs4DashSidebar bs4DashNavbar bs4DashPage sidebarMenu menuItem menuSubItem dashboardBody tabItems tabItem box
#' @importFrom shinydisconnect disconnectMessage
#' @import shinyWidgets
#'
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    bs4DashPage(
      skin = "black",
      bs4DashNavbar(title = tagList(
        tags$img(src = 'www/BIG_R_logo.png', height = '40', width = '50'),
      )
      ),
      bs4DashSidebar(
        skin="light", status = "info",
        sidebarMenu(
          flat = FALSE,
          menuItem("Home", tabName = "welcome", icon = icon("house")),
          menuItem("Dosage Calling", tabName = "dosage_calling", icon = icon("diagram-next"),
                   menuSubItem("DArT Report2VCF", tabName = "dosage2vcf", icon = icon("share-from-square")),
                   menuSubItem("Updog Dosage Calling", tabName = "updog", icon = icon("list-ol")),
                   menuSubItem("VCF Filtering", tabName = "filtering", icon = icon("filter"))),
          menuItem("Population Structure", tabName = "pop_struct", icon = icon("layer-group"),
                   menuSubItem("PCA", tabName = "pca", icon = icon("chart-simple")),
                   menuSubItem("DAPC", tabName = "dapc", icon = icon("circle-nodes"))),
          menuItem("Genomic Diversity", tabName = "diversity", icon = icon("chart-pie")),
          menuItem("GWAS", tabName = "gwas", icon = icon("think-peaks")),
          #menuItem("QTL Analysis", tabName = "qtl", icon = icon("chart-area")),
          menuItem(
            span("Genomic Prediction", bs4Badge("beta", position = "right", color = "success")),
            tabName = "prediction",
            icon = icon("right-left")),
          menuItem("Source Code", icon = icon("circle-info"), href = "https://www.github.com/Breeding-Insight/Genomics_Shiny_App"),
          menuItem(
            span("Job Queue", bs4Badge("demo", position = "right", color = "warning")),
            tabName = "slurm",
            icon = icon("clock")),
          menuItem("Help", tabName = "help", icon = icon("circle-question"))
        )
      ),
      dashboardBody(
        disconnectMessage(), #Adds generic error message for any error if not already accounted for
        tabItems(
          tabItem(tabName = "welcome",
                  # Add welcome content here
                  fluidRow(
                    box(
                      title = "General Info", width = 4,
                      "The app is currently under development",
                      style = "overflow-y: auto; height: 400px"
                    ),

                    box(
                      title = "Logo", width = 8,
                      img(src="www/BreedingInsight_Primary_RGBColor_400px.png"),
                      style = "overflow-y: auto; height: 400px")
                  )
          ),
          tabItem(
            tabName = "filtering", mod_Filtering_ui("Filtering_1")
          ),
          tabItem(
            tabName = "updog", mod_DosageCall_ui("DosageCall_1")
          ),
          tabItem(
            tabName = "dosage2vcf", mod_dosage2vcf_ui("dosage2vcf_1")
          ),
          tabItem(
            tabName = "pca", mod_PCA_ui("PCA_1")
          ),
          tabItem(
            tabName = "dapc", mod_dapc_ui("dapc_1")
          ),
          tabItem(
            tabName = "gwas", mod_gwas_ui("gwas_1")
          ),
          tabItem(
            tabName = "diversity", mod_diversity_ui("diversity_1")
          ),
          tabItem(
            tabName = "prediction", mod_GS_ui("GS_1")
          ),
          tabItem(
            tabName = "slurm", mod_slurm_ui("slurm_1")
          ),
          tabItem(
            tabName = "help", mod_help_ui("help_1")
          )
        )
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "BIGapp"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
