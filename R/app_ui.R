#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom bs4Dash bs4Badge bs4DashSidebar bs4DashNavbar bs4DashPage sidebarMenu menuItem menuSubItem dashboardBody tabItems tabItem box dashboardFooter
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
          tags$li(class = "header", style = "color: grey; margin-top: 10px; margin-bottom: 10px; padding-left: 15px;", "Menu"),
          menuItem("Home", tabName = "welcome", icon = icon("house")),
          tags$li(class = "header", style = "color: grey; margin-top: 18px; margin-bottom: 10px; padding-left: 15px;", "Genotype Processing"),
                   menuItem("DArT Report2VCF", tabName = "dosage2vcf", icon = icon("share-from-square")),
                   menuItem("Updog Dosage Calling", tabName = "updog", icon = icon("list-ol")),
                   menuItem("VCF Filtering", tabName = "filtering", icon = icon("filter")),
          tags$li(class = "header", style = "color: grey; margin-top: 18px; margin-bottom: 10px; padding-left: 15px;", "Population Structure"),
                   menuItem("PCA", tabName = "pca", icon = icon("chart-simple")),
                   menuItem("DAPC", tabName = "dapc", icon = icon("circle-nodes")),
          tags$li(class = "header", style = "color: grey; margin-top: 18px; margin-bottom: 10px; padding-left: 15px;", "Summary Metrics"),
          menuItem("Genomic Diversity", tabName = "diversity", icon = icon("chart-pie")),
          tags$li(class = "header", style = "color: grey; margin-top: 18px; margin-bottom: 10px; padding-left: 15px;", "GWAS"),
          menuItem("GWASpoly", tabName = "gwas", icon = icon("think-peaks")),
          tags$li(class = "header", style = "color: grey; margin-top: 18px; margin-bottom: 10px; padding-left: 15px;", "Genomic Selection"),
                  menuItem(
                    span("Prediction Accuracy", bs4Badge("beta", position = "right", color = "success")),
                           tabName = "prediction_accuracy",
                           icon = icon("right-left")),
                  menuItem(
                    span("Genomic Prediction", bs4Badge("beta", position = "right", color = "success")),
                      tabName = "prediction",
                      icon = icon("angles-right")),
          tags$li(class = "header", style = "color: grey; margin-top: 18px; margin-bottom: 10px; padding-left: 15px;", "Information"),
          menuItem("Source Code", icon = icon("circle-info"), href = "https://www.github.com/Breeding-Insight/Genomics_Shiny_App"),
          menuItem(
            span("Job Queue", bs4Badge("demo", position = "right", color = "warning")),
            tabName = "slurm",
            icon = icon("clock")),
          menuItem("Help", tabName = "help", icon = icon("circle-question"))
        )
      ),
      footer = dashboardFooter(
        right = div(
          style = "display: flex; align-items: center;",  # Align text and images horizontally
          div(
            style = "display: flex; flex-direction: column; margin-right: 15px; text-align: right;",
            div("© 2024 Breeding Insight"),
            div("Funded by USDA through Cornell University")
          ),
          div(
            a(
              img(src = "www/usda-logo-color.png", height = "45px"),
              style = "margin-right: 15px;"
            ),
            a(
              img(src = "www/cornell_seal_simple_web_b31b1b.png", height = "45px")
            )
          )
        ),
        left = div(
          style = "display: flex; align-items: center; height: 100%;",  # Center the version text vertically
          "v0.5.0")
      ),
      dashboardBody(
        disconnectMessage(), #Adds generic error message for any error if not already accounted for
        tags$style(
          HTML(
            ".main-footer {
            background-color: white;
            color: grey;
            height: 65px;
            padding-top: 5px;
            padding-bottom: 5px;
          }
          .main-footer a {
            color: grey;
          }"
          )
        ),
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
            tabName = "prediction_accuracy", mod_GSAcc_ui("GSAcc_1")
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
