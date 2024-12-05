#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic

    navbarPage(
      "iDEP 2.01",

      id = "navbar",
      mod_01_load_data_ui(id = "load_data"),
      mod_02_pre_process_ui(id = "pre_process"),
      mod_03_clustering_ui(id = "clustering"),
      mod_04_pca_ui(id = "pca"),
      mod_05_deg_1_ui(id = "deg"),
      mod_05_deg_2_ui(id = "deg"),
      mod_06_pathway_ui(id = "pathway"),
      mod_07_genome_ui(id = "genome"),
      mod_08_bicluster_ui(id = "bicluster"),
      mod_09_network_ui(id = "network"),
      mod_10_doc_ui(id = "doc"),

      ### Hidden Policies Tabs ###
      tabPanel(title = "Privacy Policy", value = "privacy_policy", privacy_policy_content()),
      tabPanel(title = "Terms of Use", value = "terms_of_use", terms_of_use_content())
    ),
    tags$head(includeHTML(app_sys("app/www/google_analytics_GA4.html"))),
    tags$head(includeHTML(app_sys("app/www/google_analytics_golem.html"))),

    tags$footer(
      style = "position: fixed;bottom: 0;width: 100%;background-color: #f0f0f0;
        padding: 10px;text-align: center;z-index: 99;",
      span("© 2024 Orditus LLC | "),
      actionLink(inputId = "ppolicy", "Privacy Policy"),
      span(" | "),
      actionLink(inputId = "tofu", "Terms of Use"),
      span(" | "),
      a("Orditus.com", href = "https://orditus.com/", target = "_blank")
    ) # footer
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
    "www", app_sys("app/www")
  )

  tags$head(
    favicon(
      ico = "favicon",
      rel = "shortcut icon",
      resources_path = "www",
      ext = "png"
    ),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "iDEP 2.0"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
