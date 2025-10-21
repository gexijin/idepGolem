#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  # === TIMING: UI construction start ===
  .ui_start_time <- Sys.time()
  cat(sprintf("[%s] UI: Construction started\n", format(.ui_start_time, "%H:%M:%OS3")))

  .ui_result <- tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic

    navbarPage(
      "iDEP",
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
      mod_10_doc_ui(id = "doc")
    ),
    tags$head(includeHTML(app_sys("app/www/google_analytics_GA4.html")))
  )

  .ui_end_time <- Sys.time()
  .total_ui_time <- as.numeric(difftime(.ui_end_time, .ui_start_time, units = "secs"))
  cat(sprintf("[%s] UI: Construction completed (%.3fs)\n",
              format(.ui_end_time, "%H:%M:%OS3"), .total_ui_time))
  cat(sprintf("========================================\n"))
  cat(sprintf("TOTAL UI CONSTRUCTION: %.3f seconds\n", .total_ui_time))
  cat(sprintf("========================================\n"))

  return(.ui_result)
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
      app_title = "iDEP"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
