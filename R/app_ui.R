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
      span("iDEP", style = "color: #55565B;"),
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
      lang = "en"
    ),
    mod_14_survey_ui(id = "survey"),
    tags$head(includeHTML(app_sys("app/www/google_analytics_GA4.html")))
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
      app_title = "iDEP"
    ),
    tags$style(HTML(
      "
.idep-scrollable-plot {
  overflow-x: auto;
  width: 100%;
}

.idep-scrollable-plot .shiny-plot-output {
  max-width: none;
}

.idep-scrollable-plot .shiny-plot-output img {
  display: block;
  max-width: none;
}

.navbar-default .navbar-nav > li > a {
    color: #55565B !important;
    font-weight: bold;
}
"
    )),
    tags$script(HTML(
      "
Shiny.addCustomMessageHandler('idep-set-plot-width', function(message) {
  if (!message || !message.id) { return; }
  var el = document.getElementById(message.id);
  if (!el) { return; }
  if (message.width === null || typeof message.width === 'undefined') {
    el.style.width = '';
    el.style.minWidth = '';
    return;
  }
  var width = parseInt(message.width, 10);
  if (!isFinite(width) || width <= 0) { return; }
  el.style.width = width + 'px';
  el.style.minWidth = width + 'px';
});
"
    ))
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
