#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  idep_data <- get_idep_data()
  # Tab Variable to control reactivity
  tab <- reactive(input$navbar)

  load_data <- mod_01_load_data_server(
    id = "load_data",
    idep_data = idep_data,
    tab = tab
  )
  pre_process <- mod_02_pre_process_server(
    id = "pre_process",
    load_data = load_data,
    tab = tab
  )
  mod_03_clustering_server(
    id = "clustering",
    pre_process = pre_process,
    idep_data = idep_data,
    tab = tab
  )
  mod_04_pca_server(
    id = "pca"
  )
  deg <- mod_05_deg_server(
    id = "deg",
    pre_process = pre_process
  )
  mod_06_pathway_server(
    id = "pathway",
    pre_process = pre_process,
    deg = deg,
    idep_data = idep_data,
    tab = tab
  )
}
