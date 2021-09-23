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
    idep_data = idep_data
  )
  pre_process <- mod_02_pre_process_server(
    id = "pre_process",
    load_data = load_data,
    tab = tab
  )
  mod_03_clustering_server(
    id = "clustering",
    pre_process = pre_process,
    tab = tab
  )
  mod_05_pca_server(
    id = "pca"
  )
}
