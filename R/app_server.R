#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  idep_data <- get_idep_data()
  tab <- reactive(input$navbar)

  load_data <- mod_01_load_data_server(
    id = "load_data",
    idep_data = idep_data
  )
  pre_process <- mod_02_pre_process_server(
    id = "pre_process",
    load_data = load_data,
    tab
  )
  mod_03_heatmap_server(id = "heatmap")
  mod_04_k_means_server("test")
}
