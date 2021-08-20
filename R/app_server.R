#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  load_data <- mod_01_load_data_server(id = "load_data")
  mod_02_pre_process_server(id = "pre_process", load_data)
}
