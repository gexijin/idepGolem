#' 01_load_data UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
### testing something, come back to later
mod_01_load_data_ui <- function(id) {
  ns <- shiny::NS(id)
  tagList(
   
  )
}

#' 01_load_data Server Functions
#'
#' @noRd
### testing something, come back to later
mod_01_load_data_server <- function(id, idep_data) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      
    }
  )
}

## To be copied in the UI
# mod_01_load_data_ui("01_load_data_ui_1")

## To be copied in the server
# mod_01_load_data_server("01_load_data_ui_1")
