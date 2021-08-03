#' 07_deg2 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_07_deg2_ui <- function(id) {
  ns <- NS(id)
  tagList()
}

#' 07_deg2 Server Functions
#'
#' @noRd
mod_07_deg2_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
  })
}

## To be copied in the UI
# mod_07_deg2_ui("07_deg2_ui_1")

## To be copied in the server
# mod_07_deg2_server("07_deg2_ui_1")
