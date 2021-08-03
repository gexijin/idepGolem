#' 04_k_means UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_04_k_means_ui <- function(id) {
  ns <- NS(id)
  tagList()
}

#' 04_k_means Server Functions
#'
#' @noRd
mod_04_k_means_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
  })
}

## To be copied in the UI
# mod_04_k_means_ui("04_k_means_ui_1")

## To be copied in the server
# mod_04_k_means_server("04_k_means_ui_1")
