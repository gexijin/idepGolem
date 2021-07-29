#' 06_deg1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_06_deg1_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' 06_deg1 Server Functions
#'
#' @noRd 
mod_06_deg1_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_06_deg1_ui("06_deg1_ui_1")
    
## To be copied in the server
# mod_06_deg1_server("06_deg1_ui_1")
