#' 03_heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_03_heatmap_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' 03_heatmap Server Functions
#'
#' @noRd 
mod_03_heatmap_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_03_heatmap_ui("03_heatmap_ui_1")
    
## To be copied in the server
# mod_03_heatmap_server("03_heatmap_ui_1")
