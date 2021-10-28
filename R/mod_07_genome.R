#' mod_07_genome UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_07_genome_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' mod_07_genome Server Functions
#'
#' @noRd 
mod_07_genome_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_07_genome_ui("mod_07_genome_ui_1")
    
## To be copied in the server
# mod_07_genome_server("mod_07_genome_ui_1")
