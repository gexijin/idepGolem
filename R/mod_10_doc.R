#' 10_doc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_10_doc_ui <- function(id){
  ns <- NS(id)
  fluidPage(
    p("test!"),
    plotOutput(ns("test1"))
  )

}
    
#' 10_doc Server Functions
#'
#' @noRd 
mod_10_doc_server <- function(id, pre_process, idep_data, tab){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    output$test1 <- ({
      req(!is.null(pre_process$data()))
      hist(pre_process$data()[, 2])

    })
 
  })
}
    
## To be copied in the UI
# mod_09_doc_ui("09_doc_ui_1")
    
## To be copied in the server
# mod_09_doc_server("09_doc_ui_1")
