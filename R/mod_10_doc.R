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
  tabPanel(
    "About",
    fluidPage(
      h3("Citation:"),
      p("Ge, Son & Yao, iDEP: an integrated web application for differential 
        expression and pathway analysis of RNA-Seq data, ", 
        a("BMC Bioinformatics 19:1-24, 2018", 
           href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6",
           target="_blank"
        )
      )
    )
  )

}
    
#' 10_doc Server Functions
#'
#' @noRd 
mod_10_doc_server <- function(id, pre_process, idep_data, tab){
  moduleServer( id, function(input, output, session){
    ns <- session$ns 
  })
}
    
## To be copied in the UI
# mod_09_doc_ui("10_doc_ui_1")
    
## To be copied in the server
# mod_09_doc_server("10_doc_ui_1")
