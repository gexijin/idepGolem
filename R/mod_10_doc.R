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
    "R",
    sidebarLayout(
      sidebarPanel(
        h5(
          "Test new tab",
        ),
        numericInput(
          inputId = ns("sampleID"), 
          label = h5("sample ID"), 
          min = 1, 
          max = 5, 
          value = 1
        ),
      ),
      mainPanel(
        tabsetPanel(
          id = ns("test_tabs"),
          tabPanel(
            tile = "Histogram",
            p("test!"),
            plotOutput(
              outputId = ns("test1"),
              width = "100%",
              height = "500px"
            )
          )
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

    output$test1 <- renderPlot({
      req(!is.null(pre_process$data()))
      doc_hist(
        processed_data = pre_process$data(),
        sampleID = input$sampleID )
    })
 
  })
}
    
## To be copied in the UI
# mod_09_doc_ui("09_doc_ui_1")
    
## To be copied in the server
# mod_09_doc_server("09_doc_ui_1")
