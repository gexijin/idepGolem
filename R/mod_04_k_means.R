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
  tabPanel(
    "Test",
    sidebarLayout(
      sidebarPanel(
        actionButton(ns("test"), "Show Plot")
      ),
      mainPanel(
        shinyBS::bsModal(
          id = "modal_example",
          title = "Example",
          trigger = ns("test"),
          size = "large",
          plotOutput(ns("cars_modal"))
        )
      )
    )
  )
}

#' 04_k_means Server Functions
#'
#' @noRd
mod_04_k_means_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$cars_modal <- renderPlot(
      pairs(
        ~mpg+disp+drat+wt,
        data=mtcars,
        main="Simple Scatterplot Matrix"
      )
    )
  })
}

## To be copied in the UI
# mod_04_k_means_ui("04_k_means_ui_1")

## To be copied in the server
# mod_04_k_means_server("04_k_means_ui_1")
