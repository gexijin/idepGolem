#' 05_pca UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_04_pca_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    "PCA",
    sidebarLayout(
      sidebarPanel(
        NULL
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Principal Component Analysis",
            NULL
          ),
          tabPanel(
            "Multi-Dimensional Scaling",
            NULL
          ),
          tabPanel(
            "t-SNE",
            NULL
          ),
          tabPanel(
            "Pathway Analysis of PCA",
            NULL
          )
        )
      )
    )
  )
}

#' 05_pca Server Functions
#'
#' @noRd
mod_04_pca_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
  })
}

## To be copied in the UI
# mod_05_pca_ui("05_pca_ui_1")

## To be copied in the server
# mod_05_pca_server("05_pca_ui_1")
