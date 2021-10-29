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
        # Buttons for data file format ----------
        radioButtons(
          inputId = ns("data_file_format"),
          label = "2. Choose data type",
          choices = list(
            "Read counts data (recommended)" = 1,
            "Normalized expression values (RNA-seq FPKM, microarray, etc.)" = 2,
            "Fold-changes and corrected P values from CuffDiff or any other
             program" = 3
          ),
          selected = 1
        ),

        # Conditional panel for fold changes data file ----------
        conditionalPanel(
          condition = "input.data_file_format == 3",
          checkboxInput(
            inputId = ns("no_fdr"),
            label = "Fold-changes only, no corrected P values",
            value = FALSE
          ),
          ns = ns
        )
      ),
      mainPanel(
        NULL
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
