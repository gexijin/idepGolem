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
        radioButtons(
          inputId = ns("pca_meth"),
          label = "Methods", 
          choices = list(
            "Principal Component Analysis" = 1,            
            "Multidimensional Scaling" = 3, 
            "t-SNE" = 4,
            "Pathway Analysis of PCA rotation" = 2             
          )
        ),
        conditionalPanel(
          condition = "input.pca_meth == 2",
          htmlOutput(outputId = ns("select_go_pca")),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.pca_meth != 2",
          h5("hey you suck"),
          ns = ns
        ),

        br(),
        downloadButton(
          outputId = ns("down_pca_data"),
          label = "Coordinates"
        ),
        a(h5(
          "Questions?",
          align = "right"
          ),
          href="https://idepsite.wordpress.com/pca/",
          target="_blank"
        )        
      ),
      mainPanel(
        conditionalPanel(
          condition = "input.pca_meth == 1",
          fluidRow( 
            column(
              width = 6,
              selectInput(
                inputId = ns("pca_x"),
                label = "Principal component for x-axis",
                choices = 1:5,
                selected = 1
              )
            ),
            column(
              width = 6,
              selectInput(
                inputId = ns("pca_y"),
                label = "Principal component for y-axis",
                choices = 1:5,
                selected = 2
              )
            )
          ),
          ns = ns
        ),
        plotOutput(
          outputId = ns("pca_plot"), inline=TRUE
        ),
        br(),
        conditionalPanel(
          condition = "input.pca_meth == 4", 
          actionButton(
            inputId = ns("tsne_seed"),
            label = "Re-calculate t-SNE"),
            br(),
            br(),
          ns = ns
        ),
        conditionalPanel(
          conditon = "input.pca_meth == 1 | input.pca_meth == 2 ",
          htmlOutput(ns("pca_2_factore")),
          br(),
          br(),
          ns = ns
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
