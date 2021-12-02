#' 04_pca UI Function
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
        conditionalPanel(
          condition = "input.PCA_panels == 'Principal Component Analysis'",
          fluidRow( 
            column(
              width = 12,
              selectInput(
                inputId = ns("PCAx"),
                "Principal component for x-axis",
                choices = 1:5,
                selected = 1
              )
            ),
            column(
              width = 12,
              selectInput(
                inputId = ns("PCAy"),
                "Principal component for y-axis",
                choices = 1:5,
                selected = 2
              )
            )
          ),
          ns=ns
        ),
        conditionalPanel(
          condition = "input.PCA_panels == 't-SNE'",
          fluidRow( 
            actionButton(inputId = ns("seedTSNE"), label = "Re-calculate from new seed")
          ),
          ns=ns
        )
      ),
      mainPanel(
        tabsetPanel(
          id = ns("PCA_panels"),
          tabPanel(
            title="Principal Component Analysis",
            br(),
            plotOutput(
              outputId = ns("pca_plot_obj"),
              width = "100%",
              height = "500px"
            )
          ),
          tabPanel(
            "Multi-Dimensional Scaling",
            br(),
            plotOutput(
              outputId = ns("mds_plot_obj"),
              width = "100%",
              height = "500px"
            )
          ),
          tabPanel(
            "t-SNE",
            br(),
            plotOutput(
              outputId = ns("t_sne"),
              width = "100%",
              height = "500px"
            ),
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
mod_04_pca_server <- function(id, pre_process, idep_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # PCA plot ------------
    output$pca_plot_obj <- renderPlot({
      req(!is.null(pre_process$data()))
      
      PCA_plot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info(),
        PCAx = input$PCAx,
        PCAy = input$PCAy
      )
    })
    
    # t_SNE plot -----------------
    output$t_sne <- renderPlot({
      req(!is.null(pre_process$data()))
      
      input$seedTSNE

      t_SNE_plot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info()
      )
    })
    
    
    
    # MDS plot ------------
    output$mds_plot_obj <- renderPlot({
      req(!is.null(pre_process$data()))
      
      MDS_plot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info()
        
      )
    })
    
    
    #Pathway Analysis ------------------
    
  })
}

## To be copied in the UI
# mod_05_pca_ui("05_pca_ui_1")

## To be copied in the server
# mod_05_pca_server("05_pca_ui_1")
