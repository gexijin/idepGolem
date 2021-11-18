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
        
        fluidRow( 
          column(4, selectInput(inputId = ns("PCAx"), "Principal component for x-axis", choices = 1:5, selected = 1))  
          ,column(4, selectInput(inputId = ns("PCAy"), "Principal component for y-axis", choices = 1:5, selected = 2) )
        
        ),
        
        textOutput(ns("test"))
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            title="Principal Component Analysis",
            br(),
            plotOutput(
              outputId = ns("PCA_plot_call"),
              width = "100%",
              height = "500px"
            )
          ),
          tabPanel(
            "Multi-Dimensional Scaling",
            br(),
            plotOutput(
              outputId = ns("MDS"),
              width = "100%",
              height = "500px"
            )
          ),
          tabPanel(
            "t-SNE",
            br(),
            plotOutput(
              outputId = ns("tSNE"),
              width = "100%",
              height = "500px"
            ),
            actionButton("seedTSNE", "Re-calculate using different random numbers")
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
mod_04_pca_server <- function(id, pre_process) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    #output$test <- renderText({paste0(pre_process$data()[1,1])})
    
    
    # PCA plot ------------
    output$PCA_plot_call <- renderPlot({
      req(!is.null(pre_process$data()))
      #browser()
      
      PCA_plot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info(),
        PCAx = input$PCAx,
        PCAy = input$PCAy
        
      )
    })
    
    # Update PCA Input ---------
    # observe({
    # #  req(tab() == "Principal Component Analysis")
    #   req(!is.null(pre_process$data()))
    # 
    #   updateSelectInput(
    #     inputId = "PCAx"
    #     
    #   )
    # })
    
    #t_SNE plot -----------------
    output$tSNE <- renderPlot({
      req(!is.null(pre_process$data()))

      t_SNE_plot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info()
      )
      
      
    })
    
    
    
    # MDS plot ------------
    output$MDS <- renderPlot({
      req(!is.null(pre_process$data()))
      #browser()
      
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
