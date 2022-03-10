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
            ),
          ),
          ns=ns
        ),
        conditionalPanel(
          condition = "1==1",
          fluidRow(
            column(
              width = 12,
              uiOutput(
                outputId = ns("listFactors2")
              ),
              uiOutput(
                outputId = ns("listFactors1")
              )
            ),
          ),
          ns= ns
        ),
        conditionalPanel(
          condition = "input.PCA_panels == 't-SNE'",
          fluidRow( 
            actionButton(inputId = ns("seedTSNE"), label = "Re-calculate from new seed")
          ),
          
          ns=ns
        ),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/pca/",
          target = "_blank"
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
            ),
            br(),
            shiny::textOutput(
              outputId = ns("pc_correlation")
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
          )
          # tabPanel(
          #   "Pathway Analysis of PCA",
          #   NULL
          # )
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
        PCAy = input$PCAy,
        selected_shape = input$selectFactors2,
        selected_color = input$selectFactors1
      )
    })
    
    # PC Factor Correlation ---------
    output$pc_correlation <- renderText({
      req(!is.null(pre_process$data()))
      pc_factor_correlation(
        data = pre_process$data(),
        sample_info = pre_process$sample_info()
      )
    })
    
    # t_SNE plot -----------------
    output$t_sne <- renderPlot({
      req(!is.null(pre_process$data()))
      
      input$seedTSNE

      t_SNE_plot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info(),
        selected_shape = input$selectFactors2,
        selected_color = input$selectFactors1
      )
    })
    
    
    
    # MDS plot ------------
    output$mds_plot_obj <- renderPlot({
      req(!is.null(pre_process$data()))
      
      MDS_plot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info(),
        selected_shape = input$selectFactors2,
        selected_color = input$selectFactors1
        
      )
    })
    # select color
    output$listFactors1 <- renderUI({
      req(!is.null(pre_process$data()))

      if (is.null(pre_process$sample_info()) )
      { return(HTML("Upload a sample info file to customize this plot.") ) }	 else { 
        selectInput(
          inputId = ns("selectFactors1"),
          label = "Color: ",
          choices = c( colnames(pre_process$sample_info()), "Sample_Name")
                    , selected = "Sample_Name")   } 
    })
    
    #select shape
    output$listFactors2 <- renderUI({
      req(!is.null(pre_process$data()))

      
      if (is.null(pre_process$sample_info()) )
      { return(NULL) }
      else { 
        tem <- c( colnames(pre_process$sample_info()), "Sample_Name")
        #if(length(tem)>1) { tem2 = tem[1]; tem[1] <- tem[2]; tem[1] = tem2; } # swap 2nd factor with first
        selectInput(inputId = ns("selectFactors2"),
                    label="Shape:",
                    choices=tem,
                    selected = "Sample_Name"
                    )
      } 
    })
    
    
    
    
    
    #Pathway Analysis ------------------
    
  })
}

## To be copied in the UI
# mod_05_pca_ui("05_pca_ui_1")

## To be copied in the server
# mod_05_pca_server("05_pca_ui_1")
