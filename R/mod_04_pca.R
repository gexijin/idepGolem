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
          condition = "input.PCA_panels != 'PCAtools Package Plots'",
          fluidRow(
            column(
              width = 9,
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
        #PCATools plot options
        conditionalPanel(
          condition = "input.PCA_panels == 'PCAtools Package Plots'",
          fluidRow(
            selectInput(inputId = ns("x_axis_pc"),
                        label = "X-Axis",
                        choices = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                        selected = "PC1"
            ),
            selectInput(inputId = ns("y_axis_pc"),
                        label = "Y-Axis",
                        choices = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                        selected = "PC2"
            ),
            #Dynamic Color and Shape options
            uiOutput(
              outputId = ns("pcatools_shape")
            ),
            uiOutput(
              outputId = ns("pcatools_color")
            ),
            #plot customization
            checkboxInput(inputId = ns("showLoadings"), label = "Show Loadings", value = FALSE),
            checkboxInput(inputId = ns("encircle"), label = "Encircle", value = FALSE),
            checkboxInput(inputId = ns("pointLabs"), label = "Point Labels", value = TRUE),
            numericInput(inputId = ns("pointSize"), label = "Point Size (Reccomded: 1-10)",value = 3.0, min = 1, max = 15)
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
          ),
          tabPanel(
            "PCAtools Package Plots",
            br(),
            plotOutput(
              outputId = ns("pcatools_biplot"),
              width = "100%",
              height = "500px"
            ),
            br(),
            br(),
            br(),
            plotOutput(
              outputId = ns("pcatools_scree"),
              width = "100%",
              height = "500px"
            )
            
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
    
    #PCAtools biplot  ---------------------
    output$pcatools_biplot <- renderPlot({
      req(!is.null(pre_process$data()))
      
      PCA_biplot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info(),
        selected_x = input$x_axis_pc,
        selected_y = input$y_axis_pc,
        encircle = input$encircle,
        showLoadings = input$showLoadings,
        pointlabs = input$pointLabs,
        point_size = input$pointSize,
        ui_color = input$selectColor,
        ui_shape = input$selectShape
      )
    }) 
    #PCAtools Scree Plot --------------------
    output$pcatools_scree <- renderPlot({
      req(!is.null(pre_process$data()))
      
      PCA_Scree(
        processed_data = pre_process$data()
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
        selectInput(inputId = ns("selectFactors2"),
                    label="Shape",
                    choices=tem,
                    selected = "Sample_Name"
                    )
      } 
    })
    
    # select color & shape for pcatools
    output$pcatools_color <- renderUI({
      req(!is.null(pre_process$data()))
      
      if (is.null(pre_process$sample_info()) )
      { return(HTML("Upload a sample info file to customize this plot.") ) }	 else { 
        selectInput(
          inputId = ns("selectColor"),
          label = "Color",
          choices = colnames(pre_process$sample_info()),
          selected = colnames(pre_process$sample_info())[1]
          )   } 
    })
    output$pcatools_shape <- renderUI({
      req(!is.null(pre_process$data()))
      
      if (is.null(pre_process$sample_info()) )
      { return(HTML("Upload a sample info file to customize this plot.") ) }	 else { 
        selectInput(
          inputId = ns("selectShape"),
          label = "Shape",
          choices = colnames(pre_process$sample_info()),
          selected = colnames(pre_process$sample_info())[1]
        )   } 
    })    
    
    
    
    
    #Pathway Analysis ------------------
    
  })
}

## To be copied in the UI
# mod_05_pca_ui("05_pca_ui_1")

## To be copied in the server
# mod_05_pca_server("05_pca_ui_1")
