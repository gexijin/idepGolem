#' 09_network UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_09_network_ui <- function(id){
  ns <- NS(id)
  tabPanel(
    "Network",
    sidebarLayout(
      sidebarPanel(
        h5(
          "Identify co-expression networks and sub-modules using",
          a(
            "WGCNA.",
            href = "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559",
            target = "_blank"
          ),
          "Only useful when  sample size is large(> 15)."
        ),
        numericInput(
          inputId = ns("n_genes_network"), 
          label = h5("Most variable genes to include (< 3001)"), 
          min = 10, 
          max = 3000, 
          value = 1000
        ),
        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = ns("soft_power"), 
              label = h5("Soft Threshold"), 
              min = 1, 
              max = 20, 
              value = 5
            )
          ),
          column(
            width = 6,
            numericInput(
              inputId = ns("min_module_size"),
              label = h5("Min. Module Size"),
              min = 10, 
              max = 100, 
              value = 20
            )
          )
        ),
        tags$style(
          type = "text/css",
          "#network-my_soft_power{ width:100%;   margin-top:-12px}"
        ),
        tags$style(
          type = "text/css",
          "#network-min_module_size{ width:100%;   margin-top:-12px}"
        ),
        HTML(
          "<hr style='height:1px;border:none;color:#333;background-color:#333;' />"
        ),
        htmlOutput(outputId = ns("list_wgcna_modules")),
        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = ns("edge_threshold"), 
              label = h5("Edge Threshold"), 
              min = 0, 
              max = 1, 
              value = .4, 
              step = .1
            )
          ),
          column(
            width = 6,
            numericInput(
              inputId = ns("top_genes_network"), 
              label = h5("Top genes"), 
              min = 10, 
              max = 2000, 
              value = 10, 
              step = 10
            )
          )
        ),
        tags$style(
          type = "text/css",
          "#network-edge_threshold{ width:100%;   margin-top:-12px}"
        ),
        tags$style(
          type = "text/css",
          "#network-top_genes_network{ width:100%;   margin-top:-12px}"
        ),
        actionButton(
          inputId = ns("network_layout"),
          label = "Change network layout",
          style = "float:center"
        ),
        h5("Enrichment database:"),
        htmlOutput(outputId = ns("select_go")),
        h5("The network file can be imported to", 
          a("VisANT", href = "http://visant.bu.edu/", target = "_blank"),
          " or ", 
          a("Cytoscape.", href = "http://www.cytoscape.org/", target = "_blank")
        ),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/network/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Module Plot",
            plotOutput(ns("module_plot"))
          )
        )
      )
    )
  )
}
    
#' 09_network Server Functions
#'
#' @noRd 
mod_09_network_server <- function(id, pre_process){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    wgcna <- reactive({
      req(!is.null(pre_process$data()))

      get_wgcna(
        data = pre_process$data(),
        n_genes = input$n_genes,
        soft_power = input$soft_power,
        min_module_size = input$min_module_size
      )
	  })

    output$module_plot <- renderPlot({
		})
		
  })
}
    
## To be copied in the UI
# mod_09_network_ui("09_network_ui_1")
    
## To be copied in the server
# mod_09_network_server("09_network_ui_1")
