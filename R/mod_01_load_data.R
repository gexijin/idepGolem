#' 01_load_data UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
### testing something, come back to later
mod_01_load_data_ui <- function(id) {
  ns <- shiny::NS(id)
  tabPanel("Load Data",
    sidebarLayout(
      sidebarPanel(
        
        # Test Print 
        textOutput(ns("test")),
        
        # Button to load demo dataset ----------
        # Manually namespace the goButton with id in module call
        actionButton(ns("go_button"), "Click here to load demo data"),
        tags$head(tags$style(
          "#test-go_button{color: red;
          font-size: 16px;
          font-style: italic;}"
        )),
        h5(" and just click the tabs for some magic!", style = "color:red"),
        
        # Reset Button -----------
        p(HTML(
          "<div align=\"right\"><A HREF=\"javascript:history.go(0)\"
          >Reset</A></div>" 
        )),
        
        # Species Match Drop Down ------------
        strong("1. Select or search for your species."),
        selectizeInput(
          ns('select_org'), 
          label    = NULL,
          choices  = " ",
          multiple = TRUE,
          options  = list(
            maxItems     = 1,               
            placeholder  = 'Best matching species',
            onInitialize = I('function() { this.setValue(""); }')
        )),
        
        # Conditional .GMT file input bar ----------
        conditionalPanel(
          "input.selectOrg == 'NEW'",
          fileInput(
            ns('gmt_file'), 
            label = 'Upload a geneset .GMT file for enrichment analysis 
                    (optional)',
            accept = c(
              'text/csv',
              'text/comma-separated-values',
              'text/tab-separated-values',
              'text/plain',
              '.csv',
              '.tsv',
              '.gmt'
        ))),
        
        # Buttons for data file format ----------
        radioButtons(
          ns("data_file_format"), 
          label = "2. Choose data type", 
          choices = list(
            "Read counts data (recommended)" = 1, 
            "Normalized expression values (RNA-seq FPKM, microarray, etc.)" = 2,
            "Fold-changes and corrected P values from CuffDiff or any other 
            program" = 3
          ),
          selected = 1
        )      
                         
        
        
        
      ),
      mainPanel(
        
        NULL
        
      )
    )
  )
}

#' 01_load_data Server Functions
#'
#' @noRd
### testing something, come back to later
mod_01_load_data_server <- function(id, idep_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
      
    output$test <- renderText(
      "Hello World"
      )
    
    
    }
  )
}

## To be copied in the UI
# mod_01_load_data_ui("01_load_data_ui_1")

## To be copied in the server
# mod_01_load_data_server("01_load_data_ui_1")
