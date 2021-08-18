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
      
      # Load Data Panel Sidebar -----------
      sidebarPanel(
        
        # Button to load demo dataset ----------
        # Manually namespace the goButton in tag with id in module call
        actionButton(
          inputId = ns("go_button"), 
          label = "Click here to load demo data"
        ),
        tags$head(tags$style(
          "#load_data-go_button{color: red;
          font-size: 16px;
          font-style: italic;}"
          )
        ),
        h5(" and just click the tabs for some magic!", style = "color:red"),
        
        # Reset Button -----------
        p(HTML(
          "<div align=\"right\"><A HREF=\"javascript:history.go(0)\"
           >Reset</A></div>" 
          )
        ),
        
        # Species Match Drop Down ------------
        strong("1. Select or search for your species."),
        selectizeInput(
          inputId = ns("select_org"), 
          label = NULL,
          choices = " ",
          multiple = TRUE,
          options = list(
            maxItems = 1,               
            placeholder = "Best matching species",
            onInitialize = I('function() { this.setValue(""); }')
          )
        ),
        
        # Conditional .GMT file input bar ----------
        conditionalPanel(
          'input.selectOrg == "NEW"',
          fileInput(
            inputId = ns("gmt_file"), 
            label = 
              "Upload a geneset .GMT file for enrichment analysis (optional)",
            accept = c(
              "text/csv",
              "text/comma-separated-values",
              "text/tab-separated-values",
              "text/plain",
              ".csv",
              ".tsv"
            )
          ),
          ns = ns
        ),
        
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
          "input.data_file_format == 3",
          checkboxInput(
            inputId = ns("no_fdr"), 
            label = "Fold-changes only, no corrected P values", 
            value = FALSE
          ),
          ns = ns
        ),
        
        # Expression data file input ----------           
        fileInput(
          inputId = ns("expression_file"), 
          label = "3. Upload expression data (CSV or text)",
          accept = c(
            "text/csv",
            "text/comma-separated-values",
            "text/tab-separated-values",
            "text/plain",
            ".csv",
            ".tsv"
          ) 
        ),
        
        # Link to public RNA-seq datasets ----------
        a(
          h4("Analyze public RNA-seq datasets for 9 species"), 
          href="http://bioinformatics.sdstate.edu/reads/"
        ),
        
        # Experiment design file input ----------
        fileInput(
          inputId = ns("experiment_file"), 
          label = h5("Optional: Upload an experiment design file(CSV or text)"),
          accept = c(
            "text/csv",
            "text/comma-separated-values",
            "text/tab-separated-values",
            "text/plain",
            ".csv",
            ".tsv"
          ) 
        ),
        
        # Table output for species loading progress -----------
        tableOutput(ns("species")),
        
        # Action button for Gene ID examples -----------
        h5(
          "Check this out if you want example of our gene ids, or download gene 
          mapping."
        ),
        actionButton(
          inputId = ns("gene_id_button"),
          label =  "Optional: Gene ID Examples"
        ),
        
        # Gene ID Examples Pop-up -----------
        shinyBS::bsModal(
          id = "gene_IDBs", 
          title = "Gene ID Examples",
          trigger = ns("gene_id_button"),
          size = "large",
          fluidPage(
            shinyjs::useShinyjs(),
            sidebarLayout(
              fluid = TRUE,
              sidebarPanel(
                
                # Select the user species ------------
                selectizeInput(
                  inputId = ns("user_specie"),
                  label = "What's your specie name?", 
                  choices = NULL
                ),
                
                shiny::tags$h5("Can erase and type in box"),
                
                # Select ID type for genes ------------                          
                selectizeInput(
                  inputId = ns("user_id_type"),
                  label = "What's your ID type? (Optional)", 
                  choices = NULL
                ),
                
                shiny::tags$h5("Can erase and type in box"),
                
                
                actionButton(inputId = ns("submit_id_page"), label = "submit"),
                actionButton(inputId = ns("reset_id_page"), label = "reset"),
                downloadButton(
                  outputId = ns("download_id_page"), 
                  label = "Download mapping.csv")
              ),
              mainPanel(
                
                reactable::reactableOutput(outputId = ns("table_result")),
                
                # User instructions ------------
                shiny::tags$div(
                  shiny::tags$h1("Instructions for Usage"),
                  shiny::tags$h4(
                    "This page's purpose is to give the user some interactive 
                     tools to look at our database IDs. There are two different 
                     uses for this page, see explanation below:"
                  ),
                  shiny::tags$ul(
                    shiny::tags$li(
                      shiny::tags$h4(
                        "If you only pick a species, you are receiving a table 
                         with all the different IDs related to that species. 
                         (Shown below)"
                      )
                    ),
                    shiny::tags$li(
                      shiny::tags$h4(
                        "If you pick a species and an ID type, you are receiving
                         a table with all the IDs of the ID type you pick and 
                         how they map to ensembl IDs(our preferred ID database), 
                         and you can download a csv file of the mapping ID."
                      )
                    )
                  )
                ), 
                
                reactable::reactableOutput(outputId = ns("table_default"))
                
              )
            )
          )
        ),
        
        a( 
          h5("Questions?",align = "right"), 
          href="https://idepsite.wordpress.com/data-format/",
          target="_blank"
        )
        
     
      ),
      
      # Load Data panel main -----------
      mainPanel(
        
        shinyjs::useShinyjs(),
        
        # Test Print -----------
        textOutput(ns("test")),
        
        # Table output for
        tableOutput('sampleInfoTable')
        
        
        
        
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
