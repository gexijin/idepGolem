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
  tabPanel(
    "Load Data",
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
          condition = 'input.select_org == "NEW"',
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
          condition = "input.data_file_format == 3",
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

        # Yes or no to converting IDs -------------
        checkboxInput(
          inputId = ns("no_id_conversion"),
          label = "Do not convert gene IDs to Ensembl.",
          value = FALSE
        ),

        # Link to public RNA-seq datasets ----------
        a(
          h4("Analyze public RNA-seq datasets for 9 species"),
          href = "http://bioinformatics.sdstate.edu/reads/"
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
        tableOutput(ns("species_match")),

        # Action button for Gene ID examples -----------
        h5(
          "Check this out if you want example of our gene ids, or download gene
          mapping."
        ),

        # ADD GENE ID EXAMPLE CODE FOR BUTTON ----------
        actionButton(
          inputId = ns("gene_id_button"),
          label = "Optional: COMING SOON!"
        ),
        ##################################################

        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/data-format/",
          target = "_blank"
        )
      ),

      # Load Data panel main -----------
      mainPanel(
        shinyjs::useShinyjs(),

        # Table output for sample tissue type ----------
        DT::dataTableOutput(ns("sample_info_table")),
        br(),

        # Display first 20 rows of the data ----------
        DT::dataTableOutput(ns("sample_20")),

        # Instructions and flowchart ------------
        div(
          id = ns("load_message"),
          h4("Loading R packages, please wait ... ... ...")
        ),
        htmlOutput(ns("file_format")),
        h3(
          "We found an issue with the Gene Onotology database derived from
           Ensembl Release 103, which is used in iDEP 0.93. While we are fixing
           this issue, we have reverted the database to a previous version used
           in iDEP 0.92. "
        ),
        h4("Postdoc and GRA positions available!"),
        h4(
          "If your gene IDs are not recognized, please let us know. We might be
           able to add customized gene mappings to Ensembl gene IDs."
        ),
        h3(
          "New version 0.93 released on 5/23/2021 includes upgrades to R 4.05,
           Bioconductor 3.12, larger database (5000+ species) from Ensembl
           Release 103 and STRING-db v11. Massive, manually-collected pathway
           database for 20 model organisms. Fixed KEGG pathway chart and gene
           plot.",
          style = "color:red"
        ),
        h4(
          "We recently hired Jenny Qi for database updates and user support.",
          a(
            "Email Jenny for questions.",
            href = "mailto:gelabinfo@gmail.com?Subject=iDEP"
          )
        ),
        h5(
          "iDEP has not been thoroughly tested. Please let us know if you find
           any issue/bug."
        ),
        h5("We will be happy to help prepare your data for iDEP."),
        br(),
        img(
          src = "www/flowchart.png",
          align = "center",
          width = "562",
          height = "383"
        ),
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

    # Provide species list for dropdown selection -----------
    observe({
      updateSelectizeInput(
        session = session,
        inputId = "select_org",
        choices = idep_data$species_choice,
        selected = idep_data$species_choice[1],
        server = TRUE
      )
    })

    output$file_format <- renderUI({
      shinyjs::hideElement(id = "load_message")
      i <- "<h3>Ready to load data files.</h3>"
      htmltools::HTML(paste(i, collapse = "<br/>"))
    })

    loaded_data <- reactive(input_data(
        expression_file = input$expression_file,
        experiment_file = input$experiment_file,
        go_button = input$go_button,
        demo_data_file = idep_data$demo_data_file,
        demo_metadata_file = idep_data$demo_metadata_file
    ))

    # Sample information table -----------
    output$sample_info_table <- DT::renderDataTable({
      req(!is.null(loaded_data()$sample_info))

      DT::datatable(
        isolate({
          tem <- t(loaded_data()$sample_info)
          tem <- cbind(rownames(tem), tem)
          colnames(tem)[1] <- "Study_design"
          tem
        }),
        options = list(
          pageLength = 20,
          scrollX = "400px",
          dom = "t",
          ordering = F
        ),
        rownames = FALSE)
    })

    # First 20 rows of dataset table -----------
    output$sample_20 <- DT::renderDataTable({
      req(!is.null(conversion_info()$converted_data))

      DT::datatable(
        conversion_info()$converted_data[1:20, ],
        options = list(
          pageLength = 20,
          scrollX = "400px",
          dom = "t"
        ),
        rownames = TRUE)
    })

    # Get converted IDs ----------
    conversion_info <- reactive({
      req(!is.null(loaded_data()$data))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Loading Data",
        color = "#000000"
      )

      converted <- convert_id(
        rownames(loaded_data()$data),
        idep_data = idep_data,
        select_org = input$select_org
      )

      gene_data <- gene_info(
        converted = converted,
        select_org = input$select_org,
        idep_data = idep_data
      )

      converted_data <- convert_data(
        converted = converted,
        no_id_conversion = input$no_id_conversion,
        data = loaded_data()$data
      )

      shinybusy::remove_modal_spinner()

      return(list(
        converted = converted,
        all_gene_info = gene_data,
        converted_data = converted_data
      ))
    })

    # Species match table ----------
    output$species_match <- renderTable({
      if (is.null(input$expression_file) && input$go_button == 0){
        return(NULL)
      }
      isolate({
        tem <- conversion_info()$converted$species_match
        if (is.null(tem)) {
          as.data.frame("ID not recognized.")
        } else {
          data.frame(
            "Matched Species (genes)" = tem[1,],
            check.names = FALSE)
        }
      })
      },
      digits = -1,
      spacing = "s",
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = TRUE
    )

    list(
      data_file_format = reactive(input$data_file_format),
      converted_data = reactive(conversion_info()$converted_data),
      sample_info = reactive(loaded_data()$sample_info),
      converted = reactive(conversion_info()$converted),
      all_gene_info = reactive(conversion_info()$all_gene_info),
      no_fdr = reactive(input$no_fdr),
      select_org = reactive(input$select_org)
    )
  })
}

## To be copied in the UI
# mod_01_load_data_ui("load_data") # nolint

## To be copied in the server
# mod_01_load_data_server("load_data") # nolint
