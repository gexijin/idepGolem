#' 01_pre_process UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_02_pre_process_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Pre-Process",
    sidebarLayout(

      # Pre-Process Panel Sidebar ----------
      sidebarPanel(

        # Conditional panel for read count data -----------
        conditionalPanel(
          condition = "output.data_file_format == 1",
          strong("Keep genes with minimal counts per million (CPM) in at
                  least n libraries:"),
          fluidRow(
            column(
              width = 6,
              numericInput(
                inputId = ns("min_counts"),
                label = h5("Min. CPM"),
                value = 0.5
              )
            ),
            column(
              width = 6,
              numericInput(
                inputId = ns("n_min_samples_count"),
                label = h5("n libraries"),
                value = 1
              )
            )
          ),
          tags$style(
            type = "text/css",
            "#pre_process-min_counts { width:100%;   margin-top:-12px}"
          ),
          tags$style(
            type = "text/css",
            "#pre_process-n_min_samples_count { width:100%;   margin-top:-12px}"
          ),
          radioButtons(
            inputId = ns("counts_transform"),
            label = "Transform counts data for clustering & PCA.",
            choices = c(
              "VST: variance stabilizing transform" = 2,
              "rlog: regularized log (slow) " = 3,
              "EdgeR: log2(CPM+c)" = 1
            ),
            selected = 1
          ),
          # Conditional panel for EdgeR transformation -----------
          conditionalPanel(
            condition = "input.counts_transform == 1",
            fluidRow(
              column(
                width = 5,
                h5("Pseudo count c:")
              ),
              column(
                width = 7,
                numericInput(
                  inputId = ns("counts_log_start"),
                  label = NULL,
                  value = 4
                )
              )
            ),
            ns = ns
          ),
          ns = ns
        ),

        # Conditional panel for FPKM data (2)----------
        conditionalPanel(
          condition = "output.data_file_format == 2",
          strong("Only keep genes above this level in at least n samples:"),
          fluidRow(
            column(
              width = 6,
              numericInput(
                inputId = ns("low_filter_fpkm"),
                label = h5("Min. level"),
                value = -1000
              )
            ),
            column(
              width = 6,
              numericInput(
                inputId = ns("n_min_samples_fpkm"),
                label = h5("n samples"),
                value = 1
              )
            )
          ),
          tags$style(
            type = "text/css",
            "#pre_process-low_filter { width:100%;   margin-top:-12px}"
          ),
          tags$style(
            type = "text/css",
            "#pre_process-n_min_samples_fpkm { width:100%;   margin-top:-12px}"
          ),
          radioButtons(
            inputId = ns("log_transform_fpkm"),
            label = "Log Transformation",
            choices = c("No" = FALSE, "Yes" = TRUE)
          ),
          numericInput(
            inputId = ns("log_start_fpkm"),
            label = h5("Constant c for started log: log(x+c)"),
            value = 1
          ),
          tags$style(
            type = "text/css",
            "#pre_process-log_start { width:100%;   margin-top:-12px}"
          ),
          textOutput(ns("text_transform")),
          tags$head(
            tags$style(
              "#pre_process-text_transform{color: blue;
               font-size: 16px;
               font-style: italic;}"
            )
          ),
          ns = ns
        ),

        # Select input for missing value ------------
        selectInput(
          inputId = ns("missing_value"),
          label = "Missing values imputation:",
          choices = list(
            "Gene median" = "geneMedian",
            "Treat as zero" = "treatAsZero",
            "Median within sample groups" = "geneMedianInGroup"
          ),
          selected = "geneMedian"
        ),

        # Button to trigger gene plot modal --------------
        actionButton(
          inputId = ns("gene_plot_1"),
          label = "Plot one or more genes"
        ),
        br(),
        br(),

        # Button to trigger examine data modal ----------
        actionButton(
          inputId = ns("examine_datab"),
          label = "Search processed data"
        ),
        br(),
        br(),

        # Yes or no to converting IDs -------------
        checkboxInput(
          inputId = ns("no_id_conversion"),
          label = "Do not convert gene IDs to Ensembl.",
          value = FALSE
        ),

        # Download button for processed data -----------
        downloadButton(
          outputId = ns("download_processed_data"),
          label = "Processed data"
        ),

        # Conditional panel for read count data ------------
        conditionalPanel(
          condition = "output.data_file_format == 1",
          downloadButton(
            outputId = ns("download_converted_counts"),
            label = "Converted counts data"
          ),
          ns = ns
        ),

        # Download EDA plot button -----------
        downloadButton(
          outputId = ns("download_eda_plot"),
          label = "High-resolution figure"
        ),
        br(),
        br(),
        textOutput(outputId = ns("n_genes_filter")),
        tags$head(tags$style(
          "#pre_process-n_genes_filter{color: blue;
            font-size: 16px;
            font-style: italic;}"
        )),
        textOutput(outputId = ns("readCountsBias")),
        tags$head(tags$style(
          "#pre_process-read_counts_bias{color: red;
            font-size: 16px;
            font-style: italic;}"
        )),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/pre-process/",
          target = "_blank"
        ),
        selectInput("counts_deg_method", "Method:",
          choices = list(
            "DESeq2" = 3,
            "limma-voom" = 2,
            "limma-trend" = 1
          ),
          selected = 3
        )
      ),

      # Pre-Process Panel Main -----------
      mainPanel(
        tableOutput(ns("sample_20")),
        h5(
          "Aspect ratios of figures can be adjusted by changing
           the width of browser window."
        ),

        # Conditional panel for barplot of read count data ----------
        conditionalPanel(
          condition = "output.data_file_format == 1",
          plotOutput(outputId = ns("total_counts")),
          ns = ns
        ),

        # Axis selectors -----------
        fluidRow(
          column(
            width = 4,
            selectInput(
              inputId = ns("scatter_x"),
              label = "Select a sample for x-axis",
              choices = 1:5,
              selected = 1
            )
          ),
          column(
            width = 4,
            selectInput(
              inputId = ns("scatter_y"),
              label = "Select a sample for y-axis",
              choices = 1:5,
              selected = 2
            )
          )
        ),

        # EDA scatter plot -----------
        plotOutput(outputId = "eda"),
        shinyBS::bsModal(
          id = "modalExample10",
          title = "Converted data (Most variable genes on top)",
          trigger = ns("examine_datab"),
          size = "large",
          DT::dataTableOutput(outputId = ns("examine_data"))
        ),
        shinyBS::bsModal(
          id = "modalExample1021",
          title = "Search for genes",
          trigger = "gene_plot_1",
          size = "large",
          textInput(
            inputId = ns("gene_search"),
            label = "Enter full or partial gene ID, or list of
                     genes separated by semicolon:",
            value = "HOXA1;e2f2;tp53"
          ),
          checkboxInput(
            inputId = ns("gene_plot_box"),
            label = "Show individual samples",
            value = FALSE
          ),
          plotOutput(outputId = "gene_plot"),
          conditionalPanel(
            condition = "input.gene_plot_box == 0",
            checkboxInput(
              inputId = ns("use_sd"),
              label = "Use standard deviation instead of standard error",
              value = FALSE
            ),
            ns = ns
          ),
          downloadButton(
            outputId = "downloadGenePlot",
            label = "Figure"
          )
        )
      )
    )
  )
}

#' 01_pre_process Server Functions
#'
#' @noRd
mod_02_pre_process_server <- function(id, load_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Data file format for conditional panels ----------
    # outputOptions required otherwise the value can only be used
    # if it is rendered somewhere else in the UI
    output$data_file_format <- reactive({
      load_data$data_file_format()
    })
    outputOptions(output, "data_file_format", suspendWhenHidden = FALSE)

    output$sample_20 <- renderTable({
      req(!is.null(load_data$data()))

      load_data$data()[1:20, ]
      },
      include.rownames = TRUE,
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = TRUE
    )

    # Update Variable Selection for the Scatter Plots ----------
    # sample_choice <- stats::setNames(
    # as.list(1:(dim(read_data()$data)[2])), colnames(read_data()$data)
    # )
    # observe({
    # updateSelectInput(
    # session,
    # inputId = "scatter_x",
    # choices = sample_choice,
    # selected = sample_choice[1]
    # )
    # })
    # observe({
    # updateSelectInput(
    # session,
    # inputId = "scatter_y",
    # choices = sample_choice,
    # selected = sample_choice[2]
    # )
    # })

    # Return Values -----------
    list(
      missing_value = reactive(input$missing_value),
      min_counts = reactive(input$min_counts),
      n_min_samples_count = reactive(input$n_min_samples_count),
      counts_log_start = reactive(input$counts_log_start),
      counts_transform = reactive(input$counts_transform),
      log_transform_fpkm = reactive(input$log_transform_fpkm),
      log_start_fpkm = reactive(input$log_start_fpkm),
      low_filter_fpkm = reactive(input$low_filter_fpkm),
      n_min_samples_fpkm = reactive(input$n_min_samples_fpkm),
      counts_deg_method = reactive(input$counts_deg_method)
    )
  })
}

## To be copied in the UI
# mod_02_pre_process_ui("pre_process") #nolint

## To be copied in the server
# mod_02_pre_process_server("pre_process_ui") #nolint
