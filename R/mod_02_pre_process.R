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
      ),

      # Pre-Process Panel Main -----------
      mainPanel(
        h5(
          "Aspect ratios of figures can be adjusted by changing
           the width of browser window."
        ),

        # Conditional panel for barplot of read count data ----------
        conditionalPanel(
          condition = "output.data_file_format == 1",
          plotOutput(
            outputId = ns("total_counts_gg"),
            width = "100%",
            height = "500px"
          ),
          ns = ns
        ),
        br(),

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
        br(),

        # EDA plots -----------
        plotOutput(
          outputId = ns("eda_scatter"),
          width = "100%",
          height = "500px"
        ),
        br(),
        plotOutput(
          outputId = ns("eda_boxplot"),
          width = "100%",
          height = "500px"
        ),
        br(),
        plotOutput(
          outputId = ns("eda_density"),
          width = "100%",
          height = "500px"
        ),
        br(),

        # Modal pop-ups ----------
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

    # Update Variable Selection for the Scatter Plots ----------
    observe({
      req(!is.null(load_data$data()))
      sample_choice <- stats::setNames(
        as.list(1:(dim(load_data$data())[2])), colnames(load_data$data())
      )
      updateSelectInput(
        session,
        inputId = "scatter_x",
        choices = colnames(load_data$data()),
        selected = colnames(load_data$data())[1]
      )
      updateSelectInput(
        session,
        inputId = "scatter_y",
        choices = colnames(load_data$data()),
        selected = colnames(load_data$data())[2]
    )
    })

    # Process the data with user defined criteria ----------
    processed_data <- reactive({
      req(!is.null(load_data$data()))

      pre_process(
        data = load_data$data(),
        missing_value = input$missing_value,
        data_file_format = load_data$data_file_format(),
        low_filter_fpkm = input$low_filter_fpkm,
        n_min_samples_fpkm = input$n_min_samples_fpkm,
        log_transform_fpkm = input$log_transform_fpkm,
        log_start_fpkm = input$log_start_fpkm,
        min_counts = input$min_counts,
        n_min_samples_count = input$n_min_samples_count,
        counts_transform = input$counts_transform,
        counts_log_start = input$counts_log_start,
        no_fdr = load_data$no_fdr()
      )
    })

    # Counts barplot ------------
    output$total_counts_gg <- renderPlot({
      req(!is.null(processed_data()$raw_counts))

      total_counts_ggplot(
        counts_data = processed_data()$raw_counts,
        sample_info = load_data$sample_info()
      )
    })

    # Scatter eda plot ----------
    output$eda_scatter <- renderPlot({
      req(!is.null(processed_data()$data))

      eda_scatter(
        processed_data = processed_data()$data,
        plot_xaxis = input$scatter_x,
        plot_yaxis = input$scatter_y
      )
    })

    # Box eda plot ----------
    output$eda_boxplot <- renderPlot({
      req(!is.null(processed_data()$data))

      eda_boxplot(
        processed_data = processed_data()$data,
        sample_info = load_data$sample_info()
      )
    })

    # Density eda plot ----------
    output$eda_density <- renderPlot({
      req(!is.null(processed_data()$data))

      eda_density(
        processed_data = processed_data()$data,
        sample_info = load_data$sample_info()
      )
    })

    # Return Values -----------
    list(
      NULL
    )
  })
}

## To be copied in the UI
# mod_02_pre_process_ui("pre_process") #nolint

## To be copied in the server
# mod_02_pre_process_server("pre_process_ui") #nolint
