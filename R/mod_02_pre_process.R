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
            "#pre_process-low_filter_fpkm { width:100%;margin-top:-12px}"
          ),
          tags$style(
            type = "text/css",
            "#pre_process-n_min_samples_fpkm { width:100%;margin-top:-12px}"
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
        br(),

        strong("Download Processed Data"),
        br(),
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
        br(),
        br(),
        textOutput(outputId = ns("n_genes_filter")),
        tags$head(tags$style(
          "#pre_process-n_genes_filter{color: blue;
            font-size: 16px;
            font-style: italic;}"
        )),
        textOutput(outputId = ns("read_counts_bias")),
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
           the width of browser window. (To save a plot, right-click)"
        ),
        tabsetPanel(
          id = ns("eda_tabs"),
          tabPanel(
            title = "Barplot",
            br(),
            plotOutput(
              outputId = ns("total_counts_gg"),
              width = "100%",
              height = "500px"
            )
          ),
          tabPanel(
            title = "Scatterplot",
            # Axis selectors -----------
            br(),
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
            plotOutput(
              outputId = ns("eda_scatter"),
              width = "100%",
              height = "500px"
            )
          ),
          tabPanel(
            title = "Boxplot",
            br(),
            plotOutput(
              outputId = ns("eda_boxplot"),
              width = "100%",
              height = "500px"
            )
          ),
          tabPanel(
            title = "Density Plot",
            br(),
            plotOutput(
              outputId = ns("eda_density"),
              width = "100%",
              height = "500px"
            ),
          ),
          tabPanel(
            title = "Converted Data",
            br(),
            DT::dataTableOutput(outputId = ns("examine_data"))
          ),
          tabPanel(
            title = "Individual Genes",
            br(),
            fluidRow(
              column(4,
                selectizeInput(
                  inputId = ns("selected_gene"),
                  label = "Select/Search for Genes",
                  choices = "",
                  selected = NULL,
                  multiple = TRUE
                )
              ),
              column(4,
                checkboxInput(
                  inputId = ns("gene_plot_box"),
                  label = "Show individual samples",
                  value = FALSE
                ),
                uiOutput(ns("sd_checkbox"))
              ),
              column(4,
                radioButtons(
                  inputId = ns("angle_ind_axis_lab"),
                  label = "Angle Axis Labels",
                  choices = c(0, 45, 90),
                  selected = 0
                )
              )
            ),
            plotOutput(
              outputId = ns("gene_plot"),
              width = "100%",
              height = "500px"
            )
          )
        )
      )
    )
  )
}

#' 01_pre_process Server Functions
#'
#' @noRd
mod_02_pre_process_server <- function(id, load_data, tab) {
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
      req(!is.null(load_data$converted_data()))

      updateSelectInput(
        session,
        inputId = "scatter_x",
        choices = colnames(load_data$converted_data()),
        selected = colnames(load_data$converted_data())[1]
      )
      updateSelectInput(
        session,
        inputId = "scatter_y",
        choices = colnames(load_data$converted_data()),
        selected = colnames(load_data$converted_data())[2]
      )
    })

    # Dynamic Barplot Tab ----------
    observe({
      if (load_data$data_file_format() != 1) {
        hideTab(inputId = "eda_plots", target = "Barplot")
        updateTabsetPanel(session, "eda_plots", selected = "Scatterplot")
      } else if (load_data$data_file_format() == 1) {
        showTab(inputId = "eda_plots", target = "Barplot")
        updateTabsetPanel(session, "eda_plots", selected = "Barplot")
      }
    })

    # Process the data with user defined criteria ----------
    processed_data <- reactive({
      req(!is.null(load_data$converted_data()))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Pre-Processing Data",
        color = "#000000"
      )

      processed_data <- pre_process(
        data = load_data$converted_data(),
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

      shinybusy::remove_modal_spinner()

      return(processed_data)
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

    # Merge Data Sets with Gene names ----------
    merged_processed_data <- reactive({
      req(!is.null(processed_data()$data))

      merged_data <- merge_data(
        load_data$all_gene_names(),
        processed_data()$data
      )
    })
    merged_raw_counts_data <- reactive({
      req(!is.null(processed_data()$data))

      merged_data <- merge_data(
        load_data$all_gene_names(),
        processed_data()$raw_counts
      )
    })

    # Pre-Process Data Table ----------
    output$examine_data <- DT::renderDataTable({
      req(!is.null(merged_processed_data()))

      DT::datatable(
        merged_processed_data(),
        options = list(
          pageLength = 10,
          scrollX = "400px"
        ),
        rownames = FALSE)
    })

    # Individual genes selection ----------
    observe({
      req(tab() == "Pre-Process")

      if (is.null(merged_processed_data()$symbol)) {
        choices <- merged_processed_data()$User_id
      } else {
        choices <- merged_processed_data()$symbol
      }
      updateSelectizeInput(
        session,
        inputId = "selected_gene",
        choices = choices,
        selected = NULL,
        server = TRUE
      )
    })

    # Dynamic individual gene checkbox ----------
    output$sd_checkbox <- renderUI({
      req(input$gene_plot_box == FALSE)

      checkboxInput(
        inputId = ns("use_sd"),
        label = "Use standard deviation instead of standard error",
        value = FALSE
      )
    })

    output$gene_plot <- renderPlot({
      req(!is.null(merged_processed_data()))
      req(!is.null(input$selected_gene))

      individual_plots(
        merged_data = merged_processed_data(),
        sample_info = load_data$sample_info(),
        selected_gene = input$selected_gene,
        gene_plot_box = input$gene_plot_box,
        use_sd = input$use_sd,
        lab_rotate = input$angle_ind_axis_lab
      )
    })

    # Download buttons ----------
    output$download_processed_data <- downloadHandler(
      filename = function() {
        "processed_data.csv"
      },
      content = function(file) {
        write.csv(merged_processed_data(), file)
      }
    )
    output$download_converted_counts <- downloadHandler(
      filename = function() {
        "converted_counts_data.csv"
      },
      content = function(file) {
        write.csv(merged_raw_counts_data(), file)
      }
    )

    # Text Output Information -----------
    output$n_genes_filter <- renderText({
      req(processed_data()$data_size)

      if (dim(load_data$all_gene_names())[2] == 1) {
        return(paste(
          processed_data()$data_size[1], "genes in",
          processed_data()$data_size[4], "samples.",
          processed_data()$data_size[3], " genes passed filter
          (see above). Original gene IDs used."
        ))
      } else {
        return(paste(
          processed_data()$data_size[1], "genes in",
          processed_data()$data_size[4], "samples.",
          processed_data()$data_size[3], " genes passed filter
          (see above), ", load_data$n_matched(),
          " were converted to Ensembl gene IDs in our database.
          The remaining ", processed_data()$data_size[3] -
          load_data$n_matched(), " genes were kept in the
          data using original IDs."
        ))
      }
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
