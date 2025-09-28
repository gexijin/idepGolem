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
    title = "Prep",
    sidebarLayout(

      # Pre-Process Panel Sidebar ----------
      sidebarPanel(

        # Conditional panel for read count data -----------
        conditionalPanel(
          condition = "output.data_file_format == 1",
          strong("1. Filter genes with low counts"),
          fluidRow(
            column(
              width = 6,
              numericInput(
                inputId = ns("min_counts"),
                label = tags$span("Min. CPM", style = "font-weight: normal;"),
                value = 0.5
              ),
              tippy::tippy_this(
                ns("min_counts"),
                "Counts Per Million (CPM) = (read count / total counts) * 1,000,000",
                theme = "light-border"
              )
            ),
            column(
              width = 6,

              # Min samples per row to have min CPM
              numericInput(
                inputId = ns("n_min_samples_count"),
                label = tags$span("n libraries", style = "font-weight: normal;"),
                value = 1
              ),
              tippy::tippy_this(
                ns("n_min_samples_count"),
                "Number of samples (libraries) that must have at least the min CPM",
                theme = "light-border"
              )
            )
          ),
          strong("2. Transform counts data"),
          selectInput(
            inputId = ns("counts_transform"),
            label = NULL,
            choices = c(
              "VST: variance stabilizing transform" = 2,
              "rlog: regularized log (slow) " = 3,
              "EdgeR: log2(CPM+c)" = 1
            ),
            selected = 1
          ),
          tippy::tippy_this(
            ns("counts_transform"),
            "Transformed data is used in all analyses except differential expression with DESeq2.",
            theme = "light-border"
          ),

          # Conditional panel for EdgeR transformation -----------
          conditionalPanel(
            condition = "input.counts_transform == 1",
            fluidRow(
              column(
                width = 2,
              ),              
              column(
                width = 7,
                "Pseudo count c:"
              ),
              column(
                width = 3,

                # Constant to add for a log transform
                numericInput(
                  inputId = ns("counts_log_start"),
                  label = NULL,
                  value = 4
                ),
                tippy::tippy_this(
                  ns("counts_log_start"),
                  "Constant c for log2(CPM+c). A larger c shrinks log-values towards log2(c).",
                  theme = "light-border"
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
          strong("1. Filter genes with low expression:"),
          fluidRow(
            column(
              width = 6,
              # Fold counts min (works with min samples)
              numericInput(
                inputId = ns("low_filter_fpkm"),
                label = tags$span("Min. level", style = "font-weight: normal;"),
                value = -1000
              ),
              tippy::tippy_this(
                ns("low_filter_fpkm"),
                "Minimum expression level, e.g. FPKM, RPKM, TPM, or other normalized values.",
                theme = "light-border"
              )
            ),
            column(
              width = 6,
              # Min samples per row to have the low filter
              numericInput(
                inputId = ns("n_min_samples_fpkm"),
                label = tags$span("n samples", style = "font-weight: normal;"),
                value = 1
              ),
              tippy::tippy_this(
                ns("n_min_samples_fpkm"),
                "Number of samples that must have at least the min expression level.",
                theme = "light-border"
              )
            )
          ),
          tags$style(
            type = "text/css",
            "#pre_process-low_filter_fpkm { width:100%;margin-top:-5px}"
          ),
          tags$style(
            type = "text/css",
            "#pre_process-n_min_samples_fpkm { width:100%;margin-top:-5px}"
          ),

          # Perform a log transform or not


          fluidRow(
            column(
              width = 2,
              strong("2.")
            ),
            column(
              width = 10,
              align = "left",
              style = "margin-top: -10px;",
              checkboxInput(
                inputId = ns("log_transform_fpkm"),
                label = strong("Log Transformation"),
                value = FALSE
              ),
              tippy::tippy_this(
                ns("log_transform_fpkm"),
                "Log transformed data is used in all analyses.",
                theme = "light-border"
              )
            )
          ),

          conditionalPanel(
            condition = "input.log_transform_fpkm",
            fluidRow(
              column(
                width = 1,
              ),
              column(
                width = 8,
                "constant c in log2(x+c):"
              ),
              column(
                width = 3,
                # Constant to add for a log transform
                numericInput(
                  inputId = ns("log_start_fpkm"),
                  label = NULL,
                  value = 1
                ),
                tags$style(
                  type = "text/css",
                  "#pre_process-log_start_fpkm { width:100%;   margin-top:-5px}"
                )
              )
            ),
            ns = ns
          ),
          ns = ns
        ),

        # Select input for missing value ------------
        fluidRow(
          column(
            width = 5,
            strong("3. Missing values:")
          ),
          column(
            width = 7,

            # Constant to add for a log transform
            selectInput(
              inputId = ns("missing_value"),
              label = NULL,
              choices = list(
                "Use gene median" = "geneMedian",
                "Treat as zero" = "treatAsZero",
                "Use gene median in group" = "geneMedianInGroup"
              ),
              selected = "geneMedian"
            ),
            tippy::tippy_this(
              ns("missing_value"),
              "How to handle missing values in the data matrix.",
              theme = "light-border"
            )
          )
        ),
        fluidRow(
          column(
            width = 6,
            downloadButton(
              outputId = ns("download_processed_data"),
              label = "Processed data"
            ),
            tippy::tippy_this(
              ns("download_processed_data"),
              "Download transformed data",
              theme = "light-border"
            )
          ),
          column(
            width = 6,
            # Conditional panel for read count data ------------
            conditionalPanel(
              condition = "output.data_file_format == 1",

              # Download the counts data with converted IDs
              downloadButton(
                outputId = ns("download_converted_counts"),
                label = "Converted counts"
              ),
              tippy::tippy_this(
                ns("download_converted_counts"),
                "Download counts data with converted IDs",
                theme = "light-border"
              ),
              ns = ns
            )
          )
        ),
        br(),
        downloadButton(
          outputId = ns("rds"),
          label = ".RData"
        ),
        tippy::tippy_this(
          ns("rds"),
          "Download converted data as .Rdata format",
          theme = "light-border"
        ),
        downloadButton(
          outputId = ns("report"),
          label = "Report"
        ),
        tippy::tippy_this(
          ns("report"),
          "Generate HTML report of pre-processing tab",
          theme = "light-border"
        ),
        uiOutput(ns("mapping_statistics_container")),
      ),


      # Pre-Process Panel Main -----------
      mainPanel(
        tabsetPanel(
          id = ns("eda_tabs"),

          # Barplot for read counts data ----------
          tabPanel(
            title = "Barplot",
            br(),
            plotOutput(
              outputId = ns("raw_counts_gg"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_raw_counts_gg")
            ),
            br(),
            tableOutput(
              outputId = ns("counts_table")
            )
          ),

          # Boxplot of transformed data ----------
          tabPanel(
            title = "Boxplot",
            br(),
            plotOutput(
              outputId = ns("eda_boxplot"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_eda_boxplot")
            )
          ),

          # Density plot of transformed data ---------
          tabPanel(
            title = "Density Plot",
            br(),
            plotOutput(
              outputId = ns("eda_density"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_eda_density")
            )
          ),

          # Scatterplot with interactive axes ----------
          tabPanel(
            title = "Scatterplot",
            # Axis selectors -----------
            br(),
            fluidRow(
              column(
                width = 4,
                selectInput(
                  inputId = ns("scatter_x"),
                  label = tags$span("Sample for x-axis", style = "font-weight: normal;"),
                  choices = 1:5,
                  selected = 1
                )
              ),
              column(
                width = 4,
                selectInput(
                  inputId = ns("scatter_y"),
                  label = tags$span("Sample for y-axis", style = "font-weight: normal;"),
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
            ),
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_eda_scatter")
            )
          ),

          # Density plot of transformed data ---------
          tabPanel(
            title = "Dispersion",
            br(),
            fluidRow(
              column(
                width = 4,
                selectInput(
                  inputId = ns("heat_color_select"),
                  label = "Select Heat Colors",
                  choices = NULL
                )
              ),
              column(
                width = 4,
                checkboxInput(
                  inputId = ns("rank"),
                  label = "Use rank of mean values"
                )
              ),
            ),
            plotOutput(
              outputId = ns("dev_transfrom"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_dev_transform")
            )
          ),

          # Barplot for rRNA counts ----------
          tabPanel(
            title = "QC",
            conditionalPanel(
              condition = "output.select_org == 'NEW'",
              h4("QC reports are not available for custom organisms."),
              ns = ns
            ),
            conditionalPanel(
              condition = "output.select_org != 'NEW'",

              plotOutput(
                outputId = ns("gene_counts_gg"),
                width = "100%",
                height = "500px"
              ),
              br(),
              fluidRow(
                column(
                  2, 
                  ottoPlots::mod_download_figure_ui(
                    id = ns("dl_gene_counts_gg")
                  )
                ),
                column(
                  10,
                  align = "right",
                  p(
                    "Genes types are based on ",
                    a(
                      "ENSEMBL classification.", 
                      href= "http://useast.ensembl.org/info/genome/genebuild/biotypes.html"),
                      target = "_blank"
                  )
                )
              ),
              br(),
              hr(),
            
              conditionalPanel(
                condition = "output.data_file_format == 1",
                br(),
                plotOutput(
                  outputId = ns("rRNA_counts_gg"),
                  width = "100%",
                  height = "500px"
                ),
                br(),
                fluidRow(
                  column(
                    2, 
                    ottoPlots::mod_download_figure_ui(
                      id = ns("dl_rRNA_counts_gg")
                    )
                  ),
                  column(
                    10,
                    align = "right",
                    p("Higher rRNA content may indicuate ineffective rRNA-removal.")
                  )
                ),
                br(),
                hr(),
                plotOutput(
                  outputId = ns("chr_counts_gg"),
                  width = "100%",
                  height = "2000px"
                ),
                br(),
                fluidRow(
                  column(
                    2, 
                    ottoPlots::mod_download_figure_ui(
                      id = ns("dl_chr_counts_gg")
                    ),
                  ),
                  column(
                    10,
                    align = "right",
                    p("Higher bar means more reads map to this chromosome in this sample.")
                  )
                ),
                br(),
                hr(),
                ns = ns
              ),
            
              plotOutput(
                outputId = ns("chr_normalized_gg"),
                width = "100%",
                height = "2000px"
              ),
              br(),
              fluidRow(
                column(
                  2, 
                  ottoPlots::mod_download_figure_ui(
                    id = ns("dl_chr_normalized_gg")
                  )
                ),
                column(
                  10,
                  align = "right",
                  p("A tall bar means genes on this chromosome are expressed at higher levels in a sample.")
                )
              ),
            ns = ns
            )
          ),
          # Plot panel for individual genes ---------
          tabPanel(
            title = "Gene plot",
            br(),
            fluidRow(
              column(
                4,
                # Gene ID Selection -----------
                selectizeInput(
                  inputId = ns("selected_gene"),
                  label = "Select/Search for Gene(s)",
                  choices = "",
                  selected = NULL,
                  multiple = TRUE
                )
              ),
              column(
                4,
                selectInput(
                  inputId = ns("gene_plot_box"),
                  label = "Plot Type:",
                  choices = setNames(
                    c(1,2), 
                    c("Sample Group Expression",
                      "Individual Sample Expression")
                  ),
                  selected = 1
                ),
                uiOutput(ns("sd_checkbox")),
                conditionalPanel(
                  condition = "output.data_file_format == 1",
                  checkboxInput(
                    inputId = ns("plot_raw"),
                    label = "Plot raw counts",
                    value = FALSE
                  ),
                  ns = ns
                ),
                checkboxInput(
                  inputId = ns("plot_tukey"),
                  label = "Run TukeyHSD test",
                  value = FALSE
                )
              ),
              column(
                4,
                radioButtons(
                  inputId = ns("angle_ind_axis_lab"),
                  label = "Angle Axis Labels",
                  choices = c(0, 45, 90),
                  selected = 45
                )
              )
            ),

            plotOutput(
              outputId = ns("gene_plot"),
              width = "100%",
              height = "500px"
            ),
            div(
              style = "display: flex; gap: 10px",
              ottoPlots::mod_download_figure_ui(
                id = ns("dl_gene_plot")
              ),
              downloadButton(
                outputId = ns("tukey_download"),
                label = "TukeyHSD Results"
              ),
              uiOutput(ns("signif_text"))
            )
          ),


          # Searchable table of transformed converted data ---------
          tabPanel(
            title = "Data",
            br(),
            conditionalPanel(
              condition = "output.data_file_format == 1",
              checkboxInput(
                inputId = ns("show_raw"),
                label = "Show raw counts, not transformed data",
                value = FALSE
              ),
              ns = ns
            ),
            br(),
            DT::dataTableOutput(outputId = ns("examine_data"))
          ),
          tabPanel(
            title = icon("info-circle"),
            includeHTML(app_sys("app/www/help_preprocess.html"))
          ),
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
    density_tip_shown <- reactiveVal(FALSE)
    
    # Data file format for conditional panels ----------
    # outputOptions required otherwise the value can only be used
    # if it is rendered somewhere else in the UI
    output$data_file_format <- reactive({
      load_data$data_file_format()
    })
    outputOptions(output, "data_file_format", suspendWhenHidden = FALSE)

    output$select_org <- reactive({
      load_data$select_org()
    })
    outputOptions(output, "select_org", suspendWhenHidden = FALSE)

    # Update Variable Selection for the Scatter Plots ----------
    observe({
      req(!is.null(load_data$converted_data()))

      updateSelectInput(
        session,
        inputId = "scatter_x",
        choices = colnames(processed_data()$data),
        selected = colnames(processed_data()$data)[1]
        # load_data$converted_data())[1]
      )
      updateSelectInput(
        session,
        inputId = "scatter_y",
        choices = colnames(processed_data()$data),
        selected = colnames(processed_data()$data)[2]
      )
    })

    # Dynamic Barplot Tab ----------
    observe({
      if (load_data$data_file_format() != 1) {
        hideTab(inputId = "eda_tabs", target = "Barplot")
        updateTabsetPanel(session, "eda_tabs", selected = "Boxplot")
      } else if (load_data$data_file_format() == 1) {
        showTab(inputId = "eda_tabs", target = "Barplot")
        updateTabsetPanel(session, "eda_tabs", selected = "Barplot")
      }
    })

    observeEvent(input$eda_tabs, {
      req(input$eda_tabs == "Density Plot")
      if (!density_tip_shown()) {
        showNotification(
          "Figure width can be adjusted by changing the width of browser window.",
          id = "boxplot_width_tip",
          duration = 15,
          type = "message"
        )
        density_tip_shown(TRUE)
      }
    }, ignoreNULL = TRUE)

    # Process the data with user defined criteria ----------
    processed_data <- reactive({
      req(!is.null(load_data$converted_data()))
      req(input$n_min_samples_count)
      req(input$min_counts)
      req(input$low_filter_fpkm)
      req(input$n_min_samples_fpkm)
      req(input$counts_log_start)
      req(input$log_start_fpkm)

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
    
    observe({
      req(processed_data()$data_size[3] < 1000)
      showNotification(
        "Less than 1000 genes passed through the pre-processing filter. 
         By default, all genes in the database from the selected species 
         will be used as background in enrichment analysis.",
        id = "filter_warning",
        duration = 20,
        type = "warning",
        session = session
      )
    })
    
    observe({
      req(!tab() %in% c("Clustering", "Load Data", 
                        "Network", "Bicluster", 
                        "DEG2"))
      removeNotification(id = "filter_warning",
                         session = session)
    })

    # Counts barplot ------------
    raw_counts <- reactive({
      req(!is.null(processed_data()$raw_counts))
      withProgress(message = "Generating plot ...", {
        incProgress(0.2)
        p <- total_counts_ggplot(
          counts_data = processed_data()$raw_counts,
          sample_info = load_data$sample_info(),
          type = "Raw",
          plots_color_select = load_data$plots_color_select()
        )
        refine_ggplot2(
          p = p,
          gridline = load_data$plot_grid_lines(),
          ggplot2_theme = load_data$ggplot2_theme()
        )
      })
    })
    
    output$raw_counts_gg <- renderPlot({
      print(raw_counts())
    })
    
    dl_raw_counts_gg <- ottoPlots::mod_download_figure_server(
      id = "dl_raw_counts_gg",
      filename = "raw_counts_barplot",
      figure = reactive({
        raw_counts()
      }),
      label = ""
    )
    
    output$counts_table <- renderTable({
      req(!is.null(processed_data()))
      # Sums of counts by sample
      filtered <- as.data.frame(colSums(processed_data()$raw_counts))
      original <- as.data.frame(colSums(load_data$converted_data()))
      
      # Merge and rename columns
      df <- merge(x = original, y = filtered, by = "row.names")
      colnames(df) <- c("Sample", "Total Counts", "Total Counts After Filter")
      df$`Counts Removed` <- df[,2] - df[,3] # Gene counts lost
      df <- dplyr::mutate(df, across(2:4, ~ format(.x, big.mark = ",")))
      df
    },
    bordered = TRUE,
    digits = 0,
    hover = TRUE
    )

    # gene type barplot ------------
    gene_counts <- reactive({
      req(!is.null(processed_data()$raw_counts))
      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Plotting counts by gene type",
        color = "#000000"
      )
      p <- gene_counts_ggplot(
        counts_data = load_data$converted_data(),
        sample_info = load_data$sample_info(),
        type = "Raw",
        all_gene_info = load_data$all_gene_info(),
        plots_color_select = load_data$plots_color_select()
      )
      p <- refine_ggplot2(
        p = p,
        gridline = load_data$plot_grid_lines(),
        ggplot2_theme = load_data$ggplot2_theme()
      )
      shinybusy::remove_modal_spinner()
      return(p)
    })

    output$gene_counts_gg <- renderPlot({
      print(gene_counts())
    })

    dl_gene_counts_gg <- ottoPlots::mod_download_figure_server(
      id = "dl_gene_counts_gg",
      filename = "gene_counts_barplot",
      figure = reactive({
        gene_counts()
      }),
      label = ""
    )

    # gene type barplot ------------
    rRNA_counts <- reactive({
      req(!is.null(processed_data()$raw_counts))
      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Plotting counts by gene type",
        color = "#000000"
      )
      p <- rRNA_counts_ggplot(
        counts_data = load_data$converted_data(),
        sample_info = load_data$sample_info(),
        type = "Raw",
        all_gene_info = load_data$all_gene_info(),
        plots_color_select = load_data$plots_color_select()
      )
      p <- refine_ggplot2(
        p = p,
        gridline = load_data$plot_grid_lines(),
        ggplot2_theme = load_data$ggplot2_theme()
      )
      shinybusy::remove_modal_spinner()
      return(p)
    })
    output$rRNA_counts_gg <- renderPlot({
      print(rRNA_counts())
    })

    dl_rRNA_counts_gg <- ottoPlots::mod_download_figure_server(
      id = "dl_rRNA_counts_gg",
      filename = "rRNA_counts_barplot",
      figure = reactive({
        rRNA_counts()
      }),
      label = ""
    )

    # chr counts barplot ------------
    chr_counts <- reactive({
      req(!is.null(processed_data()$raw_counts))
      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Plotting counts by Chromosome",
        color = "#000000"
      )
      p <- chr_counts_ggplot(
        counts_data = load_data$converted_data(),
        sample_info = load_data$sample_info(),
        type = "Raw",
        all_gene_info = load_data$all_gene_info()
      )
      p <- refine_ggplot2(
        p = p,
        gridline = load_data$plot_grid_lines(),
        ggplot2_theme = load_data$ggplot2_theme()
      )
      shinybusy::remove_modal_spinner()
      return(p)
    })
    output$chr_counts_gg <- renderPlot({
      print(chr_counts())
    },
    height = 2000)

    dl_chr_counts_gg <- ottoPlots::mod_download_figure_server(
      id = "dl_chr_counts_gg",
      filename = "Chr_counts_barplot",
      figure = reactive({
        chr_counts()
      }),
      label = ""
    )


    # chr normalized barplot ------------
    chr_normalized <- reactive({
      req(!is.null(processed_data()$data))
      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Pre-Processing Data",
        color = "#000000"
      )
      p <- chr_normalized_ggplot(
        counts_data = processed_data()$data,
        sample_info = load_data$sample_info(),
        type = "Raw",
        all_gene_info = load_data$all_gene_info()
      )
      p <- refine_ggplot2(
        p = p,
        gridline = load_data$plot_grid_lines(),
        ggplot2_theme = load_data$ggplot2_theme()
      )
      shinybusy::remove_modal_spinner()
      return(p)
    })
    output$chr_normalized_gg <- renderPlot({
      print(chr_normalized())
    },
    height = 2000)

    dl_chr_counts_gg <- ottoPlots::mod_download_figure_server(
      id = "dl_chr_normalized_gg",
      filename = "Chr_normalized_expression_barplot",
      figure = reactive({
        chr_normalized()
      }),
      label = ""
    )


    # Scatter eda plot ----------
    scatter <- reactive({
      req(!is.null(processed_data()$data))
      withProgress(message = "Generating plot ...", {
        incProgress(0.2)
        p <- eda_scatter(
          processed_data = processed_data()$data,
          plot_xaxis = input$scatter_x,
          plot_yaxis = input$scatter_y
        )
        refine_ggplot2(
          p = p,
          gridline = load_data$plot_grid_lines(),
          ggplot2_theme = load_data$ggplot2_theme()
        )
      })
    })
    output$eda_scatter <- renderPlot({
      print(scatter())
    })
    dl_eda_scatter <- ottoPlots::mod_download_figure_server(
      id = "dl_eda_scatter",
      filename = "scatter_plot",
      figure = reactive({
        scatter()
      }),
      label = ""
    )

    # Box eda plot ----------
    eda_box <- reactive({
      req(!is.null(processed_data()$data))
      withProgress(message = "Generating boxplot", {
        incProgress(0.2)
        p <- eda_boxplot(
          processed_data = processed_data()$data,
          sample_info = load_data$sample_info(),
          plots_color_select = load_data$plots_color_select()
        )
        refine_ggplot2(
          p = p,
          gridline = load_data$plot_grid_lines(),
          ggplot2_theme = load_data$ggplot2_theme()
        )
      })
    })
    output$eda_boxplot <- renderPlot({
      print(eda_box())
    })
    dl_eda_boxplot <- ottoPlots::mod_download_figure_server(
      id = "dl_eda_boxplot",
      filename = "transformed_boxplot",
      figure = reactive({
        eda_box()
      }),
      label = ""
    )

    # Density eda plot ----------
    density <- reactive({
      req(!is.null(processed_data()$data))
      withProgress(message = "Generating density plot", {
        incProgress(0.2)
        p <- eda_density(
          processed_data = processed_data()$data,
          sample_info = load_data$sample_info(),
          plots_color_select = load_data$plots_color_select()
        )
        refine_ggplot2(
          p = p,
          gridline = load_data$plot_grid_lines(),
          ggplot2_theme = load_data$ggplot2_theme()
        )
      })
    })
    output$eda_density <- renderPlot({
      print(density())
    })
    dl_eda_density <- ottoPlots::mod_download_figure_server(
      id = "dl_eda_density",
      filename = "density_plot",
      figure = reactive({
        density()
      }),
      label = ""
    )

    # Standard deviation vs mean plot ----------
    # Heatmap Colors ----------
    heat_colors <- list(
      "Green" = c("green"),
      "Red" = c("red"),
      "Magenta" = c("magenta"),
      "Blue" = c("blue"),
      "Brown" = c("brown")
    )
    heat_choices <- c(
      "Green",
      "Red",
      "Magenta",
      "Blue",
      "Brown"
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "heat_color_select",
        choices = heat_choices
      )
    })

    # Mean vs SD plot --------
    dev <- reactive({
      req(!is.null(processed_data()$data))
      withProgress(message = "Generating plot ...", {
        incProgress(0.2)
        p <- mean_sd_plot(
          processed_data = processed_data()$data,
          heat_cols = heat_colors[[input$heat_color_select]],
          rank = input$rank
        )
        refine_ggplot2(
          p = p,
          gridline = load_data$plot_grid_lines(),
          ggplot2_theme = load_data$ggplot2_theme()
        )
      })
    })
    output$dev_transfrom <- renderPlot({
      print(dev())
    })
    dl_dev_transform <- ottoPlots::mod_download_figure_server(
      id = "dl_dev_transform",
      filename = "transform_plot",
      figure = reactive({
        dev()
      }),
      label = ""
    )

    # Merge Data Sets with Gene names ----------
    merged_processed_data <- reactive({
      req(!is.null(processed_data()$data))

      merged_data <- merge_data(
        load_data$all_gene_names(),
        processed_data()$data,
        merge_ID = "ensembl_ID"
      )
    })
    merged_raw_counts_data <- reactive({
      req(!is.null(processed_data()$data))

      merged_data <- merge_data(
        load_data$all_gene_names(),
        processed_data()$raw_counts,
        merge_ID = "ensembl_ID"
      )
    })

    # Pre-Process Data Table ----------
    output$examine_data <- DT::renderDataTable({
      req(!is.null(merged_processed_data()))

      if(input$show_raw) {
        data_matrix <- merged_raw_counts_data()
      } else {
        data_matrix <- merged_processed_data()
      }

      DT::datatable(
        data_matrix,
        options = list(
          pageLength = 20,
          scrollX = "400px"
        ),
        rownames = FALSE
      )
    })

    # Individual plot data ------------
    individual_data <- reactive({
      req(!is.null(processed_data()$data))

      if(input$plot_raw) {
        data_matrix <- processed_data()$raw_counts
      } else {
        data_matrix <- processed_data()$data
      }


      rowname_id_swap(
        data_matrix = data_matrix,
        all_gene_names = load_data$all_gene_names(),
        select_gene_id = load_data$select_gene_id()
      )
    })

    # Individual genes selection ----------
    observe({
      req(!is.null(processed_data()$data))
      #the orders of genes stays the same when user clicks on "Plot raw counts"
      isolate({
        # Genes are sorted by SD
        sorted <- sort(
          apply( # gene SD
            individual_data(),
            MARGIN = 1,
            FUN = function(x) sd(x) #/ abs(mean(x) + 1e-10) # add small number to avoid 0
          ),
          decreasing = TRUE
        )
        # top 2 most variable genes are plotted by default
        selected <- names(sorted)[1:2]

        updateSelectizeInput(
          session,
          inputId = "selected_gene", # genes are ranked by SD
          choices = names(sorted),
          selected = selected,
          server = TRUE
        )
      })
    })

    # Dynamic individual gene checkbox ----------
    output$sd_checkbox <- renderUI({
      req(input$gene_plot_box != 2)

      checkboxInput(
        inputId = ns("use_sd"),
        label = "Use standard deviation",
        value = FALSE
      )
    })

    observe({
      shinyjs::toggle(id = "plot_tukey", 
                      condition = input$gene_plot_box != 2 &&
                        !input$plot_raw)
      shinyjs::toggle(id = "tukey_download", condition = input$plot_tukey)
      
      if (input$plot_raw == TRUE && input$plot_tukey == TRUE){
        shinyjs::reset(id = "plot_tukey")
      }
    })
    
    tukey_data <- reactive({
      req(!is.null(individual_data()))
      req(!is.na(input$selected_gene))
      
      get_tukey_data(individual_data(),
                     sample_info = load_data$sample_info(),
                     selected_gene = input$selected_gene,
                     summarized = FALSE)
    })
    
    output$tukey_download <- downloadHandler(
      filename = function(){
        req(!is.null(tukey_data()))
        req(!is.na(input$selected_gene))
        
        "TukeyHSD_results.csv"
      },
      content = function(file){
        req(!is.null(tukey_data()))
        req(!is.na(input$selected_gene))
        
        write.csv(tukey_data(), file)
      }
    )
    
    # Individual gene plot ---------
    gene_plot <- reactive({
      req(individual_data())
      req(input$selected_gene)
      req(!is.null(input$gene_plot_box))
      req(!is.null(input$use_sd))
      req(!is.null(input$plot_tukey))
      req(input$angle_ind_axis_lab)
      req(input$plot_raw != TRUE || input$plot_tukey != TRUE)
      withProgress(message = "Generating plot ...", {
        incProgress(0.2)
        p <- individual_plots(
          individual_data = individual_data(),
          sample_info = load_data$sample_info(),
          selected_gene = input$selected_gene,
          gene_plot_box = input$gene_plot_box,
          use_sd = input$use_sd,
          lab_rotate = input$angle_ind_axis_lab,
          plots_color_select = load_data$plots_color_select(),
          plot_raw = input$plot_raw,
          plot_tukey = input$plot_tukey
        )
        refine_ggplot2(
          p = p,
          gridline = load_data$plot_grid_lines(),
          ggplot2_theme = load_data$ggplot2_theme()
        )
      })
    })
  
    output$gene_plot <- renderPlot({
      req(gene_plot())
      print(gene_plot())
    })

    dl_gene_plot <- ottoPlots::mod_download_figure_server(
      id = "dl_gene_plot",
      filename = "gene_plot",
      figure = reactive({
        gene_plot()
      }),
      label = ""
    )

    output$signif_text <- renderUI({
      req(!is.null(input$plot_tukey))
      
      if (input$plot_tukey == TRUE){
        paste0('Only top 10 most significant differences displayed for each',
               ' gene (*** = pval < 0.001; ** = pval < 0.01; ',
               '* = pval < 0.05)')
      } else {NULL}
    })

    # Download buttons ----------
    output$download_processed_data <- downloadHandler(
      filename = function() {
        "processed_data.csv"
      },
      content = function(file) {
        write.csv(merged_processed_data(), file, row.names = FALSE)
      }
    )
    output$download_converted_counts <- downloadHandler(
      filename = function() {
        "converted_counts_data.csv"
      },
      content = function(file) {
        write.csv(merged_raw_counts_data(), file, row.names = FALSE)
      }
    )

    # Markdown report
    output$report <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "pre_process_report.html",
      content = function(file) {
        withProgress(message = "Generating Report (5 mins)", {
          incProgress(0.2)
          # Copy the report file to a temporary directory before processing it, in
          # case we don't have write permissions to the current working dir (which
          # can happen when deployed).
          tempReport <- file.path(tempdir(), "pre_process_workflow.Rmd")
          # tempReport
          tempReport <- gsub("\\", "/", tempReport, fixed = TRUE)

          # This should retrieve the project location on your device:
          # "C:/Users/bdere/Documents/GitHub/idepGolem"
          wd <- getwd()
          
          markdown_location <- app_sys("app/www/RMD/pre_process_workflow.Rmd")
          file.copy(from = markdown_location, to = tempReport, overwrite = TRUE)

          # Set up parameters to pass to Rmd document
          params <- list(
            loaded_data = load_data$converted_data(),
            individual_data = individual_data(),
            descr = processed_data()$descr,
            sample_info = load_data$sample_info(),
            all_gene_info = load_data$all_gene_info(),
            data_file_format = load_data$data_file_format(),
            no_id_conversion = input$no_id_conversion,
            min_counts = input$min_counts,
            n_min_samples_count = input$n_min_samples_count,
            counts_transform = input$counts_transform,
            counts_log_start = input$counts_log_start,
            log_transform_fpkm = input$log_transform_fpkm,
            log_start_fpkm = input$log_start_fpkm,
            low_filter_fpkm = input$low_filter_fpkm,
            missing_value = input$missing_value,
            scatter_x = input$scatter_x,
            scatter_y = input$scatter_y,
            sd_color = heat_colors[[input$heat_color_select]],
            rank = input$rank,
            no_fdr = load_data$no_fdr(),
            selected_gene = input$selected_gene,
            gene_plot_box = input$gene_plot_box,
            use_sd = input$use_sd,
            lab_rotate = input$angle_ind_axis_lab,
            plots_color_select = load_data$plots_color_select()
          )
          req(params)

          # Knit the document, passing in the `params` list, and eval it in a
          # child of the global environment (this isolates the code in the document
          # from the code in this app).
          rmarkdown::render(
            input = tempReport, # markdown_location,
            output_file = file,
            params = params,
            envir = new.env(parent = globalenv())
          )
        })
      }
    )
    # RDS with data and inputs
    output$rds <- downloadHandler(
      filename = paste0("idep_session_", format(Sys.time(), "%Y_%m_%d"), ".Rdata"),
      content = function(file) {
        if (load_data$data_file_format() == 1) {
          loaded_data <- load_data$converted_data()
          sample_info <- load_data$sample_info()
          data_file_format <- load_data$data_file_format()
          no_id_conversion <- input$no_id_conversion
          min_counts <- input$min_counts
          n_min_samples_count <- input$n_min_samples_count
          counts_transform <- input$counts_transform
          counts_log_start <- input$counts_log_start
          missing_value <- input$missing_value
          scatter_x <- input$scatter_x
          scatter_y <- input$scatter_y
          sd_color <- heat_colors[[input$heat_color_select]]
          rank <- input$rank

          save(loaded_data,
            sample_info,
            data_file_format,
            no_id_conversion,
            min_counts,
            n_min_samples_count,
            counts_transform,
            counts_log_start,
            missing_value,
            scatter_x,
            scatter_y,
            sd_color,
            rank,
            file = file
          )
        }
        if (load_data$data_file_format() == 2) {
          params_r <- list(
            loaded_data = load_data$converted_data(),
            sample_info = load_data$sample_info(),
            data_file_format = load_data$data_file_format(),
            no_id_conversion = input$no_id_conversion,
            log_transform_fpkm = input$log_transform_fpkm,
            log_start_fpkm = input$log_start_fpkm,
            low_filter_fpkm = input$low_filter_fpkm,
            missing_value = input$missing_value,
            scatter_x = input$scatter_x,
            scatter_y = input$scatter_y,
            sd_color = heat_colors[[input$heat_color_select]],
            rank = input$rank
          )
          save(params_r, file = file)
        }
        if (load_data$data_file_format() == 3) {
          params_r <- list(
            loaded_data = load_data$converted_data(),
            sample_info = load_data$sample_info(),
            data_file_format = load_data$data_file_format(),
            no_id_conversion = input$no_id_conversion,
            missing_value = input$missing_value,
            scatter_x = input$scatter_x,
            scatter_y = input$scatter_y,
            sd_color = heat_colors[[input$heat_color_select]],
            rank = input$rank,
            no_fdr = load_data$no_fdr()
          )
          save(params_r, file = file)
        }
      }
    )

    # Number of converted IDs ---------
    n_matched <- reactive({
      req(!is.null(processed_data))

      match_process_ids <-
        load_data$matched_ids() %in% rownames(processed_data()$data)

      return(sum(match_process_ids))
    })

    # Bias detected message -------
    read_counts_bias <- reactive({
      req(!is.null(processed_data()$raw_counts))

      counts_bias_message(
        raw_counts = processed_data()$raw_counts,
        data_file_format = load_data$data_file_format(),
        sample_info = load_data$sample_info()
      )
    })

    # Text Output Information -----------
    converted_message <- reactive({
      req(processed_data()$data_size)

      conversion_counts_message(
        data_size = processed_data()$data_size,
        all_gene_names = load_data$all_gene_names(),
        n_matched = n_matched()
      )
    })

    output$mapping_statistics_container <- renderUI({
      req(converted_message())

      tags$div(
        class = "mapping-statistics",
        br(),
        tags$p(converted_message(), style = "color: #B8860B;")
      )
    })

    # Sequencing depth warning -------
    observe({
      req(tab() == "Pre-Process")
      req(!is.null(read_counts_bias()))

      showNotification(
        ui = read_counts_bias(),
        id = "read_counts_message",
        duration = NULL,
        type = "warning"
      )
    })
    # Data type warning -------
    observe({
      req(tab() == "Pre-Process")
      req(processed_data()$data_type_warning != 0)

      message <- switch(as.character(processed_data()$data_type_warning),
        "1" = "Integers detected. Did you mean to select 'read counts'?",
        "-1" = "Non count values detected. Did you mean select 'Normalized Expression Values'?",
        "-2" = "A sample has all values as zero. it is recommended to remove that sample."
      )

      showNotification(
        ui = message,
        id = "data_type_warning",
        duration = NULL,
        type = "error"
      )
    })


    # Remove messages if the tab changes --------
    observe({
      req(tab() != "Pre-Process")

      removeNotification("read_counts_message")
      removeNotification("data_type_warning")
    })

    all_gene_info <- reactive({
      req(!is.null(load_data$converted()))

      return(
        get_gene_info(
          load_data$converted(),
          load_data$select_org(),
          gene_info_files = idep_data$gene_info_files
        )
      )
    })


    # Return Values -----------
    list(
      raw_counts = reactive(processed_data()$raw_counts),
      data = reactive(processed_data()$data),
      p_vals = reactive(processed_data()$p_vals),
      filter_size = reactive(processed_data()$data_size[3]),
      sample_info = reactive(load_data$sample_info()),
      all_gene_names = reactive(load_data$all_gene_names()),
      gmt_choices = reactive(load_data$gmt_choices()),
      converted = reactive(load_data$converted()),
      select_org = reactive(load_data$select_org()),
      gmt_file = reactive(load_data$gmt_file()),
      all_gene_info = reactive(load_data$all_gene_info()),
      data_file_format = reactive(load_data$data_file_format()),
      counts_log_start = reactive(input$counts_log_start),
      all_gene_info = reactive(all_gene_info()),
      descr = reactive(processed_data()$descr),
      heatmap_color_select = reactive(load_data$heatmap_color_select()),
      select_gene_id = reactive(load_data$select_gene_id()),
      plot_grid_lines = reactive(load_data$plot_grid_lines()),
      ggplot2_theme = reactive(load_data$ggplot2_theme())
    )
  })
}

## To be copied in the UI
# mod_02_pre_process_ui("pre_process") #nolint

## To be copied in the server
# mod_02_pre_process_server("pre_process_ui") #nolint
