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
                "Minimum counts per million. CPM adjusts for library size: (gene count / total reads) * 1,000,000.",
                theme = "light"
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
                "A gene must reach minimum CPM in this many samples to be kept for further analysis.",
                theme = "light"
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
            selected = 1,
            selectize = FALSE
          ),
          tippy::tippy_this(
            ns("counts_transform"),
            "Transformed data is used for clustering, PCA, and network analysis. Differential analysis with DESeq2 uses raw counts.",
            theme = "light"
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
                "Pseudo-count c:"
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
                  "Constant c in the log2(CPM + c) transformation; higher values reduce noise but decrease sensitivity. Typically between 1 and 10.",
                  theme = "light"
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
                "Minimum expression value (FPKM, RPKM, TPM, or similar normalized units). Defaults to -1000 to effectively bypass this filter for most datasets.",
                theme = "light"
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
                "A gene must reach the minimum expression level in at least this many samples to be kept for further analysis.",
                theme = "light"
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
                "Use log2-transformed values for all downstream analyses.",
                theme = "light"
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
            strong("Missing values:")
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
              selected = "geneMedian",
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("missing_value"),
              "You can treat missing values as zero, fill in with median values by gene in all samples, or use gene medians within each sample group.",
              theme = "light"
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
              "Download the transformed data. Includes Ensembl gene IDs and gene symbols.",
              theme = "light"
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
                "Download counts with gene IDs converted to Ensembl.",
                theme = "light"
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
          "Download the converted data as an .RData file.",
          theme = "light"
        ),
        downloadButton(
          outputId = ns("report"),
          label = tags$span(style = "color: red;", "Report")
        ),
        tippy::tippy_this(
          ns("report"),
          "Create an HTML report summarizing the preprocessing step.",
          theme = "light"
        ),
        uiOutput(ns("mapping_statistics_container")),
      ),


      # Pre-Process Panel Main -----------
      mainPanel(
        tabsetPanel(
          id = ns("eda_tabs"),

          # Barplot for read counts data ----------
          tabPanel(
            title = "Reads",
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

          # Distribution plots (combined Boxplot, Density, and Dispersion) ----------
          tabPanel(
            title = "Distribution",
            br(),
            fluidRow(
              column(
                width = 4,
                selectInput(
                  inputId = ns("distribution_plot_type"),
                  label = "Select Plot Type",
                  choices = c(
                    "Boxplot" = "boxplot",
                    "Density Plot" = "density",
                    "Dispersion" = "dispersion"
                  ),
                  selected = "boxplot",
                  selectize = FALSE
                )
              ),
              # Dispersion plot controls
              conditionalPanel(
                condition = "input.distribution_plot_type == 'dispersion'",
                column(
                  width = 4,
                  selectInput(
                    inputId = ns("heat_color_select"),
                    label = "Select Heat Colors",
                    choices = NULL,
                    selectize = FALSE
                  ),
                  tippy::tippy_this(
                    ns("heat_color_select"),
                    "Pick the color palette for the dispersion plot.",
                    theme = "light"
                  )
                ),
                column(
                  width = 4,
                  checkboxInput(
                    inputId = ns("rank"),
                    label = "Use rank of mean values"
                  ),
                  tippy::tippy_this(
                    ns("rank"),
                    "Rank genes by mean expression before plotting dispersion.",
                    theme = "light"
                  )
                ),
                ns = ns
              )
            ),
            br(),
            # Boxplot
            conditionalPanel(
              condition = "input.distribution_plot_type == 'boxplot'",
              plotOutput(
                outputId = ns("eda_boxplot"),
                width = "100%",
                height = "500px"
              ),
              ottoPlots::mod_download_figure_ui(
                id = ns("dl_eda_boxplot")
              ),
              ns = ns
            ),
            # Density Plot
            conditionalPanel(
              condition = "input.distribution_plot_type == 'density'",
              plotOutput(
                outputId = ns("eda_density"),
                width = "100%",
                height = "500px"
              ),
              ottoPlots::mod_download_figure_ui(
                id = ns("dl_eda_density")
              ),
              ns = ns
            ),
            # Dispersion Plot
            conditionalPanel(
              condition = "input.distribution_plot_type == 'dispersion'",
              plotOutput(
                outputId = ns("dev_transfrom"),
                width = "100%",
                height = "500px"
              ),
              ottoPlots::mod_download_figure_ui(
                id = ns("dl_dev_transform")
              ),
              ns = ns
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
                  selected = 1,
                  selectize = FALSE
                )
              ),
              column(
                width = 4,
                selectInput(
                  inputId = ns("scatter_y"),
                  label = tags$span("Sample for y-axis", style = "font-weight: normal;"),
                  choices = 1:5,
                  selected = 2,
                  selectize = FALSE
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
                      href = "http://useast.ensembl.org/info/genome/genebuild/biotypes.html",
                      onclick = "window.open(this.href, 'ensembl_info', 'width=900,height=700,resizable=yes,scrollbars=yes'); return false;",
                      rel = "noopener noreferrer"
                    )
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
                    div(
                      style = "display: flex; gap: 6px; align-items: center; flex-wrap: wrap;",
                      ottoPlots::mod_download_figure_ui(
                        id = ns("dl_rRNA_counts_gg")
                      )                      
                    )
                  ),
                  column(
                    3,
                    downloadButton(
                      outputId = ns("download_gene_type_data"),
                      label = "data"
                    )
                  ),
                  column(
                    7,
                    align = "right",
                    p("Unevenly high rRNA content may introduce bias.")
                  )
                ),
                tippy::tippy_this(
                  ns("download_gene_type_data"),
                  "Download read counts totals by gene type.",
                  theme = "light"
                ),
                uiOutput(ns("gene_type_warning")),
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
                    )
                  ),
                  column(
                    3,
                    uiOutput(ns("chr_boxplot_checkbox"))
                  ),
                  column(
                    7,
                    align = "right",
                    p("% reads mapped to a chromosome.")
                  )
                ),
                uiOutput(ns("chr_counts_warning")),
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
                  3,
                  uiOutput(ns("chr_normalized_boxplot_checkbox"))
                ),
                column(
                  7,
                  align = "right",
                  p("75th percentile of normalized expression by chromosomes.")
                )
              ),
              uiOutput(ns("chr_normalized_warning")),
            ns = ns
            )
          ),
          tabPanel(
            title = "Markers",
            br(),
            uiOutput(ns("markers_gene_selectors")),
            uiOutput(ns("markers_plots"))
          ),
          # Plot panel for individual genes ---------
          tabPanel(
            title = "Genes",
            br(),
            fluidRow(
              column(
                4,
                # Gene ID Selection -----------
                selectizeInput(
                  inputId = ns("selected_gene"),
                  label = "Search for Gene(s)",
                  choices = "",
                  selected = NULL,
                  multiple = TRUE
                ),
                tippy::tippy_this(
                  ns("selected_gene"),
                  "Type to search by gene ID, symbol, or name.",
                  theme = "light"
                )
              ),
              column(
                4,
                selectInput(
                  inputId = ns("gene_plot_box"),
                  label = "Plot Type",
                  choices = setNames(
                    c(1,2),
                    c("Sample Groups",
                      "Individual Samples")
                  ),
                  selected = 1,
                  selectize = FALSE
                ),
                tippy::tippy_this(
                  ns("gene_plot_box"),
                  "Choose whether to plot sample groups or individual samples.",
                  theme = "light"
                ),
                uiOutput(ns("sd_checkbox")),
                conditionalPanel(
                  condition = "output.data_file_format == 1",
                  checkboxInput(
                    inputId = ns("plot_raw"),
                    label = "Raw counts",
                    value = FALSE
                  ),
                  tippy::tippy_this(
                    ns("plot_raw"),
                    "Toggle between raw counts and transformed values.",
                    theme = "light"
                  ),
                  ns = ns
                ),
                checkboxInput(
                  inputId = ns("plot_tukey"),
                  label = "TukeyHSD test",
                  value = FALSE
                ),
                tippy::tippy_this(
                  ns("plot_tukey"),
                  "Run Tukey's post-hoc test for pairwise group comparisons (only works for transformed group plots).",
                  theme = "light"
                )
              ),
              column(
                4,
                radioButtons(
                  inputId = ns("angle_ind_axis_lab"),
                  label = "Rotate Labels",
                  choices = c(0, 45, 90),
                  selected = 45
                ),
                tippy::tippy_this(
                  ns("angle_ind_axis_lab"),
                  "Set the rotation angle of the x-axis labels.",
                  theme = "light"
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
                label = "TukeyHSD"
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
                label = "Show raw counts",
                value = FALSE
              ),
              tippy::tippy_this(
                ns("show_raw"),
                "Display the original uploaded data instead of transformed values.",
                theme = "light"
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
        hideTab(inputId = "eda_tabs", target = "Reads")
        updateTabsetPanel(session, "eda_tabs", selected = "Distribution")
      } else if (load_data$data_file_format() == 1) {
        showTab(inputId = "eda_tabs", target = "Reads")
        updateTabsetPanel(session, "eda_tabs", selected = "Reads")
      }
    })

    observeEvent(list(input$eda_tabs, input$distribution_plot_type), {
      req(input$eda_tabs == "Distribution")
      req(input$distribution_plot_type == "density")
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

    # Dynamic Markers Tab - hide if no markers found ----------
    observe({
      req(!is.null(processed_data()$data))
      payloads <- marker_payloads()

      if (length(payloads) == 0) {
        hideTab(inputId = "eda_tabs", target = "Markers")
      } else {
        showTab(inputId = "eda_tabs", target = "Markers")
      }
    })

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
      req(!tab() %in% c("Cluster", "Data",
                        "Network", "Bicluster",
                        "DEG"))
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
      req(!is.null(load_data$converted_data()))
      req(!is.null(load_data$all_gene_info()))
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
    rRNA_counts_result <- reactive({
      req(!is.null(processed_data()$raw_counts))
      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Plotting counts by gene type",
        color = "#000000"
      )
      result <- rRNA_counts_ggplot(
        counts_data = load_data$converted_data(),
        sample_info = load_data$sample_info(),
        type = "Raw",
        all_gene_info = load_data$all_gene_info(),
        plots_color_select = load_data$plots_color_select()
      )
      shinybusy::remove_modal_spinner()
      return(result)
    })

    rRNA_counts <- reactive({
      req(!is.null(rRNA_counts_result()))
      result <- rRNA_counts_result()
      p <- refine_ggplot2(
        p = result$plot,
        gridline = load_data$plot_grid_lines(),
        ggplot2_theme = load_data$ggplot2_theme()
      )
      return(p)
    })
    output$rRNA_counts_gg <- renderPlot({
      print(rRNA_counts())
    })

    output$gene_type_warning <- renderUI({
      req(!is.null(rRNA_counts_result()))
      result <- rRNA_counts_result()

      if (!is.null(result$warning)) {
        tags$div(
          style = "padding: 10px; margin: 10px 0; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 4px; color: #856404;",
          tags$strong("\u26A0 Warning: "),
          result$warning
        )
      }
    })

    dl_rRNA_counts_gg <- ottoPlots::mod_download_figure_server(
      id = "dl_rRNA_counts_gg",
      filename = "rRNA_counts_barplot",
      figure = reactive({
        rRNA_counts()
      }),
      label = ""
    )

    gene_type_totals <- reactive({
      req(load_data$data_file_format() == 1)
      req(!is.null(load_data$converted_data()))
      req(!is.null(load_data$all_gene_info()))

      counts_data <- load_data$converted_data()
      gene_info <- load_data$all_gene_info()

      df <- merge(
        counts_data,
        gene_info,
        by.x = "row.names",
        by.y = "ensembl_gene_id"
      )

      df$gene_biotype <- gsub(".*pseudogene", "Pseudogene", df$gene_biotype)
      df$gene_biotype <- gsub("TEC", "Unknown", df$gene_biotype)
      df$gene_biotype <- gsub("IG_.*", "IG", df$gene_biotype)
      df$gene_biotype <- gsub("TR_.*", "TR", df$gene_biotype)
      df$gene_biotype <- gsub("protein_coding", "Coding", df$gene_biotype)

      counts_by_type <- aggregate(
        df[, colnames(counts_data), drop = FALSE],
        by = list(df$gene_biotype),
        FUN = sum
      )
      colnames(counts_by_type)[1] <- "Gene_Type"

      sample_counts <- counts_by_type[, -1, drop = FALSE]
      sample_totals <- colSums(sample_counts, na.rm = TRUE)

      percent_matrix <- sample_counts * 0
      non_zero_samples <- sample_totals != 0

      if (any(non_zero_samples)) {
        percent_matrix[, non_zero_samples] <- sweep(
          sample_counts[, non_zero_samples, drop = FALSE],
          2,
          sample_totals[non_zero_samples],
          "/"
        ) * 100
      }

      percent_matrix[, !non_zero_samples] <- 0

      percent_matrix <- round(percent_matrix, 2)
      row_avg <- rowMeans(percent_matrix, na.rm = TRUE)
      order_idx <- order(row_avg, decreasing = TRUE)

      data.frame(
        `Gene Type` = counts_by_type$Gene_Type,
        percent_matrix,
        check.names = FALSE
      )[order_idx, ]
    })

    output$download_gene_type_data <- downloadHandler(
      filename = function() {
        "gene_type_read_totals.csv"
      },
      content = function(file) {
        req(gene_type_totals())
        utils::write.csv(gene_type_totals(), file, row.names = FALSE)
      }
    )

    # Dynamic checkbox for chr counts plot ------------
    output$chr_boxplot_checkbox <- renderUI({
      req(!is.null(load_data$converted_data()))

      # Determine default value based on number of samples
      n_samples <- ncol(load_data$converted_data())
      default_value <- n_samples > 50

      tagList(
        checkboxInput(
          inputId = ns("chr_use_boxplot"),
          label = "Use boxplot",
          value = default_value
        ),
        tippy::tippy_this(
          ns("chr_use_boxplot"),
          "Show boxplot grouped by sample groups instead of barplot.",
          theme = "light"
        )
      )
    })

    # chr counts barplot ------------
    chr_counts_result <- reactive({
      req(!is.null(processed_data()$raw_counts))
      req(!is.null(input$chr_use_boxplot))
      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Plotting counts by Chromosome",
        color = "#000000"
      )
      result <- chr_counts_ggplot(
        counts_data = load_data$converted_data(),
        sample_info = load_data$sample_info(),
        type = "Raw",
        all_gene_info = load_data$all_gene_info(),
        plots_color_select = load_data$plots_color_select(),
        use_boxplot = input$chr_use_boxplot
      )
      shinybusy::remove_modal_spinner()
      return(result)
    })

    chr_counts <- reactive({
      req(!is.null(chr_counts_result()))
      result <- chr_counts_result()
      p <- refine_ggplot2(
        p = result$plot,
        gridline = load_data$plot_grid_lines(),
        ggplot2_theme = load_data$ggplot2_theme()
      )
      return(p)
    })
    output$chr_counts_gg <- renderPlot({
      print(chr_counts())
    },
    height = 2000)

    # chr counts warning message output
    output$chr_counts_warning <- renderUI({
      req(!is.null(chr_counts_result()))
      result <- chr_counts_result()

      if (!is.null(result$warning)) {
        tags$div(
          style = "padding: 10px; margin: 10px 0; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 4px; color: #856404;",
          tags$strong("\u26A0 Warning: "),
          result$warning
        )
      }
    })

    dl_chr_counts_gg <- ottoPlots::mod_download_figure_server(
      id = "dl_chr_counts_gg",
      filename = "Chr_counts_barplot",
      figure = reactive({
        chr_counts()
      }),
      label = ""
    )

    # Dynamic checkbox for chr normalized plot ------------
    output$chr_normalized_boxplot_checkbox <- renderUI({
      req(!is.null(processed_data()$data))

      # Determine default value based on number of samples
      n_samples <- ncol(processed_data()$data)
      default_value <- n_samples > 50

      tagList(
        checkboxInput(
          inputId = ns("chr_normalized_use_boxplot"),
          label = "Use boxplot",
          value = default_value
        ),
        tippy::tippy_this(
          ns("chr_normalized_use_boxplot"),
          "Show boxplot grouped by sample groups instead of barplot.",
          theme = "light"
        )
      )
    })

    # chr normalized barplot ------------
    chr_normalized_result <- reactive({
      req(!is.null(processed_data()$data))
      req(!is.null(input$chr_normalized_use_boxplot))
      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Plotting normalized expression by Chromosome",
        color = "#000000"
      )
      result <- chr_normalized_ggplot(
        counts_data = processed_data()$data,
        sample_info = load_data$sample_info(),
        type = "Raw",
        all_gene_info = load_data$all_gene_info(),
        plots_color_select = load_data$plots_color_select(),
        use_boxplot = input$chr_normalized_use_boxplot
      )
      shinybusy::remove_modal_spinner()
      return(result)
    })

    chr_normalized <- reactive({
      req(!is.null(chr_normalized_result()))
      result <- chr_normalized_result()
      p <- refine_ggplot2(
        p = result$plot,
        gridline = load_data$plot_grid_lines(),
        ggplot2_theme = load_data$ggplot2_theme()
      )
      return(p)
    })
    output$chr_normalized_gg <- renderPlot({
      print(chr_normalized())
    },
    height = 2000)

    output$chr_normalized_warning <- renderUI({
      req(!is.null(chr_normalized_result()))
      result <- chr_normalized_result()

      if (!is.null(result$warning)) {
        tags$div(
          style = "padding: 10px; margin: 10px 0; background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 4px; color: #856404;",
          tags$strong("\u26A0 Warning: "),
          result$warning
        )
      }
    })

    dl_chr_counts_gg <- ottoPlots::mod_download_figure_server(
      id = "dl_chr_normalized_gg",
      filename = "Chr_normalized_expression_barplot",
      figure = reactive({
        chr_normalized()
      }),
      label = ""
    )

    marker_payloads <- reactive({
      expr_matrix <- processed_data()$data
      if (is.null(expr_matrix)) {
        return(list())
      }
      marker_gene_plot_payloads(
        expr_matrix = expr_matrix,
        sample_info = load_data$sample_info(),
        raw_counts = processed_data()$raw_counts,
        converted_data = load_data$converted_data(),
        all_gene_info = load_data$all_gene_info()
      )
    })

    output$markers_gene_selectors <- renderUI({
      req(tab() == "Prep")
      NULL
    })

    output$markers_plots <- renderUI({
      payloads <- marker_payloads()
      if (length(payloads) == 0) {
        return(tags$p("Marker gene plots are unavailable for the current dataset."))
      }
      tagList(
        lapply(seq_along(payloads), function(idx) {
          marker <- payloads[[idx]]
          plot_id <- paste0("marker_plot_", idx)
          section_id <- paste0("marker_section_", idx)
          local({
            local_marker <- marker
            local_plot_id <- plot_id
            local_section_id <- section_id
            local_idx <- idx

            output[[local_section_id]] <- renderUI({
              payload_list <- marker_payloads()
              if (length(payload_list) < local_idx) {
                return(NULL)
              }
              marker_payload <- payload_list[[local_idx]]
              if (is.null(marker_payload$data)) {
                return(NULL)
              }
              tagList(
                br(),
                h4(paste0(marker_payload$description)),
                mod_gene_expression_plot_ui(
                  id = ns(local_plot_id),
                  plot_height = "400px",
                  show_download = TRUE
                ),
                hr()
              )
            })

            mod_gene_expression_plot_server(
              id = local_plot_id,
              plot_data = reactive({
                payload_list <- marker_payloads()
                if (length(payload_list) < local_idx) {
                  return(NULL)
                }
                payload <- payload_list[[local_idx]]
                if (is.null(payload$data)) {
                  return(NULL)
                }
                list(
                  data = payload$data,
                  display_name = payload$display_name
                )
              }),
              palette_name = reactive(load_data$plots_color_select()),
              plot_grid_lines = reactive(load_data$plot_grid_lines()),
              ggplot2_theme = reactive(load_data$ggplot2_theme()),
              counts_are_counts = reactive(!is.null(processed_data()$raw_counts)),
              download_filename = paste0(local_marker$symbol, "_expression_barplot"),
              default_plot_type = "bar",
              default_data_type = "raw"
            )

            uiOutput(ns(local_section_id))
          })
        })
      )
    })

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

      tagList(
        checkboxInput(
          inputId = ns("use_sd"),
          label = "Standard deviation",
          value = FALSE
        ),
        tippy::tippy_this(
          ns("use_sd"),
          "Show standard deviation (SD) or standard error (SE) for grouped samples.",
          theme = "light"
        )
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
      filename = paste0(
        "pre_process_report_",
        format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
        ".html"
      ),
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

          # Persist current chromosome plot selections for report rendering
          chr_use_boxplot <- input$chr_use_boxplot
          if (is.null(chr_use_boxplot)) {
            chr_use_boxplot <- if (!is.null(load_data$converted_data())) {
              ncol(load_data$converted_data()) > 50
            } else {
              FALSE
            }
          }
          chr_normalized_use_boxplot <- input$chr_normalized_use_boxplot
          if (is.null(chr_normalized_use_boxplot)) {
            chr_normalized_use_boxplot <- if (!is.null(processed_data()$data)) {
              ncol(processed_data()$data) > 50
            } else {
              FALSE
            }
          }

          individual_data_current <- individual_data()
          
          safe_input <- function(value, default) {
            if (is.null(value) || length(value) == 0) {
              default
            } else {
              value
            }
          }
          
          selected_gene <- input$selected_gene
          if (is.null(selected_gene) || length(selected_gene) == 0) {
            if (!is.null(individual_data_current) && nrow(individual_data_current) > 0) {
              sorted <- sort(
                apply(
                  individual_data_current,
                  MARGIN = 1,
                  FUN = function(x) sd(x)
                ),
                decreasing = TRUE
              )
              if (length(sorted) > 0) {
                selected_gene <- names(sorted)[seq_len(min(2, length(sorted)))]
              } else {
                selected_gene <- character(0)
              }
            } else {
              selected_gene <- character(0)
            }
          }
          
          gene_plot_box <- safe_input(input$gene_plot_box, 1)
          use_sd <- safe_input(input$use_sd, FALSE)
          lab_rotate <- safe_input(input$angle_ind_axis_lab, 45)
          plot_raw <- safe_input(input$plot_raw, FALSE)
          plot_tukey <- safe_input(input$plot_tukey, FALSE)
          
          # Set up parameters to pass to Rmd document
          params <- list(
            loaded_data = load_data$converted_data(),
            individual_data = individual_data_current,
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
            n_min_samples_fpkm = input$n_min_samples_fpkm,
            missing_value = input$missing_value,
            scatter_x = input$scatter_x,
            scatter_y = input$scatter_y,
            sd_color = heat_colors[[input$heat_color_select]],
            rank = input$rank,
            no_fdr = load_data$no_fdr(),
            selected_gene = selected_gene,
            gene_plot_box = gene_plot_box,
            use_sd = use_sd,
            lab_rotate = lab_rotate,
            plot_raw = plot_raw,
            plot_tukey = plot_tukey,
            plots_color_select = load_data$plots_color_select(),
            plot_grid_lines = load_data$plot_grid_lines(),
            ggplot2_theme = load_data$ggplot2_theme(),
            chr_use_boxplot = chr_use_boxplot,
            chr_normalized_use_boxplot = chr_normalized_use_boxplot,
            mapping_statistics = converted_message()
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
      req(tab() == "Prep")
      req(!is.null(read_counts_bias()))

      showNotification(
        ui = read_counts_bias(),
        id = "read_counts_message",
        duration = NULL,
        type = "error"
      )
    })
    # Data type warning -------
    observe({
      req(tab() == "Prep")
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
      req(tab() != "Prep")

      removeNotification("read_counts_message")
      removeNotification("data_type_warning")
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
      descr = reactive(processed_data()$descr),
      mapping_statistics = reactive(converted_message()),
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
