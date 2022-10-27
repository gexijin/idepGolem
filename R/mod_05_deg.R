#' 05_deg1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_05_deg_1_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "DEG1",
    sidebarLayout(
      sidebarPanel(
        # Button to run DEG analysis for the specified model
        uiOutput(ns("submit_ui")),
        tags$head(tags$style(
          "#deg-submit_model_button{font-size: 16px;color: red}"
        )),
        br(),
        br(),
        # DEG analysis methods for read counts data
        conditionalPanel(
          condition = "output.data_file_format == 1",
          selectInput(
            inputId = ns("counts_deg_method"),
            label = "Method:",
            choices = list(
              "DESeq2" = 3,
              "limma-voom" = 2,
              "limma-trend" = 1
            ),
            selected = 3
          ),
          tags$style(
            type = "text/css",
            "#deg-counts_deg_method {width:100%;   margin-top:-12px}"
          ),
          ns = ns
        ),
        # Label when the limma method is selected
        conditionalPanel(
          condition = "output.data_file_format == 2",
          h5("Using the limma package"),
          ns = ns
        ),
        fluidRow(
          column(
            width = 5,
            # Adjusted significant p-value to use
            numericInput(
              inputId = ns("limma_p_val"),
              label = "FDR cutoff",
              value = 0.1,
              min = 1e-5,
              max = 1,
              step = .05
            )
          ),
          column(
            width = 7,
            # Min fold change to use
            numericInput(
              inputId = ns("limma_fc"),
              label = "Min fold-change",
              value = 2,
              min = 1,
              max = 100,
              step = 0.5
            )
          ),
          # Style both numeric inputs
          tags$style(
            type = "text/css",
            "#deg-limma_p_val { width:100%;   margin-top:-12px}"
          ),
          tags$style(
            type = "text/css",
            "#deg-limma_fc { width:100%;   margin-top:-12px}"
          )
        ),
        conditionalPanel(
          condition = "input.counts_deg_method == 3 && output.data_file_format == 1",
          checkboxInput(
            inputId = ns("threshold_wald_test"),
            label = "Threshold-based Wald Test",
            value = FALSE
          ),
          checkboxInput(
            inputId = ns("independent_filtering"),
            label = "Independent filtering of lower counts",
            value = TRUE
          ),
          ns = ns
        ),
        tags$br(),
        tags$br(),
        uiOutput(ns("download_lfc_button")),
        uiOutput(ns("note_download_lfc_button")),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/degs/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          id = ns("step_1"),
          tabPanel(
            title = "Experiment Design",
            fluidRow(
              column(
                width = 6,
                htmlOutput(outputId = ns("list_factors_deg"))
              ),
              column(
                width = 6,
                htmlOutput(outputId = ns("list_block_factors_deg"))
              )
            ),
            fluidRow(
              htmlOutput(outputId = ns("select_reference_levels"))
            ),
            htmlOutput(outputId = ns("list_interaction_terms")),
            textOutput(outputId = ns("experiment_design")),
            tags$head(tags$style(
              "#deg-experiment_design{color: red;font-size: 16px;}"
            )),
            htmlOutput(outputId = ns("list_model_comparisons")),
            h3("Use the submit button in the sidebar once the desired
               design is selected!"),
            a(
              h5("More info on DESeq2 experiment design", align = "right"),
              href = "http://rpubs.com/ge600/deseq2",
              target = "_blank"
            )
          ),
          tabPanel(
            title = "Results",
            value = ("results_tab"),
            plotOutput(
              outputId = ns("sig_gene_stats")
            ),
            br(),
            br(),
            h5(
              "Numbers of differentially expressed genes for all comparisons.
              \"B-A\" means B vs. A. Interaction terms start with \"I:\" "
            ),
            tableOutput(
              outputId = ns("sig_gene_stats_table")
            )
          ),
          tabPanel(
            title = "Venn Diagram & UpSet plot",
            checkboxInput(
              inputId = ns("up_down_regulated"),
              label = "Split gene lists by up- or down-regulation",
              value = TRUE
            ),
            htmlOutput(outputId = ns("list_comparisons_venn")),
            plotOutput(outputId = ns("venn_plot")),
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_venn")
            ),
            plotOutput(outputId = ns("upset_plot")),
            ottoPlots::mod_download_figure_ui(id = ns("dl_upset")),
            tags$p("The above graph is an UpSet plot that is an alternative to a
            venn diagram. The plot shows the intersections of the data in the
            combination matrix (bottom) and the columns show how many genes are
            in each intersection."),
            tags$a(
              h5("More info on plot", align = "right"),
              href = "https://en.wikipedia.org/wiki/UpSet_Plot#:~:text=UpSet%20plots%20are%20a%20data,sets%20(or%20vice%20versa).",
              target = "_blank"
            )
          ),
          tabPanel(
            title = "R Code",
            downloadButton(
              outputId = ns("dl_deg_code"),
              label = "Code"
            ),
            tippy::tippy_this(
              ns("dl_deg_code"),
              "Download .R file of DEG code",
              theme = "light-border"
            ),
            verbatimTextOutput(
              ns("deg_code")
            )
          )
        )
      )
    )
  )
}

mod_05_deg_2_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "DEG2",
    sidebarLayout(
      sidebarPanel(
        htmlOutput(outputId = ns("list_comparisons")),
        p("Select a comparison to examine the associated DEGs.
          \"A-B\" means A vs. B (See heatmap).
            Interaction terms start with \"I:\""),
        conditionalPanel(
          condition = "input.step_2 == 'Volcano Plot' |
            input.step_2 == 'MA Plot'",
          selectInput(
            inputId = ns("plot_color_select"),
            label = "Color scale",
            choices = "Red-Green"
          ),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.step_2 == 'Volcano Plot' ",
          mod_label_ui(ns("label_volcano")),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.step_2 == 'MA Plot' ",
          mod_label_ui(ns("label_ma")),
          ns = ns
        ),
        width = 2
      ),
      mainPanel(
        tabsetPanel(
          id = ns("step_2"),
          tabPanel(
            title = "Heatmap",
            mod_12_heatmap_ui(ns("12_heatmap_1"))
          ),
          tabPanel(
            title = "Volcano Plot",
            br(),
            plotOutput(
              outputId = ns("volcano_plot"),
              height = "500px",
              width = "100%"
            ),
            ottoPlots::mod_download_figure_ui(ns("download_volcano"))
          ),
          tabPanel(
            title = "MA Plot",
            br(),
            plotOutput(
              outputId = ns("ma_plot"),
              height = "500px",
              width = "100%"
            ),
            ottoPlots::mod_download_figure_ui(ns("download_ma"))
          ),
          tabPanel(
            title = "Scatter Plot",
            br(),
            plotOutput(
              outputId = ns("scatter_plot"),
              height = "500px",
              width = "100%"
            )
          ),
          tabPanel(
            title = "Enrichment",
            br(),
            mod_11_enrichment_ui(ns("enrichment_table_cluster")),
          )
        )
      )
    )
  )
}

#' 05_deg1 Server Functions
#'
#' @noRd
mod_05_deg_server <- function(id, pre_process, idep_data, load_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Interactive heatmap environment
    deg_env <- new.env()

    output$submit_ui <- renderUI({
      req(model_comparisons())
      tagList(
        actionButton(
          inputId = ns("submit_model_button"),
          label = "Submit",
          style = "float:right"
        ),
        tippy::tippy_this(
          ns("submit_model_button"),
          "Run DEG analysis",
          theme = "light-border"
        )
      )
    })

    # DEG STEP 1 ----------
    output$data_file_format <- reactive({
      pre_process$data_file_format()
    })
    outputOptions(
      x = output,
      name = "data_file_format",
      suspendWhenHidden = FALSE
    )

    # Experiment Design UI Elements ------------
    output$list_factors_deg <- renderUI({
      list_factors <- list_factors_ui(
        sample_info = pre_process$sample_info(),
        data_file_format = pre_process$data_file_format(),
        counts_deg_method = input$counts_deg_method
      )

      if (class(list_factors)[1] == "list") {
        return(
          checkboxGroupInput(
            inputId = ns("select_factors_model"),
            h5(list_factors$title),
            choices = list_factors$choices,
            selected = list_factors$choices[1]
          )
        )
      } else {
        return(list_factors)
      }
    })

    output$list_block_factors_deg <- renderUI({
      choices <- list_block_factors_ui(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model
      )
      req(!is.null(choices))
      return(
        checkboxGroupInput(
          inputId = ns("select_block_factors_model"),
          h5("Select a factor for batch effect or paired samples, if needed."),
          choices = choices,
          selected = NULL
        )
      )
    })

    model_comparisons <- reactive({
      req(pre_process$data() & pre_process$data_file_format() != 3)

      list_model_comparisons_ui(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model,
        processed_data = pre_process$data()
      )
    })

    output$list_model_comparisons <- renderUI({
      req(model_comparisons())
      checkboxGroupInput(
        inputId = ns("select_model_comprions"),
        label = h5(model_comparisons()$title),
        choices = model_comparisons()$choices,
        selected = model_comparisons()$choices[[1]]
      )
    })

    output$list_interaction_terms <- renderUI({
      interactions <- list_interaction_terms_ui(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model
      )
      req(!is.null(interactions))
      checkboxGroupInput(
        inputId = ns("select_interactions"),
        label = h5(
          "Interaction terms between factors(e.g. genotypes repond differently
          to treatment?):"
        ),
        choices = interactions,
        selected = NULL
      )
    })


    # Set limits for selections of factors
    observe({
      if (length(input$select_factors_model) > 6) {
        updateCheckboxGroupInput(
          session,
          ns("select_factors_model"),
          selected = tail(input$select_factors_model, 6)
        )
      }
      if (input$counts_deg_method != 3) {
        if (length(input$select_factors_model) > 2) {
          updateCheckboxGroupInput(
            session,
            ns("select_factors_model"),
            selected = tail(input$select_factors_model, 2)
          )
        }
        if (length(input$select_block_factors_model) > 1) {
          updateCheckboxGroupInput(
            session,
            ns("select_block_factors_model"),
            selected = tail(input$select_block_factors_model, 1)
          )
        }
      }

      if (length(input$select_comparisons_venn) > 5) {
        updateCheckboxGroupInput(
          session,
          ns("select_comparisons_venn"),
          selected = tail(input$select_comparisons_venn, 5)
        )
      }
    })

    output$experiment_design <- renderText({
      experiment_design_txt(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model,
        select_block_factors_model = input$select_block_factors_model,
        select_interactions = input$select_interactions
      )
    })

    output$select_reference_levels <- renderUI({
      select_choices <- select_reference_levels_ui(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model,
        data_file_format = pre_process$data_file_format(),
        counts_deg_method = input$counts_deg_method
      )
      req(!is.null(select_choices))
      lapply(names(select_choices), function(x) {
        tagList(
          column(
            width = 4,
            selectInput(
              inputId = ns(
                paste0(
                  "reference_level_factor_",
                  which(names(select_choices) == x)
                )
              ),
              label = h5(paste0("Reference/baseline level for ", x)),
              choices = setNames(
                as.list(paste0(x, ":", select_choices[[x]])),
                select_choices[[x]]
              )
            )
          )
        )
      })
    })

    factor_reference_levels <- reactive(
      return(
        c(
          input$reference_level_factor_1,
          input$reference_level_factor_2,
          input$reference_level_factor_3,
          input$reference_level_factor_4,
          input$reference_level_factor_5,
          input$reference_level_factor_6
        )
      )
    )

    # Observe submit button ------
    deg <- reactiveValues(limma = NULL)
    observeEvent(
      input$submit_model_button,
      {
        req(!is.null(pre_process$raw_counts()) |
          !is.null(pre_process$data()))
        withProgress(message = "DEG analysis...", {
          incProgress(0.4)

          # only use with DESeq2
          threshold_wald_test <- FALSE
          if (input$counts_deg_method == 3) {
            threshold_wald_test <- input$threshold_wald_test
          }

          independent_filtering <- TRUE
          if (input$counts_deg_method == 3) {
            independent_filtering <- input$independent_filtering
          }

          deg$limma <- limma_value(
            data_file_format = pre_process$data_file_format(),
            counts_deg_method = input$counts_deg_method,
            raw_counts = pre_process$raw_counts(),
            limma_p_val = input$limma_p_val,
            limma_fc = input$limma_fc,
            select_model_comprions = input$select_model_comprions,
            sample_info = pre_process$sample_info(),
            select_factors_model = input$select_factors_model,
            select_interactions = input$select_interactions,
            select_block_factors_model = input$select_block_factors_model,
            factor_reference_levels = factor_reference_levels(),
            processed_data = pre_process$data(),
            counts_log_start = pre_process$counts_log_start(),
            p_vals = pre_process$p_vals(),
            threshold_wald_test = threshold_wald_test,
            independent_filtering = independent_filtering
          )

          updateTabsetPanel(
            session = session,
            inputId = "step_1",
            selected = "results_tab"
          )
        })
      }
    )

    deg_info <- reactive({
      req(!is.null(deg$limma$results))

      deg_information(
        limma_value = deg$limma,
        gene_names = pre_process$all_gene_names(),
        processed_data = pre_process$data(),
        no_id_conversion = load_data$no_id_conversion()
      )[[1]]
    })

    deg_method <- c(
      "limma_trend",
      "limma_voom",
      "DESeq2"
    )

    name <- reactive({
      paste0(
        "deg_values_",
        deg_method[as.numeric(input$counts_deg_method)],
        ".csv"
      )
    })

    output$download_lfc <- downloadHandler(
      filename = function() {
        name()
      },
      content = function(file) {
        write.csv(deg_info(), file, row.names = FALSE)
      }
    )

    output$download_lfc_button <- renderUI({
      req(!is.null(deg_info()))
      downloadButton(
        outputId = ns("download_lfc"),
        "Results & data"
      )
    })

    output$note_download_lfc_button <- renderUI({
      req(!is.null(deg_info()))
      tippy::tippy_this(
        elementId = ns("download_lfc"),
        tooltip = "This data includes log fold change, adjusted p-value and
            processed data from Pre-Process tab.",
        theme = "light-border"
      )
    })

    output$sig_gene_stats <- renderPlot({
      req(!is.null(deg$limma$results))
      p <- sig_genes_plot(
        results = deg$limma$results
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })

    output$sig_gene_stats_table <- renderTable(
      {
        req(!is.null(deg$limma))

        genes_stat_table(limma = deg$limma)
      },
      digits = 0,
      spacing = "s",
      include.rownames = FALSE,
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = T
    )

    output$deg_code <- renderText({
      req(!is.null(deg$limma))
      deg$limma$expr
    })

    output$dl_deg_code <- downloadHandler(
      filename = function() {
        "DEG_code.R"
      },
      content = function(file) {
        if (is.null(deg$limma)) {
          writeLines(" Nothing!", file)
        } else {
          writeLines(deg$limma$expr, file)
        }
      }
    )


    output$list_comparisons_venn <- renderUI({
      req(!is.null(deg$limma))

      venn_comp <- list_comp_venn(
        limma = deg$limma,
        up_down_regulated = input$up_down_regulated
      )
      if (is.null(venn_comp$choices)) {
        selectInput(
          inputId = ns("select_comparisons_venn"),
          label = NULL,
          choices = list("All" = "All"),
          selected = "All"
        )
      } else {
        tagList(
          h4("Select comparisons:"),
          h6("The venn diagram will only display up to 5 comparisons"),
          checkboxGroupInput(
            inputId = ns("select_comparisons_venn"),
            label = NULL,
            choices = venn_comp$choices,
            selected = venn_comp$choices_first_three
          )
        )
      }
    })

    # venn diagram -----
    venn_data <- reactive({
      req(!is.null(deg$limma))
      req(!is.null(input$select_comparisons_venn))
      req(input$up_down_regulated)

      prep_venn(
        limma = deg$limma,
        up_down_regulated = input$up_down_regulated,
        select_comparisons_venn = input$select_comparisons_venn
      )
    })

    venn <- reactive({
      req(!is.null(deg$limma))
      req(!is.null(input$select_comparisons_venn))

      venn <- plot_venn(
        results = venn_data()
      )
      p <- recordPlot()
      return(p)
    })
    output$venn_plot <- renderPlot({
      print(venn())
    })
    dl_venn <- ottoPlots::mod_download_figure_server(
      id = "dl_venn",
      filename = "venn_diagram",
      figure = reactive({
        venn()
      }),
      label = ""
    )

    upset <- reactive({
      p <- plot_upset(
        results = venn_data()
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })

    output$upset_plot <- renderPlot({
      print(upset())
    })

    dl_upset <- ottoPlots::mod_download_figure_server(
      id = "dl_upset",
      filename = "upset_plot",
      figure = reactive({
        upset()
      }),
      label = ""
    )


    # DEG STEP 2 --------
    output$list_comparisons <- renderUI({
      if (is.null(deg$limma$comparisons)) {
        selectInput(
          inputId = ns("select_contrast"),
          label = NULL,
          choices = list("All" = "All"),
          selected = "All"
        )
      } else {
        selectInput(
          inputId = ns("select_contrast"),
          label = NULL,
          choices = deg$limma$comparisons
        )
      }
    })

    contrast_samples <- reactive({
      req(!is.null(input$select_contrast))

      find_contrast_samples(
        select_contrast = input$select_contrast,
        all_sample_names = colnames(pre_process$data()),
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model,
        select_model_comprions = input$select_model_comprions,
        reference_levels = factor_reference_levels(),
        counts_deg_method = input$counts_deg_method,
        data_file_format = pre_process$data_file_format()
      )
    })

    heat_data <- reactive({
      req(!is.null(deg$limma))
      req(!is.null(input$select_contrast))

      deg_heat_data(
        limma = deg$limma,
        select_contrast = input$select_contrast,
        processed_data = pre_process$data(),
        contrast_samples = contrast_samples()
      )
    })

    # bar to make the heatmap module reactive
    # otherwise, error when switching heatmap
    heatmap_bar <- reactive({
      req(!is.null(heat_data()))
      heat_data()$bar
    })

    heatmap_module <- mod_12_heatmap_server(
      id = "12_heatmap_1",
      data = reactive({
        heat_data()$genes
      }),
      bar = reactive({
        heatmap_bar()
      }),
      all_gene_names = reactive({
        pre_process$all_gene_names()
      }),
      cluster_rows = FALSE,
      heatmap_color = reactive({
        pre_process$heatmap_color_select()
      }),
      select_gene_id = reactive({
        pre_process$select_gene_id()
      })
    )

    # Plot colors -------
    plot_colors <- list(
      "Green-Red" = c("green", "grey45", "red"),
      "Red-Green" = c("red", "grey45", "green"),
      "Blue-Red" = c("blue", "grey45", "red"),
      "Green-Magenta" = c("green", "grey45", "magenta"),
      "Orange-Blue" = c("orange", "grey45", "blue")
    )

    plot_choices <- c(
      "Green-Red",
      "Red-Green",
      "Blue-Red",
      "Green-Magenta",
      "Orange-Blue"
    )

    observe({
      updateSelectInput(
        session = session,
        inputId = "plot_color_select",
        choices = plot_choices
      )
    })

    # volcano plot -----
    vol_data <- reactive({
      req(input$select_contrast, deg$limma, input$limma_p_val, input$limma_fc)

      volcano_data(
        select_contrast = input$select_contrast,
        comparisons = deg$limma$comparisons,
        top_genes = deg$limma$top_genes,
        limma_p_val = input$limma_p_val,
        limma_fc = input$limma_fc,
        processed_data = pre_process$data(),
        contrast_samples = contrast_samples()
      )
    })

    gene_labels <- mod_label_server(
      "label_volcano",
      data_list = reactive({
        vol_data()
      }),
      method = "volcano"
    )

    gene_labels_ma <- mod_label_server(
      "label_ma",
      data_list = reactive({
        vol_data()
      }),
      method = "ma"
    )

    vol_plot <- reactive({
      req(vol_data(), input$plot_color_select)

      p <- vol <- plot_volcano(
        data = vol_data()$data,
        plot_colors = plot_colors[[input$plot_color_select]],
        anotate_genes = gene_labels()
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })

    output$volcano_plot <- renderPlot({
      print(vol_plot())
    })

    download_volcano <- ottoPlots::mod_download_figure_server(
      "download_volcano",
      filename = "volcano_plot",
      figure = reactive({
        vol_plot()
      }),
      label = ""
    )

    # ma plot----------------
    ma_plot <- reactive({
      req(vol_data())

      p <- plot_ma(
        data = vol_data()$data,
        plot_colors = plot_colors[[input$plot_color_select]],
        anotate_genes = gene_labels_ma()
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })

    output$ma_plot <- renderPlot({
      print(ma_plot())
    })

    download_ma <- ottoPlots::mod_download_figure_server(
      "download_ma",
      filename = "ma_plot",
      figure = reactive({
        ma_plot()
      }),
      label = ""
    )

    output$scatter_plot <- renderPlot({
      req(!is.null(deg$limma$top_genes))

      p <- plot_deg_scatter(
        select_contrast = input$select_contrast,
        comparisons = deg$limma$comparisons,
        top_genes = deg$limma$top_genes,
        limma_p_val = input$limma_p_val,
        limma_fc = input$limma_fc,
        contrast_samples = contrast_samples(),
        processed_data = pre_process$data(),
        sample_info = pre_process$sample_info()
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })

    # Split up and down genes into two data bases ---------
    up_reg_data <- reactive({
      req(!is.null(heat_data()))

      return(
        heat_data()$genes[heat_data()$bar == 1, ]
      )
    })
    down_reg_data <- reactive({
      req(!is.null(heat_data()))

      return(
        heat_data()$genes[heat_data()$bar == -1, ]
      )
    })

    # enrichment analysis results for both up and down regulated gene
    pathway_deg <- reactive({
      req(!is.null(up_reg_data()))
      withProgress(message = "Enrichment analysis for DEGs.", {
        incProgress(0.1)
        deg_lists <- list()
        lists <- c("Upregulated", "Downregulated")

        for (direction in lists) {
          if (direction == lists[1]) {
            data <- up_reg_data()
          } else {
            data <- down_reg_data()
          }

          gene_names <- merge_data(
            all_gene_names = pre_process$all_gene_names(),
            data = data,
            merge_ID = "ensembl_ID"
          )
          # Only keep the gene names and scrap the data
          deg_lists[[direction]] <- dplyr::select_if(gene_names, is.character)
          incProgress(0.5)
        }
      })
      return(deg_lists)
    })


    enrichment_table_cluster <- mod_11_enrichment_server(
      id = "enrichment_table_cluster",
      gmt_choices = reactive({
        pre_process$gmt_choices()
      }),
      gene_lists = reactive({
        pathway_deg()
      }),
      processed_data = reactive({
        pre_process$data()
      }),
      gene_info = reactive({
        pre_process$all_gene_info()
      }),
      idep_data = idep_data,
      select_org = reactive({
        pre_process$select_org()
      }),
      converted = reactive({
        pre_process$converted()
      }),
      gmt_file = reactive({
        pre_process$gmt_file()
      }),
      plot_grid_lines = reactive({
        pre_process$plot_grid_lines()
      }),
      ggplot2_theme = reactive({
        pre_process$ggplot2_theme()
      })
    )

    list(
      limma = reactive(deg$limma),
      select_factors_model = reactive(input$select_factors_model),
      select_model_comprions = reactive(input$select_model_comprions),
      reference_levels = reactive(factor_reference_levels()),
      counts_deg_method = reactive(input$counts_deg_method)
    )
  })
}

## To be copied in the UI
# mod_06_deg_ui("06_deg_ui_1")

## To be copied in the server
# mod_06_deg_server("06_deg_ui_1")
