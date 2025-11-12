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
    title = "Stats",
    sidebarLayout(
      sidebarPanel(
        style = "height: 90vh; overflow-y: auto;", 
        # Button to run DEG analysis for the specified model
        uiOutput(ns("submit_ui")),
        tags$head(tags$style(
          "#deg-submit_model_button{font-size: 16px;color: red}"
        )),
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
            selected = 3,
            selectize = FALSE
          ),
          tippy::tippy_this(
            ns("counts_deg_method"),
            "Pick the method to detect differentially expressed genes. DESEq2 is recommended for read count data. limma-voom and limma-trend can be used for read counts or normalized expression data.",
            theme = "light"
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
          # Hide FDR cutoff for fold-change only data (type 4)
          conditionalPanel(
            condition = "output.data_file_format != 4",
            column(
              width = 5,
              # Adjusted significant p-value to use
              numericInput(
                inputId = ns("limma_p_val"),
                label = "FDR cutoff:",
                value = 0.1,
                min = 1e-5,
                max = 1,
                step = .05
              ),
              tippy::tippy_this(
                ns("limma_p_val"),
                "Adjusted p-value (FDR) threshold for differentially expressed genes.",
                theme = "light"
              )
            ),
            ns = ns
          ),
          column(
            width = 7,
            # Enforce lower limit on fold change
            tags$head(
              tags$script(HTML("
                $(document).on('shiny:connected', function(event) {
                    $('#deg-limma_fc').on('input', function() {
                        var inputElement = $(this);
                        setTimeout(function() {
                            var value = parseFloat(inputElement.val());
                            if (value < 1) {
                                inputElement.val('1.5');
                                alert('Must be greater than 1.');
                                inputElement.trigger('change'); // Trigger update
                            }
                        }, 1000); // Wait for 1 second
                    });
                });
              "))
            ),
            numericInput(
              inputId = ns("limma_fc"),
              label = "Min fold-change:",
              value = 2,
              min = 1,
              max = 100,
              step = 0.5
            ),
            tippy::tippy_this(
              ns("limma_fc"),
              "Set the minimum fold-change. For example, setting it to 2 identifys genes with at least 2-fold up- or down-regulation. ",
              theme = "light"
            )
          )
        ),
        conditionalPanel(
          condition = "input.counts_deg_method == 3 && output.data_file_format == 1",
          checkboxInput(
            inputId = ns("threshold_wald_test"),
            label = "Threshold-based Wald Test",
            value = FALSE
          ),
          tippy::tippy_this(
            ns("threshold_wald_test"),
            "Use a threshold-based Wald test in DESeq2 to test the hypothesis that the absolute log2 fold change to exceed the chosen cutoff (log2 of the minimum fold-change). Leave unchecked for the standard Wald test against zero.",
            theme = "light"
          ),
          checkboxInput(
            inputId = ns("independent_filtering"),
            label = "Independent filtering of lower counts",
            value = TRUE
          ),
          tippy::tippy_this(
            ns("independent_filtering"),
            "Let DESeq2 drop very low-count genes before p-value adjustment to improve power.",
            theme = "light"
          ),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.step_1 == 'results'",
          selectInput(
            inputId = ns("plot_color_select_1"),
            label = NULL,
            choices = "Red-Green",
            selectize = FALSE
          ),
          tippy::tippy_this(
            ns("plot_color_select_1"),
            "Change the color palette for DEG summary plots.",
            theme = "light"
          ),
          ns = ns
        ),
        tags$br(),
        tags$br(),
        fluidRow(
          column(
            width = 5,
            uiOutput(ns("download_lfc_button")),
            uiOutput(ns("note_download_lfc_button"))
          ),
          column(
            width = 7,
            uiOutput(ns("sig_genes_download_button")),
            uiOutput(ns("note_sig_genes_download"))
          )
        ),
        br(),
        downloadButton(
          outputId = ns("report"),
          label = tags$span(style = "color: red;", "Report")
        ),
        tippy::tippy_this(
          ns("report"),
          "Create an HTML report summarizing the DEG analysis in the Stats tab.",
          theme = "light"
        ),
        br(),
        br(),
        # Tip about uploading design file (only shown when no design file)
        uiOutput(ns("design_file_tip"))
      ),




      mainPanel(
        tabsetPanel(
          id = ns("step_1"),
          tabPanel(
            title = "Experiment Design",
            value = "experiment_design",
            htmlOutput(outputId = ns("list_factors_deg")),
            htmlOutput(outputId = ns("list_block_factors_deg")),
            fluidRow(
              htmlOutput(outputId = ns("select_reference_levels"))
            ),
            htmlOutput(outputId = ns("list_interaction_terms")),
            textOutput(outputId = ns("experiment_design")),
            tags$head(tags$style(
              "#deg-experiment_design{color: red;font-size: 16px;}"
            )),
            br(),
            htmlOutput(outputId = ns("list_model_comparisons"))
          ),
          tabPanel(
            title = "Results",
            value = "results",
            plotOutput(
              outputId = ns("sig_gene_stats")
            ),
            br(),
            ottoPlots::mod_download_figure_ui(ns("download_sig_gene_stats")),
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
            title = "Venn Diagram",
            value = "venn_diagram",
          checkboxInput(
            inputId = ns("up_down_regulated"),
            label = "Split gene lists by up- or down-regulation",
            value = TRUE
          ),
          tippy::tippy_this(
            ns("up_down_regulated"),
            "Separate gene lists into up- and down-regulated sets for venn diagram and upset plot.",
            theme = "light"
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
            value = "r_code",
            shinyAce::aceEditor(
              outputId = ns("deg_code"),
              value = "",
              mode = "r",
              theme = "textmate",
              readOnly = TRUE,
              height = "500px",
              fontSize = 14,
              showLineNumbers = TRUE,
              highlightActiveLine = FALSE,
              showPrintMargin = FALSE
            ),
            br(),
            downloadButton(
              outputId = ns("dl_deg_code"),
              label = "Code"
            ),
            tippy::tippy_this(
              ns("dl_deg_code"),
              "Download the R script used for DEG analysis.",
              theme = "light"
            )
          ),
          tabPanel(
            title = icon("info-circle"),
            includeHTML(app_sys("app/www/help_deg1.html"))
          )
        ),
        tags$script(
          HTML(
            sprintf(
"
Shiny.addCustomMessageHandler('%s', function(message) {
  var tabset = $('#%s');
  if (!tabset.length) { return; }
  var tabLink = tabset.find('a[data-value=\"r_code\"]');
  if (!tabLink.length) { return; }
  var tabItem = tabLink.closest('li');
  if (message.hide) {
    tabItem.hide();
  } else {
    tabItem.show();
  }
});
",
              ns("toggle_stats_r_code"),
              ns("step_1")
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
    title = "DEG",
    sidebarLayout(
      sidebarPanel(
        style = "height: 90vh; overflow-y: auto;",
        h4("Investigate DEGs"),
        br(),
        htmlOutput(outputId = ns("list_comparisons")),
        conditionalPanel("input.step_2 == 'Heatmap'",
          selectInput(
            inputId = ns("heatmap_gene_number"),
            label = "Number of genes:",
            choices = c("All DEGs"),
            selected = "All DEGs",
            selectize = FALSE
          ),
          tippy::tippy_this(
            ns("heatmap_gene_number"),
            "Pick how many genes appear in the DEG heatmap.",
            theme = "light"
          ),
          selectInput(
            inputId = ns("heatmap_fdr_fold"),
            label = "Sort by:",
            choices = c("Fold Change", "FDR"),
            selectize = FALSE
          ),
          tippy::tippy_this(
            ns("heatmap_fdr_fold"),
            "Choose whether genes are ordered by fold-change or FDR.",
            theme = "light"
          ),
          downloadButton(
            outputId = ns("download_heat_data"),
            label = "Heatmap Data"
          ),
          tippy::tippy_this(
            ns("download_heat_data"),
            "Download the data underlying the DEG heatmap.",
            theme = "light"
          ),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.step_2 == 'Volcano Plot' |
            input.step_2 == 'MA Plot' |
            input.step_2 == 'Scatter Plot'",
          selectInput(
            inputId = ns("plot_color_select"),
            label = NULL,
            choices = "Red-Green",
            selectize = FALSE
          ),
          tippy::tippy_this(
            ns("plot_color_select"),
            "Switch the color palette for the selected plot.",
            theme = "light"
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
        conditionalPanel(
          condition = "input.step_2 == 'Scatter Plot' ",
          mod_label_ui(ns("label_scatter")),
          ns = ns
        ),
        br(),
        downloadButton(
          outputId = ns("report_deg"),
          label = tags$span(style = "color: red;", "Report")
        ),
        tippy::tippy_this(
          ns("report_deg"),
          "Create an HTML report summarizing the DEG analysis in the DEG tab.",
          theme = "light"
        )
      ),


      mainPanel(
        tabsetPanel(
          id = ns("step_2"),
          tabPanel(
            title = "Genes",
            tags$div(
              style = "margin-top: 5px;",
              selectInput(
                inputId = ns("gene_direction_filter"),
                label = NULL,
                choices = c("Up-regulated genes", "Down-regulated genes"),
                selected = "Up-regulated genes",
                selectize = FALSE
              )
            ),
            tippy::tippy_this(
              ns("gene_direction_filter"),
              "Filter the gene table to view only upregulated or downregulated genes.",
              theme = "light"
            ),
            DT::dataTableOutput(ns("deg_gene_table")),
            tippy::tippy_this(
              ns("deg_gene_table"),
              "Click a gene row to open its plot. Cuick column headers to sort.",
              theme = "light"
            )
          ),
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
            ),
            ottoPlots::mod_download_figure_ui(ns("download_scatter"))
          ),
          tabPanel(
            title = "Enrichment",
            br(),
            mod_11_enrichment_ui(ns("enrichment_table_cluster")),
          ),
          tabPanel(
            title = icon("info-circle"),
            includeHTML(app_sys("app/www/help_deg2.html"))
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
    summary_only_formats <- c(3, 4)
    is_summary_format <- function(value) {
      if (is.null(value) || length(value) == 0) {
        return(FALSE)
      }
      suppressWarnings(as.numeric(value)) %in% summary_only_formats
    }

    stats_option_notice_initialized <- reactiveVal(FALSE)

    output$submit_ui <- renderUI({
      # req(model_comparisons()) # this is stopping LCF data from getting through Stats
      tagList(
        fluidRow(
          column(
            width = 7,
            h4("Identify DEGs")
          ),
          column(
            width = 5,
            actionButton(
              inputId = ns("submit_model_button"),
              label = "Submit",
              style = "float:right"
            ),
            tippy::tippy_this(
              ns("submit_model_button"),
              "Run the differential expression analysis.",
              theme = "light"
            )
          )
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

    observeEvent(
      pre_process$data_file_format(),
      {
        dfmt <- pre_process$data_file_format()
        hide_r_code <- is_summary_format(dfmt)

        session$sendCustomMessage(
          ns("toggle_stats_r_code"),
          list(hide = hide_r_code)
        )

        if (isTRUE(hide_r_code) && identical(input$step_1, "r_code")) {
          updateTabsetPanel(
            session = session,
            inputId = "step_1",
            selected = "results"
          )
        }
      },
      ignoreNULL = FALSE
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
            strong(list_factors$title),
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
          strong("Select a batch/block factor (optional)"),
          choices = choices,
          selected = NULL
        )
      )
    })

    model_comparisons <- reactive({
      req(pre_process$data() & !is_summary_format(pre_process$data_file_format()))

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
        label = strong(model_comparisons()$title),
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
        label = strong(
          "Interaction terms:"
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
              label = h5(paste0("Reference level for ", x)),
              choices = setNames(
                as.list(paste0(x, ":", select_choices[[x]])),
                select_choices[[x]]
              ),
              selectize = FALSE
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

    # Watch for changes in Stats tab options and notify user to click Submit
    observeEvent(
      list(
        input$counts_deg_method,
        input$limma_p_val,
        input$limma_fc,
        input$threshold_wald_test,
        input$independent_filtering,
        input$select_factors_model,
        input$select_block_factors_model,
        input$select_model_comprions,
        input$select_interactions,
        input$reference_level_factor_1,
        input$reference_level_factor_2,
        input$reference_level_factor_3,
        input$reference_level_factor_4,
        input$reference_level_factor_5,
        input$reference_level_factor_6
      ),
      {
        if (!stats_option_notice_initialized()) {
          stats_option_notice_initialized(TRUE)
          return()
        }

        # Only show notification after Submit has been clicked at least once
        if (is.null(input$submit_model_button) || input$submit_model_button == 0) {
          return()
        }

        if (is.null(tab()) || tab() != "Stats") {
          return()
        }

        showNotification(
          ui = "Click Submit to rerun",
          id = "stats_submit_reminder",
          type = "message",
          duration = 5
        )
      },
      ignoreNULL = FALSE
    )

    deg <- reactiveValues(limma = NULL)

    # Show design file tip in sidebar when no design file uploaded
    output$design_file_tip <- renderUI({
      sample_info <- pre_process$sample_info()

      if (is.null(sample_info)) {
        div(
          style = "background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 4px; padding: 10px; margin-bottom: 10px;",
          tags$strong("Tip:"),
          " An ",
          tags$a(
            href = "https://idepsite.wordpress.com/data-format/",
            target = "_blank",
            "experimental design file"
          ),
          " can be uploaded to build a linear model according to experiment design.",
          tags$br(),
          tags$a(
            href = "http://rpubs.com/ge600/deseq2",
            target = "_blank",
            "More info on DESeq2 experiment design"
          )
        )
      }
    })

    observe({
      method <- input$counts_deg_method
      selected_factors <- input$select_factors_model
      sample_info <- pre_process$sample_info()

      needs_factor_warning <- !is.null(method) && method %in% c(1, 2) &&
        length(selected_factors) > 3

      if (isTRUE(needs_factor_warning)) {
        shiny::showNotification(
          ui = shiny::tags$div(
            shiny::tags$strong("Model warning:"),
            " limma-trend and limma-voom support up to three main factors in iDEP."
          ),
          type = "warning",
          id = "deg1-factor-limit-warning",
          duration = NULL,
          closeButton = TRUE
        )
      } else {
        shiny::removeNotification("deg1-factor-limit-warning")
      }

      needs_rank_warning <- FALSE
      if (!is.null(method) && method %in% c(1, 2) &&
          !is.null(sample_info) && length(selected_factors) > 0) {
        # Ensure sample_info is a data frame (not a matrix)
        if (!is.data.frame(sample_info)) {
          sample_info <- as.data.frame(sample_info, stringsAsFactors = FALSE)
        }

        valid_factors <- intersect(selected_factors, colnames(sample_info))

        if (length(valid_factors) > 0) {
          level_counts <- vapply(
            valid_factors,
            function(factor_col) {
              length(unique(stats::na.omit(sample_info[[factor_col]])))
            },
            integer(1)
          )

          if (any(level_counts == 0)) {
            possible_combos <- 0
          } else {
            possible_combos <- prod(as.numeric(level_counts))
          }
          n_samples <- nrow(sample_info)
          needs_rank_warning <- possible_combos > n_samples

          if (isTRUE(needs_rank_warning)) {
            warning_text <- sprintf(
              " selected factors create up to %d combinations but only %d samples are available. Reduce factors or ensure sufficient replication.",
              possible_combos,
              n_samples
            )

            shiny::showNotification(
              ui = shiny::tags$div(
                shiny::tags$strong("Model may be rank-deficient:"),
                warning_text
              ),
              type = "warning",
              id = "deg1-rank-warning",
              duration = NULL,
              closeButton = TRUE
            )
          } else {
            shiny::removeNotification("deg1-rank-warning")
          }
        } else {
          shiny::removeNotification("deg1-rank-warning")
        }
      } else {
        shiny::removeNotification("deg1-rank-warning")
      }
    })

    warning_type <- reactiveVal(NULL)
    # Observe submit button ------
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
          
          if (is.null(input$select_model_comprions) && 
              !is_summary_format(pre_process$data_file_format())) {
            warning_type("NoComparison")
            deg$limma <- NULL
          } else {
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
              independent_filtering = independent_filtering,
              descr = pre_process$descr()
            )
            # Check for returned errors
            if ("character" %in% class(deg$limma)){
              if (grepl("the model matrix is not full rank", deg$limma)){
                warning_type("FullRankError")
                deg$limma <- NULL
              } else {
                warning_type("Unknown")
                deg$limma <- NULL
              }
            } else {
              updateTabsetPanel(
                session = session,
                inputId = "step_1",
                selected = "results"
              )
            }
          }
        })
      }
    )

    # Dynamic modal for different error types in Stats
    observe({
      req(!is.null(warning_type()))
        modal_title <- switch(
          warning_type(),
          "FullRankError" = "Please Check Experiment Design",
          "NoComparison" = "No Comparison Selected",
          "Unknown Error" = "Analysis Error Occurred"
        )
        modal_text <- switch(
          warning_type(),
          "FullRankError" = paste(
            'The model matrix generated for the analysis is not 
            "full rank." This is often due to a flaw in the experimental 
            design submitted. Check that all combinations of design factors are 
            accounted for and have entries in your data.'
          ),
          "NoComparison" = paste(
          "No comparisons selected to perform Stats analysis on. Please select 
          group comparisons (checkboxes) before submitting again."
          ), 
          "UnknownError" = paste(
            "An unexpected error occurred during analysis. 
             Please check your inputs and try again. If the error persists, 
             contact the admin")
        )
        showModal(
          modalDialog(
            title = modal_title,
            modal_text,
            easyClose = TRUE,
            footer = modalButton("Close"),
            size = "s"
          ),
          session = session
        )
        
        warning_type(NULL)
    })
    
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

    file_name <- reactive({
      dfmt <- pre_process$data_file_format()
      # If normalized data, label as generic "limma"
      if (!is.null(dfmt) && dfmt == 2) {
        method_label <- "limma"
      } else {
        idx <- as.numeric(input$counts_deg_method)
        if (is.na(idx) || is.null(idx) || idx < 1 || idx > length(deg_method)) {
          # Fallback: prefer DESeq2 when raw counts exist, otherwise use limma
          method_label <- if (!is.null(pre_process$raw_counts())) "DESeq2" else "limma"
        } else {
          method_label <- deg_method[idx]
        }
      }
      paste0(method_label, ".csv")
    })

    output$download_lfc <- downloadHandler(
      filename = function() {
        paste0("Results_LFC_Pval_", file_name())
      },
      content = function(file) {
        write.csv(deg_info(), file, row.names = FALSE)
      }
    )

    output$download_lfc_button <- renderUI({
      req(!is.null(deg_info()))
      downloadButton(
        outputId = ns("download_lfc"),
        "Results"
      )
    })

    output$note_download_lfc_button <- renderUI({
      req(!is.null(deg_info()))
      tippy::tippy_this(
        elementId = ns("download_lfc"),
        tooltip = "Download log fold-change, adjusted p-values, and the processed data.",
        theme = "light"
      )
    })

    observe({
      updateSelectInput(
        session = session,
        inputId = "plot_color_select_1",
        choices = plot_choices
      )
    })
    
    sig_genes_p <- reactive({
      req(!is.null(deg$limma$results))
      p <- sig_genes_plot(
        results = deg$limma$results,
        plot_colors = plot_colors[[input$plot_color_select_1]]
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })
    output$sig_gene_stats <- renderPlot({
      req(!is.null(deg$limma$results))
      p <- sig_genes_plot(
        results = deg$limma$results,
        plot_colors = plot_colors[[input$plot_color_select_1]]
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })

    download_sig_gene_stats <- ottoPlots::mod_download_figure_server(
      id = "download_sig_gene_stats",
      filename = "sig_gene_stats",
      figure = reactive({
        sig_genes_p()
      }),
      label = ""
    )

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

    output$sig_genes_download <- downloadHandler(
      filename = function() {
        paste0("Up_and_down_genes_", file_name())
      },
      content = function(file) {
        # Convert the matrix to a data frame and replace values
        list_genes_df <- dplyr::mutate_all(
          as.data.frame(deg$limma$results),
          ~ dplyr::case_when(
            . == -1 ~ "Down",
            . == 1 ~ "Up",
            . == 0 ~ "None"
          )
        )

        # Add rownames as a new column
        list_genes_df$gene_id <- rownames(list_genes_df)

        # Make gene_id the first column
        list_genes_df <- list_genes_df[,
          c("gene_id", setdiff(names(list_genes_df), "gene_id"))
        ]
        
        # Filter out unenriched genes (not up or down regulated)
        list_genes_df <- list_genes_df[apply(
          as.data.frame(list_genes_df[ , -1]), 1, 
          function(row) any(row %in% c("Up", "Down"))
        ), ]
        
        write.csv(list_genes_df, file, row.names = FALSE)
      }
    )

    output$sig_genes_download_button <- renderUI({
      req(!is.null(deg$limma$results))
      downloadButton(
        outputId = ns("sig_genes_download"),
        "Gene Lists"
      )
    })

    output$note_sig_genes_download <- renderUI({
      req(!is.null(deg$limma$results))
      tippy::tippy_this(
        elementId = ns("sig_genes_download"),
        tooltip = "Download the up- and down-regulated gene lists for each comparison.",
        theme = "light"
      )
    })

    # Update the Ace editor with R code when DEG results are available
    observe({
      req(!is.null(deg$limma))
      req(!is.null(deg$limma$expr))
      shinyAce::updateAceEditor(
        session = session,
        editorId = "deg_code",
        value = deg$limma$expr
      )
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
        tagList(
          selectInput(
            inputId = ns("select_comparisons_venn"),
            label = NULL,
            choices = list("All" = "All"),
            selected = "All",
            selectize = FALSE
          ),
          tippy::tippy_this(
            ns("select_comparisons_venn"),
            "Choose which comparisons to include in the Venn diagram.",
            theme = "light"
          )
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
          ),
          tippy::tippy_this(
            ns("select_comparisons_venn"),
            "Select up to five comparisons to visualize in the Venn diagram.",
            theme = "light"
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
        results = venn_data(),
        plots_color_select = load_data$plots_color_select()
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
          label = "Select a comparison:",
          choices = list("All" = "All"),
          selected = "All",
          selectize = FALSE
        )
      } else {
        tagList(
          selectInput(
            inputId = ns("select_contrast"),
            label = "Select a comparison:",
            choices = deg$limma$comparisons,
            selectize = FALSE
          ),
          tippy::tippy_this(
            ns("select_contrast"),
            "Choose a comparison to review its DEGs. For comparison labeled 'Treat vs. Ctrl', positive log fold-change means genes are upregulated in 'Treat'.  Interaction terms start with 'I:', indicating, for example, mutatant specific responses to a treatment.",
            theme = "light"
          ),
          tags$div(
            style = "margin-top: 0.5rem;",
            uiOutput(ns("selected_contrast_counts"))
          ),
          hr(),
          br()
        )
      }
    })

    output$selected_contrast_counts <- renderUI({
      limma_results <- deg$limma$results
      req(!is.null(limma_results))
      req(input$select_contrast)

      limma_results <- as.matrix(limma_results)
      col_names <- colnames(limma_results)
      if (is.null(col_names) || !input$select_contrast %in% col_names) {
        return(NULL)
      }

      calls <- limma_results[, input$select_contrast, drop = TRUE]
      if (length(calls) == 0) {
        return(NULL)
      }

      calls <- suppressWarnings(as.numeric(calls))
      up_count <- sum(calls == 1, na.rm = TRUE)
      down_count <- sum(calls == -1, na.rm = TRUE)

      tags$div(
        class = "deg-regulation-counts",
        tags$p(
          style = "margin-bottom: 0.2rem;",
          sprintf(
            "Up genes: %s; ",
            format(up_count, big.mark = ",", trim = TRUE, scientific = FALSE)
          ),
          sprintf(
            "Down genes: %s",
            format(down_count, big.mark = ",", trim = TRUE, scientific = FALSE)
          )
        )
      )
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
    

    observe({
      req(!is.null(heat_data()$genes))
      
      number_heat_genes <- nrow(heat_data()$genes)
      
      if (number_heat_genes < 5) {
        start <- 1
      } else {
        start <- 5
      }
      heat_number_vec <- seq(from = start, to = number_heat_genes, by = start)
      heat_choices <- c("All DEGs", heat_number_vec)
      updateSelectInput(
        session = session,
        inputId = "heatmap_gene_number",
        choices = heat_choices,
        selected = "All DEGs"
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

    heat_names <- reactive({
      req(!is.null(vol_data()))
      req(!is.null(input$heatmap_fdr_fold))
      req(!is.null(input$heatmap_gene_number))
      req(input$heatmap_gene_number != "All DEGs")
      
      if(input$heatmap_fdr_fold == "FDR"){
        # get ensembl id for number of genes with lowest FDR
        heat_names <- vol_data()$data |>
          dplyr::arrange(FDR) |>
          dplyr::slice(1:input$heatmap_gene_number) |>
          dplyr::pull(Row.names)
      }
      else if(input$heatmap_fdr_fold == "Fold Change"){
        # get ensembl id for number of genes with highest Fold Change
        heat_names <- vol_data()$data |>
          dplyr::arrange(dplyr::desc(abs(Fold))) |>
          dplyr::slice(1:input$heatmap_gene_number) |>
          dplyr::pull(Row.names)
      }
    })
    
    heat_names_to_ensembl <- reactive({
      req(!is.null(heat_names()))
      
      # match ensembl ids and symbols to filter
      df <- pre_process$all_gene_names() |>
        dplyr::filter(symbol %in% heat_names()) |>
        dplyr::select(ensembl_ID)
      df[['ensembl_ID']]
    })
    
    # bar to make the heatmap module reactive
    # otherwise, error when switching heatmap
    heatmap_bar <- reactive({
      req(!is.null(heat_data()))
      heat_data()$bar
    })
    
    deg2_heat_bar <- reactive({
      req(!is.null(heatmap_bar()))
      req(!is.null(input$heatmap_gene_number))
      
      if(input$heatmap_gene_number == "All DEGs"){
        heatmap_bar()
      }
      else{
        req(!is.null(heat_names_to_ensembl()))
        
        heatmap_bar()[names(heatmap_bar()) %in% heat_names_to_ensembl()]
      }
    })
    
    deg2_heat_data <- reactive({
      req(!is.null(heat_data()))
      req(!is.null(input$heatmap_gene_number))
      
      if(input$heatmap_gene_number == "All DEGs"){
        heat_data()$genes
      }
      else{
        req(!is.null(heat_names_to_ensembl()))
        # filter heat_data()$genes and heatmap_bar() for desired genes
        heat_data()$genes[rownames(heat_data()$genes) %in% heat_names_to_ensembl(),]
      }
      
    })
    
    observe({
      req(!is.null(deg2_heat_data()))
      req(!is.null(deg2_heat_bar()))
      
      heatmap_module <- mod_12_heatmap_server(
        id = "12_heatmap_1",
        data = reactive({
          deg2_heat_data()
        }),
        bar = reactive({
          deg2_heat_bar()
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
    })

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
    
    output$download_heat_data <- downloadHandler(
      filename = function() {
        "DEG_Heatmap_Data.csv"
      },
      content = function(file) {
        req(!is.null(deg2_heat_data()))
        df <- deg2_heat_data()
        # Center data to match heatmap
        df <- df - rowMeans(df, na.rm = TRUE)
        # Convert row names to gene symbols
        df <- data.frame(
          Gene_ID = rownames(df),
          rowname_id_swap(
            data_matrix = df,
            all_gene_names = pre_process$all_gene_names(),
            select_gene_id = pre_process$select_gene_id()
          )
        )
        rownames(df) <- gsub(" ", "", rownames(df))  
        
        write.csv(df, file)
      }
    )

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
        contrast_samples = contrast_samples(),
        all_gene_names = pre_process$all_gene_names(),
        select_gene_id = load_data$select_gene_id()
      )
    })

    deg_gene_table_data <- reactive({
      req(!is.null(deg$limma))
      req(!is.null(deg$limma$top_genes))
      req(input$select_contrast, input$limma_p_val, input$limma_fc)

      empty_tbl <- list(
        display = data.frame(
          `ID` = character(0),
          `Ensembl ID` = character(0),
          `log2 FC` = character(0),
          `Adj. Pval` = character(0),
          Description = character(0),
          stringsAsFactors = FALSE
        ),
        meta = data.frame(
          ensembl_ID = character(0),
          symbol = character(0),
          entrezgene_id = character(0),
          log2FC = numeric(0),
          Adjusted_P_value = numeric(0),
          description = character(0),
          stringsAsFactors = FALSE
        ),
        order_dir = "desc"
      )

      top_list <- deg$limma$top_genes
      if (length(top_list) == 0) {
        return(empty_tbl)
      }

      idx <- match(input$select_contrast, names(top_list))
      if (is.na(idx)) {
        idx <- 1
      }

      top_df <- top_list[[idx]]

      if (is.null(top_df) || nrow(top_df) == 0) {
        return(empty_tbl)
      }

      top_df <- as.data.frame(top_df)

      if (ncol(top_df) < 2) {
        return(empty_tbl)
      }

      top_df <- top_df[, 1:2, drop = FALSE]
      colnames(top_df) <- c("log2FC", "Adjusted_P_value")
      top_df$ensembl_ID <- rownames(top_df)

      top_df <- top_df |>
        dplyr::filter(
          is.finite(log2FC),
          is.finite(Adjusted_P_value)
        )

      fc_cutoff <- input$limma_fc
      if (is.null(fc_cutoff) || !is.finite(fc_cutoff) || fc_cutoff <= 0) {
        fc_cutoff <- 1
      }

      top_df <- top_df |>
        dplyr::filter(
          Adjusted_P_value <= input$limma_p_val,
          abs(log2FC) >= log2(fc_cutoff)
        )

      direction <- input$gene_direction_filter
      if (is.null(direction)) {
        direction <- "Up-regulated genes"
      }

      if (identical(direction, "Up-regulated genes")) {
        top_df <- top_df |>
          dplyr::filter(log2FC > 0)
      } else if (identical(direction, "Down-regulated genes")) {
        top_df <- top_df |>
          dplyr::filter(log2FC < 0)
      }

      order_dir <- if (identical(direction, "Down-regulated genes")) "asc" else "desc"

      if (nrow(top_df) == 0) {
        empty_tbl$order_dir <- order_dir
        return(empty_tbl)
      }

      gene_names <- pre_process$all_gene_names()

      if (!is.null(gene_names) && "ensembl_ID" %in% colnames(gene_names)) {
        gene_names <- gene_names |>
          dplyr::distinct(ensembl_ID, .keep_all = TRUE)
        top_df <- top_df |>
          dplyr::left_join(
            gene_names |>
              dplyr::select(ensembl_ID, symbol),
            by = "ensembl_ID"
          )
      } else {
        top_df$symbol <- NA_character_
      }

      gene_info <- pre_process$all_gene_info()
      if (!is.null(gene_info) &&
        "ensembl_gene_id" %in% colnames(gene_info)) {
        keep_cols <- c("ensembl_gene_id", "entrezgene_id", "description")
        keep_cols <- keep_cols[keep_cols %in% colnames(gene_info)]
        gene_info <- gene_info |>
          dplyr::select(dplyr::all_of(keep_cols)) |>
          dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
        top_df <- top_df |>
          dplyr::left_join(
            gene_info,
            by = c("ensembl_ID" = "ensembl_gene_id")
          )
      }

      if (!"entrezgene_id" %in% colnames(top_df)) {
        top_df$entrezgene_id <- NA_character_
      }
      if (!"description" %in% colnames(top_df)) {
        top_df$description <- NA_character_
      }

      top_df <- top_df |>
        dplyr::mutate(
          symbol = dplyr::coalesce(symbol, ensembl_ID),
          entrezgene_id = as.character(entrezgene_id),
          description = as.character(description),
          description = gsub(";.*|\\[.*", "", description),
          description = dplyr::na_if(trimws(description), "")
        )

      if (identical(order_dir, "asc")) {
        top_df <- top_df |>
          dplyr::arrange(
            log2FC,
            Adjusted_P_value
          )
      } else {
        top_df <- top_df |>
          dplyr::arrange(
            dplyr::desc(log2FC),
            Adjusted_P_value
          )
      }

      meta_df <- top_df |>
        dplyr::transmute(
          ensembl_ID,
          symbol,
          entrezgene_id,
          log2FC,
          Adjusted_P_value,
          description
        )

      display_df <- top_df |>
        dplyr::mutate(
          symbol_display = ifelse(
            is.na(symbol),
            "",
            htmltools::htmlEscape(symbol)
          ),
          ensembl_display = ifelse(
            is.na(ensembl_ID),
            "",
            htmltools::htmlEscape(ensembl_ID)
          ),
          description_display = ifelse(
            is.na(description),
            NA_character_,
            htmltools::htmlEscape(description)
          ),
          `ID` = dplyr::if_else(
            !is.na(entrezgene_id) & entrezgene_id != "",
            sprintf(
              "<a href='https://www.ncbi.nlm.nih.gov/gene/%s' target='_blank'>%s</a>",
              entrezgene_id,
              symbol_display
            ),
            symbol_display
          ),
          `Ensembl ID` = dplyr::if_else(
            !is.na(ensembl_ID) & grepl("^ENS", ensembl_ID),
            sprintf(
              "<a href='https://www.ensembl.org/id/%s' target='_blank'>%s</a>",
              ensembl_ID,
              ensembl_display
            ),
            ensembl_display
          ),
          `log2 FC` = sprintf("%.3f", log2FC),
          `Adj. Pval` = {
            formatted <- formatC(Adjusted_P_value, format = "e", digits = 2)
            formatted <- gsub("e([+-])0*(\\d+)", "e\\1\\2", formatted)
            formatted
          },
          Description = dplyr::coalesce(description_display, "")
        ) |>
        dplyr::transmute(
          `ID`,
          `Ensembl ID`,
          `log2 FC`,
          `Adj. Pval`,
          Description
        ) |>
        as.data.frame(stringsAsFactors = FALSE)

      list(
        display = display_df,
        meta = meta_df,
        order_dir = order_dir
      )
    })

    output$deg_gene_table <- DT::renderDataTable({
      req(deg_gene_table_data())
      data <- deg_gene_table_data()
      shiny::validate(
        shiny::need(
          nrow(data$display) > 0,
          "No genes meet the significance thresholds for this comparison."
        )
      )
      order_dir <- data$order_dir
      if (is.null(order_dir)) {
        order_dir <- "desc"
      }
      DT::datatable(
        data$display,
        options = list(
          pageLength = 100,
          lengthChange = FALSE,
          dom = "frtip",
          scrollX = TRUE,
          order = list(list(2, order_dir))
        ),
        selection = list(mode = "single"),
        escape = FALSE,
        class = "cell-border stripe",
        rownames = FALSE
      )
    })

    observeEvent(deg_gene_table_data(), {
      table_data <- deg_gene_table_data()
      req(nrow(table_data$display) > 0)
      selected <- input$deg_gene_table_rows_selected
      if (!is.null(selected) && length(selected) > 0) {
        valid_sel <- selected[selected >= 1 & selected <= nrow(table_data$display)]
        if (!identical(valid_sel, selected)) {
          DT::selectRows(
            DT::dataTableProxy("deg_gene_table", session = session),
            valid_sel
          )
        }
      }
    }, ignoreNULL = FALSE)

    selected_gene_meta <- reactive({
      table_data <- deg_gene_table_data()
      req(nrow(table_data$meta) > 0)
      sel <- input$deg_gene_table_rows_selected
      if (is.null(sel) || length(sel) == 0) {
        return(NULL)
      }
      sel <- sel[1]
      if (sel < 1 || sel > nrow(table_data$meta)) {
        return(NULL)
      }
      table_data$meta[sel, , drop = FALSE]
    })

    gene_expression_data <- reactive({
      meta <- selected_gene_meta()
      if (is.null(meta)) {
        return(NULL)
      }
      expr_matrix <- pre_process$data()
      req(!is.null(expr_matrix))

      gene_id <- meta$ensembl_ID
      if (is.null(gene_id) || is.na(gene_id) || gene_id == "" ||
        !gene_id %in% rownames(expr_matrix)) {
        alt_id <- meta$symbol
        if (!is.null(alt_id) && !is.na(alt_id) && alt_id %in% rownames(expr_matrix)) {
          gene_id <- alt_id
        } else {
          return(NULL)
        }
      }

      sample_idx <- contrast_samples()
      if (is.null(sample_idx) || length(sample_idx) == 0) {
        sample_names <- colnames(expr_matrix)
      } else if (is.numeric(sample_idx)) {
        sample_idx <- sample_idx[sample_idx >= 1 & sample_idx <= ncol(expr_matrix)]
        sample_names <- colnames(expr_matrix)[sample_idx]
      } else {
        sample_names <- intersect(sample_idx, colnames(expr_matrix))
        if (length(sample_names) == 0) {
          sample_names <- colnames(expr_matrix)
        }
      }

      if (!gene_id %in% rownames(expr_matrix) || length(sample_names) == 0) {
        return(NULL)
      }

      expr_values <- as.numeric(expr_matrix[gene_id, sample_names, drop = TRUE])
      if (all(is.na(expr_values))) {
        return(NULL)
      }

      sample_groups <- detect_groups(
        sample_names,
        pre_process$sample_info(),
        preserve_original = TRUE
      )
      sample_groups[is.na(sample_groups) | sample_groups == ""] <- sample_names[is.na(sample_groups) | sample_groups == ""]

      raw_data_values <- rep(NA_real_, length(sample_names))
      raw_data_matrix <- NULL
      if (!is.null(pre_process$raw_counts())) {
        raw_data_matrix <- pre_process$raw_counts()
      } else if (!is.null(load_data$converted_data())) {
        raw_data_matrix <- load_data$converted_data()
      }
      if (!is.null(raw_data_matrix)) {
        if (inherits(raw_data_matrix, "SummarizedExperiment")) {
          if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
            raw_data_matrix <- SummarizedExperiment::assay(raw_data_matrix)
          } else {
            raw_data_matrix <- NULL
          }
        }
      }
      if (!is.null(raw_data_matrix)) {
        if (is.data.frame(raw_data_matrix)) {
          raw_data_matrix <- as.matrix(raw_data_matrix)
        }
        if (is.matrix(raw_data_matrix)) {
          available_samples <- intersect(sample_names, colnames(raw_data_matrix))
          if (length(available_samples) > 0) {
            candidate_values <- as.character(c(meta$ensembl_ID, meta$symbol, gene_id))
            raw_data_candidates <- unique(Filter(
              function(x) !is.null(x) && !is.na(x) && nzchar(x),
              candidate_values
            ))
            for (candidate in raw_data_candidates) {
              if (candidate %in% rownames(raw_data_matrix)) {
                raw_vec <- as.numeric(raw_data_matrix[candidate, available_samples, drop = TRUE])
                match_idx <- match(available_samples, sample_names)
                raw_data_values[match_idx] <- raw_vec
                break
              }
            }
          }
        }
      }

      data.frame(
        sample = sample_names,
        expression = expr_values,
        group = factor(sample_groups, levels = unique(sample_groups)),
        raw_data = raw_data_values,
        stringsAsFactors = FALSE
      )
    })

    selected_gene_plot_data <- reactiveVal(NULL)

    mod_gene_expression_plot_server(
      id = "deg_gene_expression",
      plot_data = reactive(selected_gene_plot_data()),
      palette_name = reactive(load_data$plots_color_select()),
      plot_grid_lines = reactive(pre_process$plot_grid_lines()),
      ggplot2_theme = reactive(pre_process$ggplot2_theme()),
      counts_are_counts = reactive(!is.null(pre_process$raw_counts())),
      download_filename = "gene_expression_plot",
      default_plot_type = "bar",
      default_data_type = "raw"
    )

    genes_tab_notice_shown <- reactiveVal(FALSE)

    observeEvent(list(input$step_2, tab()), {
      req(!is.null(input$step_2))
      if (!identical(tab(), "DEG")) {
        return()
      }
      if (identical(input$step_2, "Genes") && !genes_tab_notice_shown()) {
        showNotification(
          ui = tags$div(
            tags$strong("Tip:"),
            " Click on a row for a gene plot."
          ),
          type = "message",
          duration = 8,
          closeButton = TRUE
        )
        genes_tab_notice_shown(TRUE)
      }
    }, ignoreNULL = TRUE)

    observeEvent(input$deg_gene_table_rows_selected, {
      sel <- input$deg_gene_table_rows_selected
      if (is.null(sel) || length(sel) == 0) {
        selected_gene_plot_data(NULL)
        return()
      }

      expr_df <- gene_expression_data()
      if (is.null(expr_df) || nrow(expr_df) == 0) {
        selected_gene_plot_data(NULL)
        return()
      }

      meta <- selected_gene_meta()
      if (is.null(meta)) {
        selected_gene_plot_data(NULL)
        return()
      }

      symbol_value <- meta$symbol
      if (is.null(symbol_value) || is.na(symbol_value) || symbol_value == "") {
        symbol_value <- meta$ensembl_ID
      }
      description_value <- meta$description
      if (is.null(description_value) || is.na(description_value)) {
        description_value <- ""
      }
      display_name <- symbol_value
      if (nzchar(description_value)) {
        display_name <- paste0(symbol_value, ": ", description_value)
      }

      selected_gene_plot_data(
        list(
          data = expr_df,
          display_name = display_name
        )
      )

      showModal(
        modalDialog(
          mod_gene_expression_plot_ui(
            id = ns("deg_gene_expression"),
            plot_height = "500px",
            show_download = TRUE
          ),
          easyClose = TRUE,
          footer = modalButton("Close"),
          size = "l"
        )
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
    
    gene_labels_scat <- mod_label_server(
      "label_scatter",
      data_list = reactive({
        vol_data()
      }),
      method = "scatter"
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
        anotate_genes = gene_labels_ma()              ### RESTART HERE
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

    ## scatter plot-----------
    scatter_plot <- reactive({
      req(!is.null(deg$limma$top_genes))

      p <- plot_deg_scatter(
        select_contrast = input$select_contrast,
        comparisons = deg$limma$comparisons,
        top_genes = deg$limma$top_genes,
        limma_p_val = input$limma_p_val,
        limma_fc = input$limma_fc,
        contrast_samples = contrast_samples(),
        processed_data = pre_process$data(),
        sample_info = pre_process$sample_info(),
        plot_colors = plot_colors[[input$plot_color_select]],
        all_gene_names = pre_process$all_gene_names(),
        anotate_genes = gene_labels_scat()
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })

    output$scatter_plot <- renderPlot({
      print(scatter_plot())
    })

    download_scatter <- ottoPlots::mod_download_figure_server(
      "download_scatter",
      filename = "scatter_plot",
      figure = reactive({
        scatter_plot()
      }),
      label = ""
    )

    # Split up and down genes into two data bases ---------
    up_reg_data <- reactive({
      req(!is.null(heat_data()))

      return(
        heat_data()$genes[heat_data()$bar == 1, , drop = FALSE]
      )
    })
    down_reg_data <- reactive({
      req(!is.null(heat_data()))

      return(
        heat_data()$genes[heat_data()$bar == -1, , drop = FALSE]
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
          
          if (nrow(data) != 0){
          gene_names <- merge_data(
            all_gene_names = pre_process$all_gene_names(),
            data = data,
            merge_ID = "ensembl_ID"
          )
          deg_lists[[direction]] <- dplyr::select_if(gene_names, is.character)
          } else {
            deg_lists[[direction]] <- NULL
          }
          # Only keep the gene names and scrap the data
          
          incProgress(0.5)
        }
      })
      return(deg_lists)
    })


    # Report generation for Stats tab ----
    output$report <- downloadHandler(
      filename = paste0(
        "deg_stats_report_",
        format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
        ".html"
      ),
      content = function(file) {
        withProgress(message = "Generating Report", {
          incProgress(0.2)

          # Copy the report file to a temporary directory
          tempReport <- file.path(tempdir(), "deg_stats_workflow.Rmd")
          tempReport <- gsub("\\", "/", tempReport, fixed = TRUE)

          markdown_location <- app_sys("app/www/RMD/deg_stats_workflow.Rmd")
          file.copy(from = markdown_location, to = tempReport, overwrite = TRUE)

          # Get plot colors - must match the format in the main app
          plot_colors_list <- list(
            "Green-Red" = c("green", "grey45", "red"),
            "Red-Green" = c("red", "grey45", "green"),
            "Blue-Red" = c("blue", "grey45", "red"),
            "Green-Magenta" = c("green", "grey45", "magenta"),
            "Orange-Blue" = c("orange", "grey45", "blue")
          )
          selected_colors <- plot_colors_list[[input$plot_color_select_1]]
          if (is.null(selected_colors)) {
            selected_colors <- plot_colors_list[["Red-Green"]]
          }

          # Compute the stats table
          stats_table <- if (!is.null(deg$limma)) {
            genes_stat_table(limma = deg$limma)
          } else {
            NULL
          }

          # Get R code if available
          r_code <- if (!is.null(deg$limma) && !is.null(deg$limma$expr)) {
            deg$limma$expr
          } else {
            NULL
          }

          # Get venn diagram data if available
          venn_data_param <- NULL
          venn_comparisons <- NULL
          if (!is.null(deg$limma) && !is.null(input$select_comparisons_venn)) {
            tryCatch({
              venn_data_param <- prep_venn(
                limma = deg$limma,
                up_down_regulated = input$up_down_regulated,
                select_comparisons_venn = input$select_comparisons_venn
              )
              venn_comparisons <- input$select_comparisons_venn
            }, error = function(e) {
              # Silently fail if venn data not available
            })
          }

          # Set up parameters to pass to Rmd document
          params <- list(
            pre_processed_descr = pre_process$descr(),
            mapping_statistics = pre_process$mapping_statistics(),
            sample_info = pre_process$sample_info(),
            deg_results = if (!is.null(deg$limma)) deg$limma$results else NULL,
            stats_table = stats_table,
            limma_p_val = input$limma_p_val,
            limma_fc = input$limma_fc,
            counts_deg_method = input$counts_deg_method,
            data_file_format = pre_process$data_file_format(),
            threshold_wald_test = if (!is.null(input$threshold_wald_test)) input$threshold_wald_test else FALSE,
            independent_filtering = if (!is.null(input$independent_filtering)) input$independent_filtering else TRUE,
            plot_colors = selected_colors,
            all_gene_names = pre_process$all_gene_names(),
            select_gene_id = pre_process$select_gene_id(),
            plot_grid_lines = pre_process$plot_grid_lines(),
            ggplot2_theme = pre_process$ggplot2_theme(),
            plots_color_select = load_data$plots_color_select(),
            r_code = r_code,
            venn_data = venn_data_param,
            venn_comparisons = venn_comparisons,
            up_down_regulated = input$up_down_regulated
          )

          req(params)

          # Render the document
          rmarkdown::render(
            input = tempReport,
            output_file = file,
            params = params,
            envir = new.env(parent = globalenv())
          )
        })
      }
    )

    # Report generation for DEG tab ----
    output$report_deg <- downloadHandler(
      filename = paste0(
        "deg_report_",
        format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
        ".html"
      ),
      content = function(file) {
        withProgress(message = "Generating DEG Report", {
          incProgress(0.2)

          # Copy the report file to a temporary directory
          tempReport <- file.path(tempdir(), "deg_workflow.Rmd")
          tempReport <- gsub("\\", "/", tempReport, fixed = TRUE)

          markdown_location <- app_sys("app/www/RMD/deg_workflow.Rmd")
          file.copy(from = markdown_location, to = tempReport, overwrite = TRUE)

          # Get plot colors - must match the format in the main app
          plot_colors_list <- list(
            "Green-Red" = c("green", "grey45", "red"),
            "Red-Green" = c("red", "grey45", "green"),
            "Blue-Red" = c("blue", "grey45", "red"),
            "Green-Magenta" = c("green", "grey45", "magenta"),
            "Orange-Blue" = c("orange", "grey45", "blue")
          )
          selected_colors <- plot_colors_list[[input$plot_color_select]]
          if (is.null(selected_colors)) {
            selected_colors <- plot_colors_list[["Red-Green"]]
          }

          # Calculate up and down gene counts
          up_count <- 0
          down_count <- 0
          if (!is.null(heat_data()$bar)) {
            up_count <- sum(heat_data()$bar == 1, na.rm = TRUE)
            down_count <- sum(heat_data()$bar == -1, na.rm = TRUE)
          }

          # Generate gene tables for up and down regulated genes
          up_genes_table <- NULL
          down_genes_table <- NULL

          if (!is.null(deg$limma$top_genes) && !is.null(input$select_contrast)) {
            top_list <- deg$limma$top_genes
            idx <- match(input$select_contrast, names(top_list))
            if (!is.na(idx) && !is.null(top_list[[idx]])) {
              top_df <- as.data.frame(top_list[[idx]])
              if (ncol(top_df) >= 2) {
                colnames(top_df)[1:2] <- c("log2FC", "Adjusted_P_value")
                top_df$ensembl_ID <- rownames(top_df)

                # Filter by significance
                fc_cutoff <- input$limma_fc
                if (is.null(fc_cutoff) || !is.finite(fc_cutoff) || fc_cutoff <= 0) {
                  fc_cutoff <- 1
                }

                top_df <- top_df |>
                  dplyr::filter(
                    is.finite(log2FC),
                    is.finite(Adjusted_P_value),
                    Adjusted_P_value <= input$limma_p_val,
                    abs(log2FC) >= log2(fc_cutoff)
                  )

                # Add gene names and info
                gene_names <- pre_process$all_gene_names()
                if (!is.null(gene_names) && "ensembl_ID" %in% colnames(gene_names)) {
                  gene_names <- gene_names |>
                    dplyr::distinct(ensembl_ID, .keep_all = TRUE)
                  top_df <- top_df |>
                    dplyr::left_join(
                      gene_names |> dplyr::select(ensembl_ID, symbol),
                      by = "ensembl_ID"
                    )
                } else {
                  top_df$symbol <- NA_character_
                }

                gene_info <- pre_process$all_gene_info()
                if (!is.null(gene_info) && "ensembl_gene_id" %in% colnames(gene_info)) {
                  keep_cols <- c("ensembl_gene_id", "entrezgene_id", "description")
                  keep_cols <- keep_cols[keep_cols %in% colnames(gene_info)]
                  gene_info <- gene_info |>
                    dplyr::select(dplyr::all_of(keep_cols)) |>
                    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
                  top_df <- top_df |>
                    dplyr::left_join(
                      gene_info,
                      by = c("ensembl_ID" = "ensembl_gene_id")
                    )
                }

                if (!"entrezgene_id" %in% colnames(top_df)) {
                  top_df$entrezgene_id <- NA_character_
                }
                if (!"description" %in% colnames(top_df)) {
                  top_df$description <- NA_character_
                }

                top_df <- top_df |>
                  dplyr::mutate(
                    symbol = dplyr::coalesce(symbol, ensembl_ID),
                    entrezgene_id = as.character(entrezgene_id),
                    description = as.character(description),
                    description = gsub(";.*|\\[.*", "", description),
                    description = dplyr::na_if(trimws(description), "")
                  )

                # Split into up and down
                up_df <- top_df |>
                  dplyr::filter(log2FC > 0) |>
                  dplyr::arrange(dplyr::desc(log2FC), Adjusted_P_value) |>
                  dplyr::slice(1:50)

                down_df <- top_df |>
                  dplyr::filter(log2FC < 0) |>
                  dplyr::arrange(log2FC, Adjusted_P_value) |>
                  dplyr::slice(1:50)

                # Create display tables
                if (nrow(up_df) > 0) {
                  up_genes_table <- up_df |>
                    dplyr::transmute(
                      Symbol = symbol,
                      `Ensembl ID` = ensembl_ID,
                      `log2 FC` = sprintf("%.3f", log2FC),
                      `Adj. P-value` = formatC(Adjusted_P_value, format = "e", digits = 2),
                      Description = dplyr::coalesce(description, "")
                    ) |>
                    as.data.frame(stringsAsFactors = FALSE)
                }

                if (nrow(down_df) > 0) {
                  down_genes_table <- down_df |>
                    dplyr::transmute(
                      Symbol = symbol,
                      `Ensembl ID` = ensembl_ID,
                      `log2 FC` = sprintf("%.3f", log2FC),
                      `Adj. P-value` = formatC(Adjusted_P_value, format = "e", digits = 2),
                      Description = dplyr::coalesce(description, "")
                    ) |>
                    as.data.frame(stringsAsFactors = FALSE)
                }
              }
            }
          }

          # Set up parameters to pass to Rmd document
          params <- list(
            pre_processed_descr = pre_process$descr(),
            mapping_statistics = pre_process$mapping_statistics(),
            sample_info = pre_process$sample_info(),
            select_contrast = input$select_contrast,
            limma_p_val = input$limma_p_val,
            limma_fc = input$limma_fc,
            counts_deg_method = input$counts_deg_method,
            data_file_format = pre_process$data_file_format(),
            threshold_wald_test = if (!is.null(input$threshold_wald_test)) input$threshold_wald_test else FALSE,
            independent_filtering = if (!is.null(input$independent_filtering)) input$independent_filtering else TRUE,
            plot_colors = selected_colors,
            all_gene_names = pre_process$all_gene_names(),
            select_gene_id = pre_process$select_gene_id(),
            plot_grid_lines = pre_process$plot_grid_lines(),
            ggplot2_theme = pre_process$ggplot2_theme(),
            plots_color_select = load_data$plots_color_select(),
            heatmap_color_select = pre_process$heatmap_color_select(),
            vol_data = if (!is.null(vol_data())) vol_data() else NULL,
            heat_data = if (!is.null(heat_data()$genes)) heat_data()$genes else NULL,
            heat_bar = if (!is.null(heat_data()$bar)) heat_data()$bar else NULL,
            up_count = up_count,
            down_count = down_count,
            up_genes_table = up_genes_table,
            down_genes_table = down_genes_table,
            gene_labels_volcano = if (!is.null(gene_labels())) gene_labels() else NULL,
            gene_labels_ma = if (!is.null(gene_labels_ma())) gene_labels_ma() else NULL,
            gene_labels_scatter = if (!is.null(gene_labels_scat())) gene_labels_scat() else NULL,
            contrast_samples = if (!is.null(contrast_samples())) contrast_samples() else NULL,
            processed_data = pre_process$data(),
            top_genes = if (!is.null(deg$limma$top_genes)) deg$limma$top_genes else NULL,
            comparisons = if (!is.null(deg$limma$comparisons)) deg$limma$comparisons else NULL
          )

          req(params)

          # Render the document
          rmarkdown::render(
            input = tempReport,
            output_file = file,
            params = params,
            envir = new.env(parent = globalenv())
          )
        })
      }
    )

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
      filter_size = reactive({
        pre_process$filter_size()
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
      }),
      heat_colors = reactive({
        strsplit(load_data$heatmap_color_select(), "-")[[1]][c(1,3)]
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
