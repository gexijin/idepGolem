#' 03_heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_03_clustering_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "Cluster",
    # JavaScript to reset brush
    tags$script(HTML("
      Shiny.addCustomMessageHandler('resetBrush', function(message) {
        Shiny.setInputValue(message.brushId, null);
      });
    ")),
    # Change the style of radio labels
    # Note that the name https://groups.google.com/g/shiny-discuss/c/ugNEaHizlck
    # input IDs should be defined by namespace
    tags$style(
      type = "text/css",
      paste(
        paste0("#", ns("cluster_meth"), " .radio label { font-weight: bold; color: red;}"),
        ".more-options { display: block; width: 100%; }",
        ".more-options summary { display: flex; align-items: center; cursor: pointer; font-weight: 600; margin: 0; }",
        ".more-options summary::marker, .more-options summary::-webkit-details-marker { display: none; }",
        ".more-options summary::before { content: '+'; margin-right: 6px; font-size: 14px; line-height: 1; }",
        ".more-options[open] summary::before { content: '\\2212'; }",
        ".more-options-body { margin-top: 10px; padding: 10px 0; background-color: #f7f9fc; border-radius: 0; width: 100%; }",
        sep = "\n"
      )
    ),
    sidebarLayout(

      # Heatmap Panel Sidebar ----------
      sidebarPanel(
        width = 3,
        # Select Clustering Method ----------
        conditionalPanel(
          condition = "input.cluster_panels == 'Heatmap' |
            input.cluster_panels == 'sample_tab'",
          radioButtons(
            inputId = ns("cluster_meth"),
            label = NULL,
            choices = list(
              "Hierarchical Clustering" = 1,
              "k-Means Clustering" = 2
            ),
            selected = 1
          ),
          ns = ns
        ),
        HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />'),
        conditionalPanel(
          condition = "input.cluster_panels == 'Heatmap' |
          input.cluster_panels == 'word_cloud' |
          input.cluster_panels == 'Gene SD Distribution' ",
          fluidRow(
            column(width = 6, p("Top Genes:")),
            column(
              width = 6,
              numericInput(
                inputId = ns("n_genes"),
                label = NULL,
                min = 10,
                max = 12000,
                value = 1500,
                step = 100
              ),
              tippy::tippy_this(
                ns("n_genes"),
                "Using the transformed data, genes are ranked by standard deviation across all samples. This highlights overall expression changes without using sample groups, enabling exploratory analysis. Normally 1000-5000.",
                theme = "light"
              )
            )
          ),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.cluster_panels == 'Heatmap'",
          fluidRow(
            column(width = 4, p("Dendrogram")),
            column(
              width = 8,
              selectInput(
                inputId = ns("dendrogram_display"),
                label = NULL,
                choices = c(
                  "Row" = "row",
                  "Column" = "column",
                  "Both" = "both",
                  "None" = "none"
                ),
                selected = "row",
                selectize = FALSE
              ),
              tippy::tippy_this(
                ns("dendrogram_display"),
                "Choose whether to display row, column, or both dendrograms.",
                theme = "light"
              )
            )
          ),
          ns = ns
        ),
        conditionalPanel(
          condition = "(input.cluster_panels == 'Heatmap' |
          input.cluster_panels == 'word_cloud' |
          input.cluster_panels == 'sample_tab') &&  input.cluster_meth == 2",

          # k- means slidebar -----------

          sliderInput(
            inputId = ns("k_clusters"),
            label = "Number of Clusters:",
            min = 2,
            max = 20,
            value = 6,
            step = 1
          ),
          tippy::tippy_this(
            ns("k_clusters"),
            "Set how many clusters for k-means.",
            theme = "light"
          ),

          # Re-run k-means with a different seed
          actionButton(
            inputId = ns("k_means_re_run"),
            label = "New Seed"
          ),
          tippy::tippy_this(
            ns("k_means_re_run"),
            "Re-run k-means with a new random seed.",
            theme = "light"
          ),
          # Elbow plot pop-up
          actionButton(
            inputId = ns("elbow_pop_up"),
            label = "How many clusters?"
          ),
          tippy::tippy_this(
            ns("elbow_pop_up"),
            "Show the k-means elbow plot to help choose the number of clusters.",
            theme = "light"
          ),
          # Line break ---------
          HTML(
            '<hr style="height:1px;border:none;
           color:#333;background-color:#333;" />'
          ),
          ns = ns
        ),

        # Clustering methods for Heatmap ----------
        conditionalPanel(
          condition = "input.cluster_meth == 1 &&
            (input.cluster_panels == 'Heatmap' |
            input.cluster_panels == 'word_cloud' |
            input.cluster_panels == 'sample_tab')",
          fluidRow(
            column(width = 4, p("Distance")),
            column(
              width = 8,
              selectInput(
                inputId = ns("dist_function"),
                label = NULL,
                choices = NULL,
                width = "100%",
                selectize = FALSE
              ),
              tippy::tippy_this(
                ns("dist_function"),
                "Pick the distance metric for hierarchical clustering.",
                theme = "light"
              )
            )
          ),
          fluidRow(
            column(width = 4, p("Linkage")),
            column(
              width = 8,
              selectInput(
                inputId = ns("hclust_function"),
                label = NULL,
                choices = c(
                  "average", "complete", "single",
                  "median", "centroid", "mcquitty"
                ),
                width = "100%",
                selectize = FALSE
              ),
              tippy::tippy_this(
                ns("hclust_function"),
                "Choose how clusters are merged when building the dendrogram. Average linkage uses the average distance between all pairs of items in two clusters. Complete linkage uses the maximum distance between items in the two clusters, while single linkage uses the minimum distance. Ward's method minimizes the total within-cluster variance. Guidance: choose 'single' for chaining-sensitive shapes, 'complete' for compact/furthest linkage, 'average' for a balance between extremes, and 'ward.D'/'ward.D2' when you want variance-based, spherical clusters (requires a Euclidean-like distance).",
                theme = "light"
              )
            )
          ),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.cluster_panels == 'Heatmap'",
          fluidRow(
            column(width = 4, p("Samples color")),
            column(
              width = 8,
              htmlOutput(ns("list_factors_heatmap"))
            )
          ),
          fluidRow(
            column(width = 4, p("Label Genes:")),
            column(
              width = 8,
              div(
                id = ns("selected_genes_wrapper"),
                selectizeInput(
                  inputId = ns("selected_genes"),
                  label = NULL,
                  choices = c("Top 5", "Top 10", "Top 15"),
                  multiple = TRUE
                ),
                tippy::tippy_this(
                  ns("selected_genes_wrapper"),
                  "Add labels for the top genes on the heatmap.",
                  theme = "light"
                )
              )
            )
          ),
          tags$details(
            class = "more-options",
            tags$summary("More options"),
            div(
              class = "more-options-body",
              checkboxInput(
                inputId = ns("gene_centering"),
                label = "Center genes (substract mean)",
                value = TRUE
              ),
              tippy::tippy_this(
                ns("gene_centering"),
                "Subtract each gene's mean value before clustering.",
                theme = "light"
              ),
              checkboxInput(
                inputId = ns("gene_normalize"),
                label = "Normalize genes (divide by SD)",
                value = FALSE
              ),
              tippy::tippy_this(
                ns("gene_normalize"),
                "Substract mean and scale by standard deviation before clustering.",
                theme = "light"
              ),
              checkboxInput(
                inputId = ns("letter_overlay"),
                label = "Letter Overlay",
                value = TRUE
              ),
              tippy::tippy_this(
                ns("letter_overlay"),
                "Toggle letter overlays on the sample color bar above the heatmap.",
                theme = "light"
              ),
              checkboxInput(
                inputId = ns("show_color_key"),
                label = "Color Key",
                value = FALSE
              ),
              tippy::tippy_this(
                ns("show_color_key"),
                "Show or hide the color key legend on the main heatmap.",
                theme = "light"
              ),

              fluidRow(
                column(width = 4, p("Sample Colors")),
                column(
                  width = 8,
                  selectInput(
                    inputId = ns("sample_color"),
                    label = NULL,
                    choices = c("Pastel 1", "Dark 2", "Dark 3", 
                                "Set 2", "Set 3", "Warm",
                                "Cold", "Harmonic", "Dynamic"),
                    selected = "Dynamic",
                    selectize = FALSE
                  ),
                  tippy::tippy_this(
                    ns("sample_color"),
                    "Color palette for sample annotations.",
                    theme = "light"
                  )
                )
              ),
              fluidRow(
                column(width = 4, p("Max Z score")),
                column(
                  width = 8,
                  numericInput(
                    inputId = ns("heatmap_cutoff"),
                    label = NULL,
                    value = 3,
                    min = 2,
                    step = 1
                  ),
                  tippy::tippy_this(
                    ns("heatmap_cutoff"),
                    "Cap extremely large or small values on the heatmap. This improves color contrast for most values. Normally 2-4. ",
                    theme = "light"
                  )
                )
              ),
              tags$hr()
            )
          ),
          ns = ns
        ),

        conditionalPanel(
          condition = "input.cluster_panels == 'word_cloud'",
          uiOutput(
            outputId = ns("cloud_ui")
          ),
          ns = ns
        ),
        br(),
        fluidRow(
          conditionalPanel(
            condition = "input.cluster_panels == 'Heatmap' ",
            column(
              width = 4,
              downloadButton(
                outputId = ns("download_heatmap_data"),
                label = "Data"
              ),
              tippy::tippy_this(
                ns("download_heatmap_data"),
                "Download the data currently displayed in the heatmap.",
                theme = "light"
              )
            ),
            ns = ns
          ),
          column(
            width = 4,
            downloadButton(
              outputId = ns("report"),
              label = tags$span(style = "color: red;", "Report")
            ),
            tippy::tippy_this(
              ns("report"),
              "Create an HTML report summarizing the Clustering tab.",
              theme = "light"
            )
          )
        )
      ),




      #########################################################################
      # Main Panel
      #########################################################################

      mainPanel(
        width = 9,
        tabsetPanel(
          id = ns("cluster_panels"),
          # Heatmap panel ----------
          tabPanel(
            title = "Heatmap",
            br(),
            fluidRow(
              column(
                width = 5,
                plotOutput(
                  outputId = ns("heatmap_main"),
                  height = "600px",
                  width = "100%",
                  brush = brushOpts(id = ns("ht_brush"),
                                    delayType = "debounce",
                                    clip = TRUE)
                ),
                tippy::tippy_this(
                  ns("heatmap_main"),
                  "Drag over any region of the heatmap to zoom in.",
                  theme = "light"
                ),
                uiOutput(ns("dl_heatmap_main_download_ui"))
              ),
              column(
                width = 7,
                div(
                  style = "max-height: calc(120vh - 300px); overflow-y: auto; -webkit-overflow-scrolling: touch; padding-right: 10px;",
                  conditionalPanel(
                    condition = paste0(
                      "input.cluster_meth == 2 || ",
                      "(input.cluster_meth == 1 && input.ht_brush != null)"
                    ),
                    checkboxInput(
                      inputId = ns("cluster_enrichment"),
                      label = strong("Show enrichment"),
                      value = FALSE
                    ),
                    tippy::tippy_this(
                      ns("cluster_enrichment"),
                      "Run GO enrichment for the selected genes (hierarchical clustering). For k-means, all clusters are analyzed.",
                      theme = "light"
                    ),
                    ns = ns
                  ),
                  conditionalPanel(
                    condition = "input.cluster_enrichment == 1 ",
                    mod_11_enrichment_ui(ns("enrichment_table_cluster")),
                    ns = ns
                  ),
                  plotOutput(
                    outputId = ns("sub_heatmap"),
                    height = "100%",
                    width = "100%"
                  ),
                  uiOutput(ns("dl_heatmap_sub_download_ui"))
                )
              )
            )
          ),
          tabPanel(
            br(),
            div('Generate a word cloud of pathways that contain genes from the 
                selected cluster (Must run clustering with heatmap first). 
                Words are ranked by frequency.'),
            uiOutput(
              outputId = ns("cloud_error")
            ),
            title = "Word Cloud",
            value = "word_cloud",
            wordcloud2::wordcloud2Output(
              outputId = ns("word_cloud"),
              height = "600px"
            ),
            downloadButton(
              outputId = ns("cloud_download"),
              label = "Data Download"
            )
          ),
          # Gene Standard Deviation Distribution ----------
          tabPanel(
            title = "Gene SD Distribution",
            br(),
            plotOutput(
              outputId = ns("sd_density_plot"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(ns("dl_gene_dist"))
          ),
          # Sample Tree -----------------
          tabPanel(
            title = "Sample Tree",
            value = "sample_tab",
            h5(
              "Using genes with maximum expression level at the top 75%.
               Data is transformed and clustered as specified in the sidebar."
            ),
            br(),
            plotOutput(
              outputId = ns("sample_tree"),
              width = "100%",
              height = "400px"
            ),
            ottoPlots::mod_download_figure_ui(ns("dl_sample_tree"))
          ),
          tabPanel(
            title = icon("info-circle"),
            includeHTML(app_sys("app/www/help_clustering.html"))
          )
        )
      )
    )
  )
}








#########################################################################
# Server function
#########################################################################

#' 03_heatmap Server Functions
#'
#' @noRd
mod_03_clustering_server <- function(id, pre_process, load_data, idep_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Interactive heatmap environment
    shiny_env <- new.env()

    observeEvent(pre_process$data(), {
      data_mat <- pre_process$data()
      req(!is.null(data_mat))
      sample_count <- ncol(data_mat)
      if (!is.null(sample_count) && sample_count > 30 && isTRUE(input$letter_overlay)) {
        updateCheckboxInput(
          session = session,
          inputId = "letter_overlay",
          value = FALSE
        )
      }
    }, ignoreNULL = TRUE)

    # Reset to Heatmap whenever the Cluster tab becomes active again
    observeEvent(tab(), {
      req(tab())
      if (tab() == "Cluster") {
        updateTabsetPanel(
          session = session,
          inputId = "cluster_panels",
          selected = "Heatmap"
        )
      }
    })


    # Update Slider Input ---------
    observe({
      req(tab() == "Cluster")
      req(!is.null(pre_process$data()))
      if (nrow(pre_process$data()) > 12000) {
        max_genes <- 12000
      } else {
        max_genes <- round(nrow(pre_process$data()), -2)
      }
      updateNumericInput(
        inputId = "n_genes",
        max = max_genes
      )
    })
    # Distance functions -----------
    dist_funs <- dist_functions()
    dist_choices <- setNames(
      1:length(dist_funs),
      names(dist_funs)
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "dist_function",
        choices = dist_choices
      )
    })

    dendrogram_selection <- reactive({
      selection <- req(input$dendrogram_display)
      list(
        sample = selection %in% c("column", "both"),
        row = selection %in% c("row", "both")
      )
    })

    # Hclust Functions ----------
    hclust_funs <- hcluster_functions()

    # Sample color bar selector ----------
    output$list_factors_heatmap <- renderUI({
      data_mat <- pre_process$data()
      sample_info <- pre_process$sample_info()
      sample_count <- NULL

      if (!is.null(data_mat)) {
        sample_count <- ncol(data_mat)
      } else if (!is.null(sample_info)) {
        sample_count <- nrow(sample_info)
      }

      choices <- "None"

      include_names <- TRUE
      if (!is.null(sample_count) && sample_count > 0 && !is.null(data_mat) && !is.null(colnames(data_mat))) {
        detected_groups <- tryCatch(
          detect_groups(colnames(data_mat)),
          error = function(e) NULL
        )
        if (!is.null(detected_groups)) {
          unique_groups <- unique(detected_groups)
          include_names <- length(unique_groups) < sample_count
        }
      }

      if (include_names) {
        choices <- c(choices, "Names")
      }

      valid_factors <- character()
      if (!is.null(sample_info) && ncol(sample_info) > 0 && !is.null(sample_count) && sample_count > 0) {
        factor_names <- colnames(sample_info)
        keep_factor <- vapply(
          factor_names,
          function(factor_name) {
            values <- sample_info[, factor_name] # each column
            values <- values[!is.na(values)]
            if(length(unique(values)) >= 20) {
              showNotification(
                paste0("Factor '", factor_name, "' has too many unique values and will be ignored."),
                type = "warning",
                duration = 5
              )
              return(FALSE)
            }
            # ignore if every value is unique or too many unique
            length(unique(values)) < sample_count && length(unique(values)) >= 20
          },
          logical(1)
        )
        valid_factors <- factor_names[keep_factor]
      }

      if (length(valid_factors) > 0) {
        choices <- c(choices, valid_factors)
      }

      if (!is.null(sample_info) && ncol(sample_info) > 0) {
        choices <- c(choices, "All factors")
      }

      # Ensure choices are unique while preserving order
      choices <- unique(choices)

      default_selection <- if ("All factors" %in% choices) {
        "All factors"
      } else if ("Names" %in% choices) {
        "Names"
      } else {
        "None"
      }

      selected <- default_selection
      existing_selection <- isolate(input$select_factors_heatmap)
      if (!is.null(existing_selection) && !is.na(existing_selection) && existing_selection %in% choices) {
        selected <- existing_selection
      }

      tagList(
        selectInput(
          inputId = ns("select_factors_heatmap"),
          label = NULL,
          choices = choices,
          selected = selected,
          selectize = FALSE
        ),
        tippy::tippy_this(
          ns("select_factors_heatmap"),
          "Color bar for sample groups. Choose which sample annotation to display on the top of heatmap.",
          theme = "light"
        )
      )
    })


    # Standard Deviation Density Plot ----------
    sd_density_plot <- reactive({
      req(!is.null(pre_process$data()))

      p <- sd_density(
        data = pre_process$data(),
        n_genes_max = input$n_genes
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })

    output$sd_density_plot <- renderPlot({
      print(sd_density_plot())
    })

    dl_gene_dist <- ottoPlots::mod_download_figure_server(
      id = "dl_gene_dist",
      filename = "sd_density_plot",
      figure = reactive({
        sd_density_plot()
      }),
      label = ""
    )



    # Heatmap Data -----------
    heatmap_data <- reactive({
      req(!is.null(pre_process$data()))

      process_heatmap_data(
        data = pre_process$data(),
        n_genes_max = input$n_genes,
        gene_centering = input$gene_centering,
        gene_normalize = input$gene_normalize,
        sample_centering = FALSE,
        sample_normalize = FALSE,
        all_gene_names = pre_process$all_gene_names(),
        select_gene_id = pre_process$select_gene_id()
      )
    })


    # Heatmap Click Value ---------
    observe({
      req(!is.null(pre_process$all_gene_names()))
      req(!is.null(pre_process$data()))
      req(!is.null(heatmap_data()))

      updateSelectizeInput(
        session = session,
        inputId = "selected_genes",
        label = NULL,
        choices = c("Top 5", "Top 10", "Top 15", row.names(heatmap_data())),
        selected = input$selected_genes
      )
    })

    # split "green-white-red" to c("green", "white", "red")
    heatmap_color_select <- reactive({
      req(pre_process$heatmap_color_select())
      unlist(strsplit(pre_process$heatmap_color_select(), "-"))
    })

    # Cache row/column dendrograms so we only rebuild when inputs that affect
    # clustering actually change.
    row_dendrogram <- eventReactive(
      list(
        heatmap_data(),
        input$cluster_meth,
        input$dist_function,
        input$hclust_function,
        dendrogram_selection()$row
      ),
      {
        if (isolate(input$cluster_meth) != 1 ||
            !isolate(dendrogram_selection()$row)) {
          return(NULL)
        }
        mat <- isolate(heatmap_data())
        if (is.null(mat) || nrow(mat) < 2) {
          return(NULL)
        }
        dist_fun <- dist_funs[[as.numeric(isolate(input$dist_function))]]
        h_fun <- hclust_funs[[isolate(input$hclust_function)]]
        stats::as.dendrogram(h_fun(dist_fun(mat)))
      },
      ignoreNULL = FALSE
    )

    column_dendrogram <- eventReactive(
      list(
        heatmap_data(),
        input$cluster_meth,
        input$dist_function,
        input$hclust_function,
        dendrogram_selection()$sample
      ),
      {
        if (isolate(input$cluster_meth) != 1 ||
            !isolate(dendrogram_selection()$sample)) {
          return(NULL)
        }
        mat <- isolate(heatmap_data())
        if (is.null(mat) || ncol(mat) < 2) {
          return(NULL)
        }
        dist_fun <- dist_funs[[as.numeric(isolate(input$dist_function))]]
        h_fun <- hclust_funs[[isolate(input$hclust_function)]]
        stats::as.dendrogram(h_fun(dist_fun(t(mat))))
      },
      ignoreNULL = FALSE
    )

    # HEATMAP -----------
    # Information on interactivity
    # https://jokergoo.github.io/2020/05/15/interactive-complexheatmap/
    output$heatmap_main <- renderPlot(
      {
        req(!is.null(heatmap_data()))
        req(!is.null(input$select_factors_heatmap))

        withProgress(message = "Creating Heatmap", value = 0, {
          incProgress(0.3, detail = "Processing data")

          shiny_env$ht <- heatmap_main_object()

          incProgress(0.4, detail = "Calculating positions")

          # Ensure heatmap object is valid before getting positions
          if (!is.null(shiny_env$ht)) {
            tryCatch({
              shiny_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht)
            }, error = function(e) {
              # If position detection fails, set to NULL and continue
              shiny_env$ht_pos_main <- NULL
            })
          }

          incProgress(0.3, detail = "Finalizing")
        })

        return(shiny_env$ht)
      }
      #,width = 300 # , # this avoids the heatmap being redraw # no longer needed when removed clicking
      #, height = 600
    )
    
    # Color palette for experiment groups on heatmap
    group_pal <- reactive({
      req(!is.null(pre_process$sample_info()))
      req(!is.na(input$sample_color))

      groups <- as.vector(as.matrix(pre_process$sample_info()))
      pal <- setNames(
        colorspace::qualitative_hcl(length(unique(groups)),
                                    palette = input$sample_color,
                                    c = 70),
        unique(groups)
      )
      sample_list <- as.list(as.data.frame(pre_process$sample_info()))

      lapply(sample_list, function(x){
        setNames(
          pal[unique(x)],
          unique(x)
        )
      })
    })

    main_heatmap_group_info <- reactive({
      req(!is.null(heatmap_data()))
      req(!is.null(input$select_factors_heatmap))
      req(input$select_factors_heatmap != "")
      req(!is.null(input$letter_overlay))
      req(!is.null(input$sample_color))

      if (identical(input$select_factors_heatmap, "None")) {
        groups <- rep("None", ncol(heatmap_data()))
        return(list(
          groups = groups,
          group_colors = NULL
        ))
      }

      group_pal_val <- NULL
      if (input$select_factors_heatmap == "All factors") {
        req(!is.null(pre_process$sample_info()))
        req(!is.null(group_pal()))
        group_pal_val <- group_pal()
      }

      info <- sub_heat_ann(
        data = heatmap_data(),
        sample_info = pre_process$sample_info(),
        select_factors_heatmap = input$select_factors_heatmap,
        group_pal = group_pal_val,
        sample_color = input$sample_color,
        use_letter_overlay = isTRUE(input$letter_overlay)
      )

      list(
        groups = info$groups,
        group_colors = info$group_colors
      )
    })
    
    # Reactive for heatmap generation with double-render prevention
    heatmap_main_object <- reactive({
      req(!is.null(heatmap_data()))
      req(!is.null(input$select_factors_heatmap))
      req(input$select_factors_heatmap != "")
      req(!is.null(input$letter_overlay))

      # Wait for enrichment results if enrichment is enabled for k-means
      if (isTRUE(input$cluster_enrichment) && input$cluster_meth == 2) {
        req(!is.null(cluster_pathway_labels()))
      }

      # Ensure stable state before proceeding
      if (!is.null(pre_process$sample_info())) {
        # For "All factors", ensure group_pal is ready
        if (input$select_factors_heatmap == "All factors") {
          req(!is.null(group_pal()))
        }
      }

      # Assign heatmap to be used in multiple components
      try(
        obj <- heatmap_main(
          data = heatmap_data(),
          cluster_meth = input$cluster_meth,
          heatmap_cutoff = input$heatmap_cutoff,
          sample_info = pre_process$sample_info(),
          select_factors_heatmap = input$select_factors_heatmap,
          dist_funs = dist_funs,
          dist_function = input$dist_function,
          hclust_function = input$hclust_function,
          sample_clustering = dendrogram_selection()$sample,
          heatmap_color_select = heatmap_color_select(),
          row_dend = dendrogram_selection()$row,
          k_clusters = input$k_clusters,
          re_run = input$k_means_re_run,
          selected_genes = input$selected_genes,
          group_pal = if (input$select_factors_heatmap == "All factors") {
            group_pal()
          } else {
            NULL
          },
          sample_color = input$sample_color,
          show_column_names = TRUE,
          show_heatmap_legend = isTRUE(input$show_color_key),
          row_dend_obj = if (input$cluster_meth == 1 && dendrogram_selection()$row) {
            row_dendrogram()
          } else {
            NULL
          },
          col_dend_obj = if (input$cluster_meth == 1 && dendrogram_selection()$sample) {
            column_dendrogram()
          } else {
            NULL
          },
          use_letter_overlay = isTRUE(input$letter_overlay),
          show_cluster_labels = isTRUE(input$cluster_enrichment) && input$cluster_meth == 2,
          custom_cluster_labels = if (isTRUE(input$cluster_enrichment) && input$cluster_meth == 2) {
            tryCatch(cluster_pathway_labels(), error = function(e) NULL)
          } else {
            NULL
          }
        )
      )

      return(obj)
    })


    # Replace default download button with actionLink/tooltip for heatmaps
    setup_download_link <- function(
        ui_id,
        trigger_id,
        figure,
        filename,
        default_width,
        default_height,
        label_tag = NULL,
        icon_tag = NULL,
        tooltip_text = "Click to download plot in preferred format and size.") {
      min_size <- 2
      max_size <- 30
      width_id <- paste0(trigger_id, "_width")
      height_id <- paste0(trigger_id, "_height")
      pdf_id <- paste0(trigger_id, "_pdf")
      png_id <- paste0(trigger_id, "_png")
      svg_id <- paste0(trigger_id, "_svg")

      draw_download_plot <- function(plot_obj) {
        if (inherits(plot_obj, "idep_heatmap_bundle")) {
          ComplexHeatmap::draw(
            plot_obj$heatmap,
            annotation_legend_list = plot_obj$legends,
            annotation_legend_side = "top"
          )
        } else {
          print(plot_obj)
        }
      }

      output[[ui_id]] <- renderUI({
        req(figure())
        tagList(
          actionLink(
            inputId = ns(trigger_id),
            label = label_tag,
            icon = icon_tag,
            title = tooltip_text,
            `aria-label` = tooltip_text
          ),
          tippy::tippy_this(
            ns(trigger_id),
            tooltip_text,
            theme = "light"
          )
        )
      })

      width_value <- reactive({
        value <- input[[width_id]]
        if (is.numeric(value)) {
          return(max(min_size, min(max_size, value, na.rm = TRUE)))
        }
        default_width
      })

      height_value <- reactive({
        value <- input[[height_id]]
        if (is.numeric(value)) {
          return(max(min_size, min(max_size, value, na.rm = TRUE)))
        }
        default_height
      })

      observeEvent(input[[trigger_id]], {
        req(figure())
        showModal(modalDialog(
          numericInput(
            inputId = ns(width_id),
            label = "Width (in)",
            value = default_width,
            min = min_size,
            max = max_size
          ),
          numericInput(
            inputId = ns(height_id),
            label = "Height (in)",
            value = default_height,
            min = min_size,
            max = max_size
          ),
          h5("The plot will be rendered differently depending on size.\n            When the dimensions are too small, error or blank plot\n               will be generated."),
          downloadButton(outputId = ns(pdf_id), label = "PDF"),
          downloadButton(outputId = ns(png_id), label = "PNG"),
          downloadButton(outputId = ns(svg_id), label = "SVG"),
          size = "s"
        ))
      })

      output[[pdf_id]] <- downloadHandler(
        filename = paste0(filename, ".pdf"),
        content = function(file) {
          plot_obj <- figure()
          req(plot_obj)
          on.exit(removeModal(), add = TRUE)
          pdf(file, width = width_value(), height = height_value())
          draw_download_plot(plot_obj)
          dev.off()
        }
      )

      output[[png_id]] <- downloadHandler(
        filename = paste0(filename, ".png"),
        content = function(file) {
          plot_obj <- figure()
          req(plot_obj)
          on.exit(removeModal(), add = TRUE)
          png(
            filename = file,
            res = 360,
            width = width_value(),
            height = height_value(),
            units = "in"
          )
          draw_download_plot(plot_obj)
          dev.off()
        }
      )

      output[[svg_id]] <- downloadHandler(
        filename = paste0(filename, ".svg"),
        content = function(file) {
          plot_obj <- figure()
          req(plot_obj)
          on.exit(removeModal(), add = TRUE)
          svg(file, width = width_value(), height = height_value())
          draw_download_plot(plot_obj)
          dev.off()
        }
      )
    }

    setup_download_link(
      ui_id = "dl_heatmap_main_download_ui",
      trigger_id = "dl_heatmap_main_download",
      figure = heatmap_main_object,
      filename = "heatmap_main",
      default_width = 6,
      default_height = 16,
      label_tag = NULL,
      icon_tag = icon("download")
    )

    # Heatmap Click Value - DISABLED
    # Click functionality has been disabled to only show sub heatmap on brush
    # output$ht_click_content <- renderUI({
    #   req(!is.null(current_method()))
    #
    #   main_ready <- !is.null(shiny_env$ht) &&
    #     !is.null(shiny_env$ht_pos_main) &&
    #     length(shiny_env$ht@ht_list) >= 1
    #
    #   if (!is.null(input$ht_main_click) && main_ready) {
    #     group_info <- tryCatch(main_heatmap_group_info(), error = function(e) NULL)
    #     if (!is.null(group_info)) {
    #       matrix_obj <- shiny_env$ht@ht_list[[1]]
    #       if (!is.null(matrix_obj@matrix)) {
    #         click_data_main <- matrix_obj@matrix
    #         groups_main <- group_info$groups
    #         if (is.factor(groups_main)) {
    #           groups_main <- as.character(groups_main)
    #         }
    #         return(cluster_heat_click_info(
    #           click = input$ht_main_click,
    #           ht_sub = shiny_env$ht,
    #           ht_sub_obj = matrix_obj,
    #           ht_pos_sub = shiny_env$ht_pos_main,
    #           sub_groups = groups_main,
    #           group_colors = group_info$group_colors,
    #           cluster_meth = current_method(),
    #           click_data = click_data_main
    #         ))
    #       }
    #     }
    #   }
    #
    #   return(NULL)
    # })

    # depending on the number of genes selected
    # change the height of the sub heatmap
    height_sub_heatmap <- reactive({
      if (is.null(input$ht_brush)) {
        return(400)
      }

      # Ensure we have valid heatmap objects before trying to access positions
      if (is.null(shiny_env$ht) || is.null(shiny_env$ht_pos_main)) {
        return(400)
      }

      # Get the row ids of selected genes
      tryCatch({
        lt <- InteractiveComplexHeatmap::getPositionFromBrush(input$ht_brush)
        pos1 <- lt[[1]]
        pos2 <- lt[[2]]
        pos <- InteractiveComplexHeatmap::selectArea(
          shiny_env$ht,
          mark = FALSE,
          pos1 = pos1,
          pos2 = pos2,
          verbose = FALSE,
          ht_pos = shiny_env$ht_pos_main
        )

        row_index_list <- tryCatch(as.list(pos$row_index),
          error = function(e) NULL
        )
        total_rows <- 0
        if (!is.null(row_index_list)) {
          total_rows <- sum(vapply(row_index_list, length, integer(1)))
        }
        if (total_rows == 0) {
          total_rows <- length(unlist(pos[1, "row_index"]))
        }

        # convert to height, pixels
        height1 <- max(
          400, # minimum
          min(
            2000000, # maximum
            12 * total_rows
          )
        )
        return(height1)
      }, error = function(e) {
        # If there's an error (e.g., viewport not found), return default height
        return(400)
      })
    })

    # Subheatmap creation ---------
    output$sub_heatmap <- renderPlot(
      {
        if (is.null(input$ht_brush)) {
          grid::grid.newpage()
          # Only show message if main heatmap data is available
          if (!is.null(heatmap_data()) && !is.null(input$select_factors_heatmap)) {
            grid::grid.text(
              "Drag over a region on the heatmap to zoom in.",
              x = 0.5,
              y = 0.5,
              gp = grid::gpar(
                fontsize = 14,
                col = "#666666",
                fontface = "italic"
              )
            )
          }
          return(invisible(NULL))
        }

        withProgress(message = "Creating sub-heatmap", value = 0, {

          submap_return <- heatmap_sub_object_calc()
          if (is.null(submap_return)) {
            grid::grid.newpage()
            grid::grid.text(
              "Unable to create sub-heatmap. Please try selecting a different region.",
              x = 0.5,
              y = 0.5,
              gp = grid::gpar(
                fontsize = 14,
                col = "#666666",
                fontface = "italic"
              )
            )
            return(invisible(NULL))
          }


          shiny_env$submap_data <- submap_return$submap_data

          ht_static <- ComplexHeatmap::draw(
            submap_return$ht_select,
            annotation_legend_list = submap_return$lgd,
            annotation_legend_side = "top"
          )


          ht_static
        })
      },
      # adjust height of the zoomed in heatmap dynamically based on selection
      height = reactive(height_sub_heatmap())
#      ,width = 500 # this avoids the heatmap being redraw
    )

    # Reactive input versions to store values every submit press
    selected_factors_heatmap <- reactive({
      req(!is.na(input$select_factors_heatmap))
      input$select_factors_heatmap
    })
    
    submitted_pal <- reactive({
      input$sample_color
    })

    current_method <- reactive({
      input$cluster_meth
    })
    
    heatmap_sub_object_calc <- reactive({
      req(!is.null(submitted_pal()))
      req(!is.null(selected_factors_heatmap()))
      req(!is.null(input$ht_brush))

      # Ensure we have valid heatmap objects before processing brush
      # These are set when the main heatmap renders
      if (is.null(shiny_env$ht) || is.null(shiny_env$ht_pos_main)) {
        return(NULL)
      }

      submap_return <- tryCatch({ # tolerates error; otherwise stuck with spinner
        heat_sub(
          ht_brush = input$ht_brush,
          ht = shiny_env$ht,
          ht_pos_main = shiny_env$ht_pos_main,
          heatmap_data = heatmap_data(),
          sample_info = pre_process$sample_info(),
          select_factors_heatmap = selected_factors_heatmap(),
          cluster_meth = current_method(),
          group_pal = if (selected_factors_heatmap() == "All factors") {
            group_pal()
          } else {
            NULL
          },
          sample_color = submitted_pal(),
          use_letter_overlay = isTRUE(input$letter_overlay)
        )},
        error = function(e) {e$message}
      )

      if ("character" %in% class(submap_return)){
        submap_return <- NULL
      }

      if (!is.null(dim(submap_return$ht_select))){
        if (nrow(submap_return$ht_select) == 0 ||
            ncol(submap_return$ht_select) == 0) {
          submap_return <- NULL
        }
      }

      if (!is.null(submap_return)) {
        shiny_env$submap_data <- submap_return$submap_data
      } else {
        shiny_env$submap_data <- NULL
      }

      return(submap_return)
    })
    # Subheatmap creation ---------
    heatmap_sub_object <- reactive({
      if (is.null(input$ht_brush)) {
        return(NULL)
      }

      submap_return <- heatmap_sub_object_calc()
      if (is.null(submap_return)) {
        return(NULL)
      }

      shiny_env$submap_data <- submap_return$submap_data
      structure(
        list(
          heatmap = submap_return$ht_select,
          legends = submap_return$lgd
        ),
        class = "idep_heatmap_bundle"
      )
    })

    setup_download_link(
      ui_id = "dl_heatmap_sub_download_ui",
      trigger_id = "dl_heatmap_sub_download",
      figure = heatmap_sub_object,
      filename = "heatmap_zoom",
      default_width = 8,
      default_height = 12,
      label_tag = icon("download"),
      icon_tag = NULL
    )

    # gene lists for enrichment analysis - Hierarchical clustering
    # This reactive depends on the brush selection for hierarchical clustering
    hierarchical_gene_lists <- reactive({
      req(!is.null(pre_process$select_gene_id()))
      req(!is.null(input$ht_brush))
      req(input$cluster_meth == 1)
      req(!is.null(shiny_env$submap_data))
      req(is.data.frame(shiny_env$submap_data) || is.matrix(shiny_env$submap_data))
      req(nrow(shiny_env$submap_data) > 0)

      gene_lists <- list()

      gene_names <- merge_data(
        all_gene_names = pre_process$all_gene_names(),
        data = shiny_env$submap_data,
        merge_ID = pre_process$select_gene_id()
      )

      # Only keep the gene names and scrap the data
      gene_lists[["Selection"]] <- dplyr::select_if(gene_names, is.character)

      return(gene_lists)
    })

    # gene lists for enrichment analysis - k-means clustering
    # This reactive does NOT depend on brush selection, only on cluster assignments
    kmeans_gene_lists <- reactive({
      req(!is.null(pre_process$select_gene_id()))
      req(input$cluster_meth == 2)
      req(heatmap_data())
      req(input$k_clusters)
      req(pre_process$select_gene_id())
      req(shiny_env$ht)

      gene_lists <- list()

      row_ord <- ComplexHeatmap::row_order(shiny_env$ht)

      req(!is.null(names(row_ord)))

      for (i in 1:length(row_ord)) {
        if (i == 1) {
          clusts <- data.frame(
            "cluster" = rep(names(row_ord[i]), length(row_ord[[i]])),
            "row_order" = row_ord[[i]]
          )
        } else {
          tem <- data.frame(
            "cluster" = rep(names(row_ord[i]), length(row_ord[[i]])),
            "row_order" = row_ord[[i]]
          )
          clusts <- rbind(clusts, tem)
        }
      }
      clusts$id <- rownames(heatmap_data()[clusts$row_order, ])

      req(length(unique(clusts$cluster)) == input$k_clusters)
      # disregard user selection use clusters for enrichment
      for (i in 1:input$k_clusters) {
        cluster_data <- subset(clusts, cluster == i)
        row.names(cluster_data) <- cluster_data$id

        gene_names <- merge_data(
          all_gene_names = pre_process$all_gene_names(),
          data = cluster_data,
          merge_ID = pre_process$select_gene_id()
        )

        # Only keep the gene names and scrap the data
        gene_lists[[paste0("", i)]] <-
          dplyr::select_if(gene_names, is.character)
      }

      return(gene_lists)
    })

    k_means_list <- reactive({
      req(!is.null(kmeans_gene_lists()))
      kmeans_gene_lists()
    })

    gene_list_clust <- reactive({
      req(!is.null(input$cluster_meth))

      if (current_method() == 1){
        # For hierarchical clustering, only return gene lists if we have valid data
        req(!is.null(input$ht_brush))
        req(!is.null(shiny_env$submap_data))
        hierarchical_gene_lists()
      } else {
        # For k-means clustering, gene lists are based on clusters
        req(!is.null(k_means_list()))
        k_means_list()
      }
    })

    output$cloud_ui <- renderUI({
      req(!is.null(k_means_list()))
      tagList(
        selectInput(
          label = "Select Cluster:",
          inputId = ns("select_cluster"),
          choices = unique(names(k_means_list())),
          selected = unique(names(k_means_list()))[1],
          selectize = FALSE
        ),
        selectInput(
          label = "Select GO:",
          inputId = ns("cloud_go"),
          choices = setNames(
            c( "KEGG", "GOBP", "GOCC", "GOMF"),
            c("KEGG",
              "GO Biological Process",
              "GO Cellular Component",
              "GO Molecular Function")
          ),
          selectize = FALSE
        )
      )
    })
    
    # Sample Tree ----------
    sample_tree <- reactive({
      req(!is.null(pre_process$data()), input$cluster_meth == 1)

      draw_sample_tree(
        tree_data = pre_process$data(),
        gene_centering = input$gene_centering,
        gene_normalize = input$gene_normalize,
        sample_centering = FALSE,
        sample_normalize = FALSE,
        hclust_funs = hclust_funs,
        hclust_function = input$hclust_function,
        dist_funs = dist_funs,
        dist_function = input$dist_function
      )
      p <- recordPlot()
      return(p)
    })

    output$sample_tree <- renderPlot({
      print(sample_tree())
    })

    dl_sample_tree <- ottoPlots::mod_download_figure_server(
      id = "dl_sample_tree",
      filename = "sample_tree",
      figure = reactive({
        sample_tree()
      }),
      label = ""
    )

    # Handle clustering method changes ----------
    observeEvent(input$cluster_meth, {
      if (input$cluster_meth == 1) {
        # Hierarchical clustering
        showTab(
          inputId = "cluster_panels",
          target = "sample_tab"
        )
        hideTab(
          inputId = "cluster_panels",
          target = "word_cloud"
        )
      } else if (input$cluster_meth == 2) {
        # K-means clustering
        hideTab(
          inputId = "cluster_panels",
          target = "sample_tab"
        )
        showTab(
          inputId = "cluster_panels",
          target = "word_cloud"
        )
      }

      # Reset enrichment checkbox when switching methods
      updateCheckboxInput(
        session = session,
        inputId = "cluster_enrichment",
        value = FALSE
      )

      # Clear the brush selection by sending JavaScript to reset it
      session$sendCustomMessage(
        type = "resetBrush",
        message = list(brushId = ns("ht_brush"))
      )

      # Clear sub-heatmap environment variables to prevent viewport errors
      # when switching between clustering methods
      # NOTE: We do NOT clear shiny_env$ht or shiny_env$ht_pos_main here
      # because the new heatmap will overwrite them when it renders
      shiny_env$submap_data <- NULL
    })

    # Auto-uncheck enrichment checkbox when it should be hidden ----------
    observe({
      # If hierarchical clustering is selected and no region is brushed,
      # uncheck the enrichment checkbox to prevent it from staying checked
      # while hidden
      if (input$cluster_meth == 1 && is.null(input$ht_brush)) {
        if (!is.null(input$cluster_enrichment) && input$cluster_enrichment) {
          updateCheckboxInput(
            session = session,
            inputId = "cluster_enrichment",
            value = FALSE
          )
        }
      }
    })

    # k-Cluster elbow plot ----------
    output$k_clusters <- renderPlot({
      req(!is.null(heatmap_data()))

      k_means_elbow(
        heatmap_data = heatmap_data()
      )
    })
    # pop-up modal
    observeEvent(input$elbow_pop_up, {
      showNotification(
        ui = "Generating plot. May take 5-10 seconds...",
        id = "elbow_pop_up_message",
        duration = 6,
        type = "message"
      )
      showModal(modalDialog(
        plotOutput(ns("k_clusters")),
        footer = NULL,
        easyClose = TRUE,
        title = tags$h5(
          "Following the elbow method, one should choose k so that adding
          another cluster does not substantially reduce the within groups sum of squares.",
          tags$a(
            "Wikipedia",
            href = "https://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set",
            target = "_blank"
          )
        ),
      ))
    })

    # Heatmap Download Data -----------
    heatmap_data_download <- reactive({
      req(!is.null(pre_process$all_gene_names()))
      req(!is.null(heatmap_data()))

      data <- prep_download(
        heatmap = heatmap_main_object(),
        heatmap_data = heatmap_data(),
        cluster_meth = input$cluster_meth
      )

      merged_data <- merge_data(
        pre_process$all_gene_names(),
        data,
        merge_ID = pre_process$select_gene_id()
      )
    })

    output$download_heatmap_data <- downloadHandler(
      filename = function() {
        "heatmap_data.csv"
      },
      content = function(file) {
        write.csv(heatmap_data_download(), file)
      }
    )
    
    enrichment_table_cluster <- mod_11_enrichment_server(
      id = "enrichment_table_cluster",
      gmt_choices = reactive({
        pre_process$gmt_choices()
      }),
      gene_lists = reactive({
        gene_list_clust()
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

    # Extract top pathway names for each cluster for heatmap labels
    cluster_pathway_labels <- reactive({
      # Only compute labels when enrichment is enabled and using k-means
      req(input$cluster_enrichment == TRUE)
      req(input$cluster_meth == 2)
      req(!is.null(enrichment_table_cluster$pathway_table()))

      pathway_table <- enrichment_table_cluster$pathway_table()

      # Extract top pathway for each cluster
      labels <- sapply(1:input$k_clusters, function(i) {
        cluster_name <- as.character(i)
        if (cluster_name %in% names(pathway_table)) {
          cluster_data <- pathway_table[[cluster_name]]
          # Check if data frame has rows and Pathway column
          if (is.data.frame(cluster_data) && nrow(cluster_data) > 0 && "Pathway" %in% colnames(cluster_data)) {
            # Get the first (most significant) pathway
            top_pathway <- cluster_data$Pathway[1]
            # Truncate long pathway names
            if (nchar(top_pathway) > 30) {
              top_pathway <- paste0(substr(top_pathway, 1, 27), "...")
            }
            return(top_pathway)
          }
        }
        # Fallback to default label if no pathway found
        return(paste("Cluster", i))
      })

      return(labels)
    })
    
    # Generate word/frequency data for word cloud
    word_cloud_data <- reactive({
      req(!is.na(input$select_cluster))
      req(!is.null(input$cloud_go))
      req(!is.null(k_means_list()))

      withProgress(message = "Creating Word Cloud", value = 0.5, {
        prep_cloud_data(gene_lists = k_means_list(),
                        cluster = input$select_cluster,
                        cloud_go = input$cloud_go,
                        select_org = pre_process$select_org(),
                        converted = pre_process$converted(),
                        gmt_file = pre_process$gmt_file(),
                        idep_data = idep_data,
                        gene_info = pre_process$all_gene_info())
      })
    })

    output$word_cloud <- wordcloud2::renderWordcloud2({
      req(!is.null(word_cloud_data()))
      
      if ("character" %in% class(word_cloud_data())){
        NULL
      } else {
        
        wordcloud2::wordcloud2(word_cloud_data(),
                               shape = "circle",
                               rotateRatio = 0,
                               color = "random-dark",
                               shuffle = FALSE)
      }
    })
    
    # Error message UI for word cloud
    output$cloud_error <- renderUI({
      req(!is.null(word_cloud_data()))
      
      if ("character" %in% class(word_cloud_data())){
        div(style = "color:red;",
            "Pathways Not Found for selected cluster!")
      } else {NULL}
    })
    
    output$cloud_download <- downloadHandler(
      filename = "word_cloud_data.csv",
      content = function(file) {
        req(!is.null(word_cloud_data()))
        
        write.csv(word_cloud_data(), file)
      }
    )
    
    # Markdown report------------
    output$report <- downloadHandler(

      # For PDF output, change this to "report.pdf"
      filename = paste0(
        "clustering_report_",
        format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
        ".html"
      ),
      content = function(file) {
        withProgress(message = "Generating report", {
          incProgress(0.2)
          # Copy the report file to a temporary directory before processing it, in
          # case we don't have write permissions to the current working dir (which
          # can happen when deployed).
          tempReport <- file.path(tempdir(), "clustering_workflow.Rmd")
          # tempReport
          tempReport <- gsub("\\", "/", tempReport, fixed = TRUE)

          # This should retrieve the project location on your device:
          # "C:/Users/bdere/Documents/GitHub/idepGolem"
          wd <- getwd()

          markdown_location <- app_sys("app/www/RMD/clustering_workflow.Rmd")
          file.copy(from = markdown_location, to = tempReport, overwrite = TRUE)

          # Set up parameters to pass to Rmd document
          params <- list(
            pre_processed_data = pre_process$data(),
            sample_info = pre_process$sample_info(),
            descr = pre_process$descr(),
            mapping_statistics = pre_process$mapping_statistics(),
            all_gene_names = pre_process$all_gene_names(),
            n_genes = input$n_genes,
            k_clusters = input$k_clusters,
            cluster_meth = input$cluster_meth,
            select_gene_id = pre_process$select_gene_id(),
            list_factors_heatmap = input$list_factors_heatmap,
            heatmap_color_select = heatmap_color_select(),
            dist_function = input$dist_function,
            hclust_function = input$hclust_function,
            heatmap_cutoff = input$heatmap_cutoff,
            gene_centering = input$gene_centering,
            gene_normalize = input$gene_normalize,
            sample_clustering = dendrogram_selection()$sample,
            show_row_dend = dendrogram_selection()$row,
            selected_genes = input$selected_genes,
            submap_data = if(!is.null(shiny_env$submap_data)) shiny_env$submap_data else NULL,
            select_factors_heatmap = input$select_factors_heatmap,
            sample_color = input$sample_color
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
  })
}

## To be copied in the UI
# mod_03_heatmap_ui("03_heatmap_ui_1")

## To be copied in the server
# mod_03_heatmap_server("03_heatmap_ui_1")
