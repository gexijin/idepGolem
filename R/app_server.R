#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # file size is 5MB by default. This changes it to 30MB
  # options(shiny.maxRequestSize = 30*1024^2)
  options(warn = -1) # turn off warning
  pdf(file = NULL)

  # Database configuration and idep_data are now loaded globally in run_app()
  # This means all user sessions share the same species database (13,391 species)
  # saving ~6.4 MB of memory and ~0.3 seconds of startup time per user

  # Tab Variable to control reactivity
  tab <- reactive(input$navbar)

  load_data <- mod_01_load_data_server(
    id = "load_data",
    idep_data = idep_data,
    tab = tab
  )
  pre_process <- mod_02_pre_process_server(
    id = "pre_process",
    load_data = load_data,
    tab = tab
  )
  single_column_notice_shown <- reactiveVal(FALSE)
  single_column_notice_id <- "single_column_notice"
  observe({
    data_format <- load_data$data_file_format()
    data_format_numeric <- suppressWarnings(as.numeric(data_format))
    has_data_format <- length(data_format_numeric) > 0 && !is.na(data_format_numeric)
    summary_formats <- c(3, 4)
    is_summary_format <- has_data_format && data_format_numeric %in% summary_formats
    prep_data <- pre_process$data()
    placeholder_present <- !is.null(prep_data) && "Ctrl" %in% colnames(prep_data) && ncol(prep_data) == 2
    placeholder_is_zero <- placeholder_present && all(prep_data[, "Ctrl"] == 0)

    # Hide tabs for type 3 data (fold-change & P-values) regardless of number of comparisons
    # This includes both single comparison (with dummy Ctrl) and multiple comparisons
    hide_prep <- has_data_format && data_format_numeric == 3 && placeholder_is_zero
    hide_analysis_tabs <- isTRUE(is_summary_format)  # Hide PCA, Bicluster, Network for summary-level data
    is_type_3 <- has_data_format && data_format_numeric == 3  # Fold-change + FDR
    is_type_4 <- has_data_format && data_format_numeric == 4  # Fold-change only (no FDR)

    top_targets <- c("Prep", "Cluster", "PCA", "Bicluster", "Network")
    analysis_only_targets <- c("PCA", "Bicluster", "Network")  # Tabs to hide for type 3 with multiple comparisons
    venn_input_id <- "deg-step_1"
    venn_target <- "venn_diagram"
    deg_plots_input_id <- "deg-step_2"
    # Type 3: hide MA Plot and Scatter Plot
    deg_plots_type3_targets <- c("MA Plot", "Scatter Plot")
    # Type 4: hide all plot tabs that require expression data or FDR
    deg_plots_type4_targets <- c("Volcano Plot", "MA Plot", "Scatter Plot", "R Code")
    deg_plots_fallback <- "Genes"

    if (hide_prep) {
      # Single comparison: hide all tabs (original behavior)
      lapply(top_targets, function(tab_name) hideTab(inputId = "navbar", target = tab_name))
      hideTab(inputId = venn_input_id, target = venn_target)
      lapply(deg_plots_type3_targets, function(tab_name) hideTab(inputId = deg_plots_input_id, target = tab_name))

      if (isolate(input$navbar) %in% top_targets) {
        updateNavbarPage(session, inputId = "navbar", selected = "Data")
      }
      if (identical(isolate(input[[venn_input_id]]), venn_target)) {
        updateTabsetPanel(session, inputId = venn_input_id, selected = "results")
      }
      if (isolate(input[[deg_plots_input_id]]) %in% deg_plots_type3_targets) {
        updateTabsetPanel(session, inputId = deg_plots_input_id, selected = deg_plots_fallback)
      }
      if (!single_column_notice_shown()) {
        showNotification(
          ui = HTML(paste(
            "Use the Stats tab to choose a cutoff to define differentially expressed genes,",
            "then visit the DEG tab for enrichment analysis based on those genes.",
            "Switch to the Pathway tab for pathway analysis using the fold-changes of all genes."
          )),
          id = single_column_notice_id,
          duration = 30,
          type = "message"
        )
        single_column_notice_shown(TRUE)
      }
    } else if (hide_analysis_tabs) {
      # Multiple comparisons: hide only PCA, Bicluster, Network (keep Prep and Cluster visible)
      lapply(analysis_only_targets, function(tab_name) hideTab(inputId = "navbar", target = tab_name))

      # Show Prep and Cluster tabs
      showTab(inputId = "navbar", target = "Prep")
      showTab(inputId = "navbar", target = "Cluster")

      # Show venn diagram
      showTab(inputId = venn_input_id, target = venn_target)

      # Handle DEG plot tabs based on data type
      if (is_type_4) {
        # Type 4 (fold-change only): hide all plot tabs
        lapply(deg_plots_type4_targets, function(tab_name) hideTab(inputId = deg_plots_input_id, target = tab_name))
        # Show Heatmap (not in hide list)
        showTab(inputId = deg_plots_input_id, target = "Heatmap")
        if (isolate(input[[deg_plots_input_id]]) %in% deg_plots_type4_targets) {
          updateTabsetPanel(session, inputId = deg_plots_input_id, selected = deg_plots_fallback)
        }
      } else if (is_type_3) {
        # Type 3 (fold-change + FDR): hide MA Plot and Scatter Plot only
        lapply(deg_plots_type3_targets, function(tab_name) hideTab(inputId = deg_plots_input_id, target = tab_name))
        # Show Volcano Plot, Heatmap, and R Code
        showTab(inputId = deg_plots_input_id, target = "Volcano Plot")
        showTab(inputId = deg_plots_input_id, target = "Heatmap")
        showTab(inputId = deg_plots_input_id, target = "R Code")
        if (isolate(input[[deg_plots_input_id]]) %in% deg_plots_type3_targets) {
          updateTabsetPanel(session, inputId = deg_plots_input_id, selected = "Volcano Plot")
        }
      }

      if (isolate(input$navbar) %in% analysis_only_targets) {
        updateNavbarPage(session, inputId = "navbar", selected = "Data")
      }
    } else {
      # Not summary data (types 1 and 2): show all tabs
      lapply(top_targets, function(tab_name) showTab(inputId = "navbar", target = tab_name))
      showTab(inputId = venn_input_id, target = venn_target)
      # Show all DEG plot tabs for regular data types
      all_deg_tabs <- c("Heatmap", "Volcano Plot", "MA Plot", "Scatter Plot", "R Code")
      lapply(all_deg_tabs, function(tab_name) showTab(inputId = deg_plots_input_id, target = tab_name))
      if (single_column_notice_shown()) {
        removeNotification(id = single_column_notice_id)
        single_column_notice_shown(FALSE)
      }
    }
  })

  mod_03_clustering_server(
    id = "clustering",
    pre_process = pre_process,
    load_data = load_data,
    idep_data = idep_data,
    tab = tab
  )
  mod_04_pca_server(
    id = "pca",
    load_data = load_data,
    pre_process = pre_process,
    idep_data = idep_data
  )
  deg <- mod_05_deg_server(
    id = "deg",
    pre_process = pre_process,
    idep_data = idep_data,
    load_data = load_data,
    tab = tab
  )
  mod_06_pathway_server(
    id = "pathway",
    pre_process = pre_process,
    deg = deg,
    idep_data = idep_data,
    tab = tab
  )
  mod_07_genome_server(
    id = "genome",
    pre_process = pre_process,
    deg = deg,
    idep_data = idep_data
  )
  mod_08_bicluster_server(
    id = "bicluster",
    pre_process = pre_process,
    idep_data = idep_data,
    tab = tab
  )
  mod_09_network_server(
    id = "network",
    pre_process = pre_process,
    idep_data = idep_data,
    tab = tab
  )
  mod_10_doc_server(
    id = "doc",
    pre_process = pre_process,
    idep_data = idep_data,
    tab = tab
  )
}
