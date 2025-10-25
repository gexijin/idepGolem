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
    prep_data <- pre_process$data()
    placeholder_present <- !is.null(prep_data) && "Ctrl" %in% colnames(prep_data) && ncol(prep_data) == 2
    placeholder_is_zero <- placeholder_present && all(prep_data[, "Ctrl"] == 0)

    # Hide tabs for type 3 data (fold-change & P-values) regardless of number of comparisons
    # This includes both single comparison (with dummy Ctrl) and multiple comparisons
    hide_prep <- isTRUE(data_format == 3) && placeholder_is_zero
    hide_analysis_tabs <- isTRUE(data_format == 3)  # Hide PCA, Bicluster, Network for all type 3 data

    top_targets <- c("Prep", "Cluster", "PCA", "Bicluster", "Network")
    analysis_only_targets <- c("PCA", "Bicluster", "Network")  # Tabs to hide for type 3 with multiple comparisons
    venn_input_id <- "deg-step_1"
    venn_target <- "venn_diagram"
    deg_plots_input_id <- "deg-step_2"
    deg_plots_targets <- c("Heatmap", "MA Plot", "Scatter Plot", "R Code")
    deg_plots_fallback <- "Volcano Plot"

    if (hide_prep) {
      # Single comparison: hide all tabs (original behavior)
      lapply(top_targets, function(tab_name) hideTab(inputId = "navbar", target = tab_name))
      hideTab(inputId = venn_input_id, target = venn_target)
      lapply(deg_plots_targets, function(tab_name) hideTab(inputId = deg_plots_input_id, target = tab_name))

      if (isolate(input$navbar) %in% top_targets) {
        updateNavbarPage(session, inputId = "navbar", selected = "Data")
      }
      if (identical(isolate(input[[venn_input_id]]), venn_target)) {
        updateTabsetPanel(session, inputId = venn_input_id, selected = "results")
      }
      if (isolate(input[[deg_plots_input_id]]) %in% deg_plots_targets) {
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

      # Show venn diagram and deg plot tabs
      showTab(inputId = venn_input_id, target = venn_target)
      lapply(deg_plots_targets, function(tab_name) showTab(inputId = deg_plots_input_id, target = tab_name))

      if (isolate(input$navbar) %in% analysis_only_targets) {
        updateNavbarPage(session, inputId = "navbar", selected = "Data")
      }
    } else {
      # Not type 3 data: show all tabs
      lapply(top_targets, function(tab_name) showTab(inputId = "navbar", target = tab_name))
      showTab(inputId = venn_input_id, target = venn_target)
      lapply(deg_plots_targets, function(tab_name) showTab(inputId = deg_plots_input_id, target = tab_name))
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
