#' 01_load_data UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_01_load_data_ui <- function(id) {
  ns <- shiny::NS(id)
  tabPanel(
    title = "Data",
    # move notifications and progress bar to the center of the screen
    tags$head(
      tags$style(
        HTML(".shiny-notification {
              width: 300px;
              position:fixed;
              top: calc(85%);
              left: calc(5%);
               }

             #shiny-notification-load_prompt {
               width: 600px;
               max-width: 60vw;
               top: calc(10%);
               right: 20px;
               left: auto;
               background-color: #f8f9fa;
               opacity: 1;
             }

             .dis_gray { background-color: gray; }
             .more-options { display: block; width: 100%; }
             .more-options summary { display: flex; align-items: center; cursor: pointer; font-weight: 600; margin: 0; }
             .more-options summary::marker, .more-options summary::-webkit-details-marker { display: none; }
             .more-options summary::before { content: '+'; margin-right: 6px; font-size: 14px; line-height: 1; }
             .more-options[open] summary::before { content: '\u2212'; }
             .more-options-body { margin-top: 10px; padding: 10px 0; background-color: #f7f9fc; border-radius: 0; width: 100%; }
             .more-options-body .shiny-input-container > label { font-weight: 400; }
             .load-data-info-icon { padding: 0; border: none; background: transparent; color: #6c757d; font-size: 18px; line-height: 1; }
             .load-data-info-icon:hover { color: #343a40; text-decoration: none; }
             .load-data-info-icon:focus { outline: none; box-shadow: none; }
              ")
      )
    ),
    sidebarLayout(
      ##################################################################
      #       Load Data sidebar panel ----
      ##################################################################
      sidebarPanel(
        uiOutput(ns("reset_button")),
        # Species Match Drop Down ------------
        # Species Match Drop Down ------------
        fluidRow(
          column(
            width = 5,
            strong("1. Species"),
          ),
          column(
            width = 7,
            textOutput(ns("selected_species"))
          )
        ),
        tags$head(tags$style("#load_data-selected_species{color: blue;
                                 font-size: 12px;
                                 font-style: italic;
                                 }")),

        # .GMT file input bar ----------
        fluidRow(
          column(
            width = 6,
            align = "center",
            # Species list and genome assemblies ----------
            actionButton(
              inputId = ns("genome_assembl_button"),
              label = "Select"
            ),
            tippy::tippy_this(
              ns("genome_assembl_button"),
              "Select your study organism. Human is the default.",
              theme = "light"
            )
          ),
          column(
            width = 6,
            align = "center",
            checkboxInput(
              inputId = ns("new_species"),
              label = "New",
              value = FALSE
            ),
            tippy::tippy_this(
              ns("new_species"),
              "Analyze data for a species not in our list. You can still do most analyses. For pathway analysis, upload a custom GMT file (below).",
              theme = "light"
            )
          )
        ),

        # Conditional GMT file upload for new species ----------
        conditionalPanel(
          condition = "input.new_species == true",
          br(),
          p("Custom pathway file (optional)"),
          fileInput(
            inputId = ns("gmt_file"),
            label = NULL,
            accept = c(
              "text/csv",
              "text/comma-separated-values",
              "text/tab-separated-values",
              "text/plain",
              ".csv",
              ".tsv",
              ".gmt"
            ),
            placeholder = "Upload a .GMT file for pathway analysis"
          ),
          tippy::tippy_this(
            ns("gmt_file"),
            "Optional. Upload a pathway .GMT file tailored to your species; skip this if you do not have one.",
            theme = "light"
          ),
          ns = ns
        ),

        # Dropdown for data file format ----------
        selectInput(
          inputId = ns("data_file_format"),
          label = "2. Data type",
          choices = list(
            "..." = 0,
            "Read counts data" = 1,
            "Normalized expression data" = 2,
            "Fold changes & adjusted p-values" = 3,
            "Fold changes only" = 4
          ),
          selected = 0,
          selectize = FALSE
        ),
        tippy::tippy_this(
          ns("data_file_format"),
          paste(
            "Raw read counts enable DESeq2.",
            "Choose normalized expression for TPM/FPKM, microarray, or proteomics.",
            "Select fold change + adjusted p-values when both statistics are available,",
            "or choose Fold changes only if you just have log-fold changes."
          ),
          theme = "light"
        ),
        # Load expression data options ----------
        # Includes load demo action button, demo data dropdown, and expression
        # file upload box
        conditionalPanel(
          condition = "input.data_file_format != 0",
          uiOutput(ns("load_data_ui")),
          ns = ns
        ),

        # Experiment design file input and interactive builder button ----------
        conditionalPanel(
          condition = "input.data_file_format != 0",
          uiOutput(ns("design_file_ui")),
          uiOutput(ns("wizard_btn_ui")),
          ns = ns
        ),
        tags$details(
          class = "more-options",
          tags$summary(span(strong("Settings"), id = ns("global_settings_summary"))),
          tippy::tippy_this(
            ns("global_settings_summary"),
            "Show shared display and ID-conversion settings for the entire app.",
            theme = "light"
          ),
          div(
            class = "more-options-body",
            div(
              id = ns("heatmap_color_select_container"),
              style = "cursor: help;",
              selectInput(
                inputId = ns("heatmap_color_select"),
                label = "Heatmap color:",
                choices = c(
                  "Green-Black-Red",
                  "Red-Black-Green",
                  "Blue-White-Red",
                  "Green-Black-Magenta",
                  "Blue-Yellow-Red",
                  "Blue-White-Brown",
                  "Orange-White-Blue"
                ),
                selected = "Green-Black-Red",
                width = "100%",
                selectize = FALSE
              )
            ),
            tippy::tippy_this(
              ns("heatmap_color_select_container"),
              "Select the color palette for heatmaps, sample trees, and network plots.",
              theme = "light"
            ),
            selectInput(
              inputId = ns("multiple_map"),
              label = "Multiple mapped IDs:",
              choices = list(
                "Sum" = "sum",
                "Average" = "mean",
                "Median" = "median",
                "Max" = "max",
                "Max SD" = "max_sd"
              ),
              selected = "Sum",
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("multiple_map"),
              "When several IDs match one gene, choose how to combine them (sum, average, median, max, or most variable). For transcript counts, select sum to get gene-level totals.",
              theme = "light"
            ),
            div(
              id = ns("plots_color_select_container"),
              style = "cursor: help;",
              selectInput(
                inputId = ns("plots_color_select"),
                label = "Plots Color:",
                choices = c(
                  "Set1",
                  "Set2",
                  "Set3",
                  "Paired",
                  "Dark2",
                  "Accent",
                  "Pastel1",
                  "Pastel2",
                  "Spectral"
                ),
                selected = "Set1",
                width = "100%",
                selectize = FALSE
              )
            ),
            tippy::tippy_this(
              ns("plots_color_select_container"),
              "Pick the color set for PCA and QC plots so sample groups are easy to see.",
              theme = "light"
            ),
            selectInput(
              inputId = ns("select_gene_id"),
              label = "Gene ID type for plots:",
              choices = c("symbol", "ensembl_ID", "User_ID"),
              selected = "symbol",
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("select_gene_id"),
              "Choose which gene identifier appears in plots and tables.",
              theme = "light"
            ),
            selectInput(
              inputId = ns("ggplot2_theme"),
              label = "ggplot2 theme:",
              choices = c(
                "default", # no change
                "gray",
                "bw",
                "light",
                "dark",
                "classic",
                "minimal",
                "linedraw"
              ),
              selected = "default",
              width = "100%",
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("ggplot2_theme"),
              "Apply a different overall style to every ggplot2 plot.",
              theme = "light"
            ),
            checkboxInput(
              inputId = ns("plot_grid_lines"),
              label = "Add grid lines to plots",
              value = FALSE
            ),
            tippy::tippy_this(
              ns("plot_grid_lines"),
              "Add light grid lines so sample values are easier to compare.",
              theme = "light"
            ),
            checkboxInput(
              inputId = ns("no_id_conversion"),
              label = "Do not convert gene IDs",
              value = FALSE
            ),
            tippy::tippy_this(
              ns("no_id_conversion"),
              "Keep your original gene IDs. By default we convert to Ensembl IDs used by pathway databases.",
              theme = "light"
            ),
            numericInput(
              inputId = ns("max_groups"),
              label = "Max groups in legends:",
              value = 12,
              min = 5,
              max = 30,
              step = 1
            ),
            tippy::tippy_this(
              ns("max_groups"),
              "Maximum number of groups to show in plot legends. Groups beyond this limit are recoded as 'Other'.",
              theme = "light"
            ),
            numericInput(
              inputId = ns("max_group_name_length"),
              label = "Max group name length:",
              value = 30,
              min = 10,
              max = 50,
              step = 5
            ),
            tippy::tippy_this(
              ns("max_group_name_length"),
              "Maximum character length for group names in legends. Longer names are truncated.",
              theme = "light"
            ),
          )
        ),
        br(),
        fluidRow(
          column(
            width = 4,
            align = "left",
            # Link to public RNA-seq datasets ----------
            a(
              "Public Data",
              href = "http://bioinformatics.sdstate.edu/reads/",
              style = "color: #1357A6;"
            ),
          ),
          column(
            width = 4,
            align = "center",
            a(
              "Cite iDEP",
              href = "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6#citeas",
              target = "_blank",
              style = "color: #1357A6;"
            ),
          ),
          column(
            width = 4,
            align = "right",
            # Action link to reveal Gene ID examples -----------
            actionLink(
              inputId = ns("gene_ids_link"),
              label = "Gene IDs",
              style = "color: #1357A6;"
            ),
            tippy::tippy_this(
              ns("gene_ids_link"),
              "See examples of the accepted gene identifiers for each species.",
              theme = "light"
            )
          )
        )
      ),

      ##################################################################
      #       Load Data panel main ----
      ##################################################################
      mainPanel(
        shinyjs::useShinyjs(),
        # Default species notification ----------
        uiOutput(ns("default_species_message")),

        # Gene ID conversion statistics ----------
        uiOutput(ns("conversion_stats_message")),
        div(
          style = "overflow-x: auto;",
          tableOutput(ns("sample_info_table"))
        ),

        # Display first 20 rows of the data ----------
        div(
          style = "overflow-x: auto;",
          tableOutput(ns("sample_20"))
        ),
        div(
          id = ns("load_message"),
          h1("From data to discoveries", style = "color: #B30000; font-weight: 700;"),
          br(),
          h3("Visualize, analyze, & unveil pathways — in minutes!", style = "color: #10652E;"),
          br(),
          br(),
          br(),
          h4("Loading R packages, please wait ... ... ..."),
        ),
        uiOutput(ns("show_gmt")),
        # Display file format help html document when prompted ----
        uiOutput(ns("format_help_ui")),
        # Hide welcome screen after data is loaded -----
        uiOutput(ns("welcome_ui"))
      )
    )
  )
}


#' 01_load_data Server Functions
#'
#' @noRd
### testing something, come back to later
mod_01_load_data_server <- function(id, idep_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    go_button_count <- reactive({
      if (is.null(input$go_button)) {
        0
      } else {
        input$go_button
      }
    })

    # Reactive value to store selected species
    select_org <- reactiveVal(NULL)

    selected_demo <- reactiveVal(NULL)
    demo_preview_content <- reactiveVal(NULL)

    # Stores a sample_info matrix built via the interactive design builder GUI
    gui_design <- reactiveVal(NULL)

    n_guided_factors <- reactiveVal(2L)
    guided_factor_defs <- reactiveVal(NULL)
    # Named list: factor index (character) → character vector of level per sample
    sample_assignments <- reactiveVal(list())

    demo_preview_tables <- list(
      `1` = data.frame(
        Ensembl = c(
          "ENSG00000198888",
          "ENSG00000198763",
          "ENSG00000198804",
          "ENSG00000198712",
          "ENSG00000228253"
        ),
        ctrl_1 = c(2, 90, 234, 3222, 9304),
        ctrl_2 = c(11, 66, 765, 2999, 11160),
        treat_1 = c(245, 54, 32, 3450, 12608),
        treat_2 = c(308, 73, 77, 287, 13041),
        stringsAsFactors = FALSE
      ),
      `2` = data.frame(
        symbol = c("A2M", "A4GALT", "AAAS", "AACS", "AADAC"),
        wt_1 = c(7.7757, 1.5048, 8.9854, 2.6064, 2.9547),
        wt_2 = c(7.8172, 1.6626, 9.0822, 3.7765, 2.3480),
        mutant_1 = c(8.5988, 3.9233, 8.8083, 3.9611, 3.8708),
        mutant_2 = c(8.6342, 2.8120, 8.6981, 2.9102, 4.3300),
        stringsAsFactors = FALSE
      ),
      `3` = data.frame(
        genes = c("Gnai3", "Cdc45", "Scml2", "Narf", "Cav2"),
        treatA_LFC = c(-0.0981, 0.510, 0.545, -0.229, -0.592),
        treatA_Pval = c(0.99, 0.001, 0.021, 0.120, 1e-10),
        treatB_LFC = c(-1.81, 0.10, 1.98, -0.02, 0.13),
        treatB_Pval = c(0.001, 0.25, 0.0001, 0.02, 0.001),
        stringsAsFactors = FALSE
      ),
      `4` = data.frame(
        genes = c("Gnai3", "Cdc45", "Scml2", "Narf", "Cav2"),
        WTvsMu_LFC = c(-0.098, 0.510, 0.545, -0.230, -0.592),
        IRvsMock_LFC = c(0.760, 1.638, 0.929, -0.973, -0.176),
        Interaction_LFC = c(-0.659, -1.499, 1.149, 0.162, 0.970),
        stringsAsFactors = FALSE
      )
    )

    build_preview_table <- function(df) {
      df[] <- lapply(df, as.character)

      tags$table(
        style = "width: 100%; border-collapse: collapse; margin-top: 10px;",
        tags$thead(
          tags$tr(
            lapply(names(df), function(col) {
              tags$th(
                style = "border: 1px solid #ddd; padding: 8px; text-align: center; background-color: #f2f2f2;",
                col
              )
            })
          )
        ),
        tags$tbody(
          lapply(seq_len(nrow(df)), function(i) {
            tags$tr(
              lapply(df[i, ], function(val) {
                tags$td(
                  style = "border: 1px solid #ddd; padding: 8px; text-align: center;",
                  val
                )
              })
            )
          })
        )
      )
    }

    get_data_type_details <- function(type) {
      switch(as.character(type),
        `1` = list(
          title = "Read Counts",
          body = tagList(
            p("Raw gene-by-sample count matrix. Values indicate the number of sequencing reads assigned to each gene. Counts are typically integers, but estimated counts (e.g., from kallisto or Salmon) may be non-integer and should still be treated as count data."),
            tags$ul(
              tags$li("First column contains gene IDs such as Ensembl (recommended), Entrez, symbols, ..."),
              tags$li("Column headers: sample names; append '_1', '_2', etc, for replicates; avoid spaces and '-' characters. ")
            )
          )
        ),
        `2` = list(
          title = "Normalized Expression Matrix",
          body = tagList(
            p("Provide a gene-by-sample matrix with normalized values such as TPM/FPKM, microarray intensities, proteomics, etc."),
            tags$ul(
              tags$li("First column contains gene IDs such as Ensembl (recommended), Entrez, symbols, ..."),
              tags$li("Column headers: sample names; append '_1', '_2', etc, for replicates; avoid spaces and '-' characters.")
            )
          )
        ),
        `3` = list(
          title = "Fold Change & Adjusted P-values",
          body = tagList(
            p("Upload summary statistics for one or more contrasts."),
            tags$ul(
              tags$li("First column contains gene IDs such as Ensembl (recommended), Entrez, symbols, ..."),
              tags$li("Log2 fold-change columns paired with adjusted P-value/FDR columns."),
              tags$li("Order columns as LFC, FDR, LFC, FDR for each contrast.")
            )
          )
        ),
        `4` = list(
          title = "Fold Changes Only",
          body = tagList(
            p("Upload log2 fold-changes when p-values or FDRs are unavailable."),
            tags$ul(
              tags$li("First column contains gene IDs such as Ensembl (recommended), Entrez, symbols, ..."),
              tags$li("One column per contrast containing log2 fold-change values.")
            )
          )
        ),
        list(
          title = "Upload Data",
          body = p("Upload the appropriate file for the selected data type.")
        )
      )
    }

    # increase max input file size
    options(shiny.maxRequestSize = 200 * 1024^2)

    # Initialize species selection with first species as default
    observe({
      if (!is.null(idep_data$org_info) && nrow(idep_data$org_info) > 0) {
        # Get first species from org_info (ordered by 'top' field)
        first_species_id <- idep_data$org_info$id[1]
        first_species_name <- idep_data$org_info$name2[1]

        # Only set the selected value with minimal choices (just the selected one)
        # The selectInput is hidden and only used to store the selection state
        # Users select via the modal DataTable, not the dropdown
        # We provide just the selected choice to ensure Shiny can set the value
        select_org(first_species_id)

        # Set initial selected species name
        selected_species_name(first_species_name)
      }
    })

    # Pop-up modal for gene assembl information ----
    observeEvent(input$genome_assembl_button, {
      # Create species count summary by source
      org_info_df <- idep_data$org_info
      source_counts <- table(org_info_df$group)

      # Categorize sources according to requirements
      ensembl_sources <- source_counts[!grepl("^Custom|^STRING", names(source_counts))]
      custom_count <- sum(source_counts[grepl("^Custom", names(source_counts))])
      string_count <- sum(source_counts[grepl("^STRING", names(source_counts))])

      # Create Ensembl breakdown
      ensembl_sources <- ensembl_sources[ensembl_sources > 1] # remove singletons 'Bacteria'
      ix <- which(names(ensembl_sources) == "ENSEMBL")
      names(ensembl_sources)[ix] <- "main"

      names(ensembl_sources) <- gsub("Ensembl", "", names(ensembl_sources))
      ensembl_parts <- paste0(names(ensembl_sources), "(", ensembl_sources, ")", collapse = ", ")

      # Create count summary text
      all_parts <- character(0)
      if (length(ensembl_parts) > 0) {
        all_parts <- c(all_parts, paste("Ensembl:", paste(ensembl_parts, collapse = ", ")))
      }
      if (string_count > 0) {
        all_parts <- c(all_parts, paste0("STRING-db (", string_count, ")"))
      }
      if (custom_count > 0) {
        all_parts <- c(all_parts, paste0("Custom (", custom_count, ")."))
      }

      count_text <- paste(paste(all_parts, collapse = "; "), " Use STRING-db annotations as a last resort.")

      shiny::showModal(
        shiny::modalDialog(
          size = "l",
          h3("Click on a species to select"),
          p(count_text, style = "font-style: italic; color: #666;"),
          easyClose = TRUE,
          DT::renderDataTable({
            df <- idep_data$org_info[
              ,
              c("ensembl_dataset", "name", "academicName", "taxon_id", "group")
            ]
            colnames(df) <- c(
              "Ensembl/STRING-db ID",
              "Name (Assembly)",
              "Academic Name",
              "Taxonomy ID",
              "Source"
            )
            row.names(df) <- NULL
            DT::datatable(
              df,
              selection = "single",
              options = list(
                lengthChange = FALSE,
                pageLength = 10,
                scrollY = "400px",
                columnDefs = list(list(visible = FALSE, targets = 0)),
                dom = "frtip",
                initComplete = DT::JS(
                  "function(settings, json) {",
                  "  $('.dataTables_filter input')",
                  "    .css('width', '300px')",
                  "    .attr('placeholder', 'Human, Homo sapiens, 9606');",
                  "}"
                )
              ),
              callback = DT::JS(
                paste0(
                  "table.on('click', 'tr', function() {
                    var data = table.row(this).data();
                    if (data) {
                      Shiny.setInputValue('", id, "-clicked_row', data[0]);
                    }
                  });"
                )
              ),
              rownames = FALSE
            )
          })
        )
      )
    })

    selected_species_name <- reactiveVal()

    # Handle new species checkbox ----
    observeEvent(input$new_species, {
      if (input$new_species) {
        select_org("NEW")
        updateCheckboxInput(
          session = session,
          inputId = "no_id_conversion",
          label = "Do not convert gene IDs",
          value = TRUE
        )
        selected_species_name("NEW")
      } else {
        # Reset to first species when unchecked
        first_species_id <- idep_data$org_info$id[1]
        first_species_name <- idep_data$org_info$name2[1]

        select_org(first_species_id)
        updateCheckboxInput(
          session = session,
          inputId = "no_id_conversion",
          label = "Do not convert gene IDs",
          value = FALSE
        )
        selected_species_name(first_species_name)
      }
    })

    # Handle GMT file upload ----
    observeEvent(input$gmt_file, {
      req(!is.null(input$gmt_file))
      # Automatically check the new species checkbox when GMT is uploaded
      updateCheckboxInput(
        session = session,
        inputId = "new_species",
        value = TRUE
      )
    })

    output$show_gmt <- renderUI({
      req(!is.null(input$gmt_file))

      in_file <- input$gmt_file
      in_file <- in_file$datapath
      lines <- scan(in_file, what = "", sep = "\n")
      total <- length(lines)
      if (length(lines) > 3) {
        lines <- lines[1:3]
      }
      tagList(
        h4("The uploaded GMT file has ", total, " gene-sets.
        Data will be analyzed in a custom mode, not using our database.
        The gene IDs in the GMT file must match those in the expression data."),
        br(),
        tags$pre(
          style = "font-size: 12px;",
          lines[1],
          if (length(lines) > 1) lines[2] else "",
          if (length(lines) > 2) lines[3] else ""
        )
      )
    })
    observeEvent(input$clicked_row, {
      # find species ID from ensembl_dataset
      selected_id <- find_species_id_by_ensembl(
        input$clicked_row,
        idep_data$org_info
      )
      selected_name <- find_species_by_id_name(selected_id, idep_data$org_info)

      # Update species selection with just the selected choice
      select_org(selected_id)

      # Uncheck new species when selecting from database
      updateCheckboxInput(
        session = session,
        inputId = "new_species",
        value = FALSE
      )

      # Reset no_id_conversion to FALSE for database species
      updateCheckboxInput(
        session = session,
        inputId = "no_id_conversion",
        label = "Do not convert gene IDs",
        value = FALSE
      )

      selected_species_name(selected_name)

      # Close the modal after selection
      removeModal()
    })

    output$selected_species <- renderText({
      selected_species_name()
    })

    # Default species message when user hasn't clicked Select ----------
    output$default_species_message <- renderUI({
      req(!is.null(loaded_data()$data))
      req(!is.null(idep_data$org_info) && nrow(idep_data$org_info) > 0)
      req(input$genome_assembl_button == 0 || is.null(input$genome_assembl_button))
      # Check if user is still using the first species (default)
      first_species_id <- idep_data$org_info$id[1]
      first_species_name <- idep_data$org_info$name2[1]

      # Only show message if using default species and not in NEW mode
      if (select_org() == first_species_id && select_org() != "NEW") {
        div(
          style = "background-color: #d1ecf1; border: 1px solid #bee5eb;
                   border-radius: 4px; padding: 10px; margin: 10px 0;
                   color: #0c5460;",
          h5(
            icon("info-circle"),
            paste0(" Using ", first_species_name, " genome annotations and pathways.")
          )
        )
      }
    })

    # Gene ID conversion statistics message ----------
    output$conversion_stats_message <- renderUI({
      req(!is.null(conversion_info()$converted))
      req(select_org() != "NEW")
      req(!is.null(loaded_data()$data))

      original_count <- length(conversion_info()$converted$origninal_ids)
      converted_count <- length(conversion_info()$converted$ids)
      conversion_rate <- round((converted_count / original_count) * 100, 1)

      # Determine styling and icon based on conversion rate
      if (conversion_rate < 80) {
        bg_color <- "#fff3cd"
        border_color <- "#ffeaa7"
        text_color <- "#856404"
        icon_name <- "exclamation-triangle"
        message_type <- "Warning: Low conversion rate!"
      } else {
        bg_color <- "#d4edda"
        border_color <- "#c3e6cb"
        text_color <- "#155724"
        icon_name <- "check-circle"
        message_type <- "Gene ID conversion:"
      }

      div(
        style = paste0(
          "background-color: ", bg_color, "; border: 1px solid ",
          border_color, "; border-radius: 4px; padding: 10px; margin: 10px 0;
                       color: ", text_color, ";"
        ),
        h5(
          icon(icon_name),
          paste0(
            " ", message_type, " ", converted_count, " out of ",
            original_count, " genes (", conversion_rate,
            "%) converted to Ensembl/STRING IDs."
          )
        )
      )
    })

    observeEvent(input$data_file_format,
      {
        req(input$data_file_format != 0)

        # Update dropdown text
        updateSelectInput(
          session = session,
          inputId = "data_file_format",
          choices = list(
            "..." = 0,
            "Read counts data" = 1,
            "Normalized expression data" = 2,
            "Fold changes & adjusted p-values" = 3,
            "Fold changes only" = 4
          ),
          selected = input$data_file_format
        )

        # Show notification only when on Data tab
        req(tab() == "Data")

        details <- get_data_type_details(input$data_file_format)
        preview <- demo_preview_tables[[as.character(input$data_file_format)]]

        demo_preview_content(preview)

        # Build notification UI with details and preview
        notification_ui <- if (!is.null(preview)) {
          table_html <- build_preview_table(preview)
          tagList(
            tags$strong(details$title),
            tags$br(),
            details$body,
            div(
              style = "max-height: 300px; overflow-y: auto; margin-top: 10px;",
              table_html
            )
          )
        } else {
          tagList(
            tags$strong(details$title),
            tags$br(),
            details$body
          )
        }

        showNotification(
          ui = notification_ui,
          duration = 20,
          type = "message",
          id = "load_prompt"
        )
      },
      ignoreNULL = TRUE
    )

    # Dismiss notification when user interacts with file input or demo button
    observeEvent(input$expression_file, {
      removeNotification("load_prompt")
    })

    observeEvent(input$demo_modal_button, {
      removeNotification("load_prompt")
    })

    observe({
      req(tab() == "Data")
      req(input$data_file_format == 0)

      removeNotification("load_prompt")
      demo_preview_content(NULL)
      removeModal()
    })

    observeEvent(tab(), {
      if (tab() != "Data") {
        removeNotification("load_prompt")
        demo_preview_content(NULL)
        removeModal()
      }
    })

    # Available demo data files for current data format ----
    demo_choices <- reactive({
      req(input$data_file_format)

      if (input$data_file_format <= 0) {
        return(NULL)
      }

      files <- idep_data$demo_file_info
      files <- files[files$type == input$data_file_format, , drop = FALSE]

      selected_species <- selected_species_name()
      if (is.null(selected_species) || !nzchar(selected_species)) {
        return(NULL)
      }
      if (identical(selected_species, "NEW")) {
        return(NULL)
      }
      has_species_column <- "species" %in% names(files)

      if (
        has_species_column &&
          !is.null(selected_species) &&
          nzchar(selected_species)
      ) {
        normalize_species <- function(x) {
          tolower(trimws(as.character(x)))
        }

        species_match <- normalize_species(files$species) ==
          normalize_species(selected_species)

        files <- files[species_match, , drop = FALSE]
      }

      if (!nrow(files)) {
        return(NULL)
      }

      setNames(as.list(files$ID), files$name)
    })

    observeEvent(demo_choices(),
      {
        choices <- demo_choices()
        values <- unlist(choices, use.names = FALSE)

        if (length(values) == 0) {
          selected_demo(NULL)
        } else {
          current <- selected_demo()
          if (is.null(current) || !(current %in% values)) {
            selected_demo(values[[1]])
          }
        }
      },
      ignoreNULL = FALSE
    )

    observeEvent(input$select_demo,
      {
        req(!is.null(input$select_demo))
        selected_demo(input$select_demo)
      },
      ignoreNULL = TRUE
    )

    # UI elements for load demo action button, demo data drop down, and -----
    # expression file upload
    output$load_data_ui <- renderUI({
      req(go_button_count() == 0)
      req(input$data_file_format)

      choices <- demo_choices()
      has_demo_datasets <- !is.null(choices) && length(choices) > 0

      tagList(
        strong("3. Expression matrix (CSV or text)"),
        fluidRow(
          column(
            width = 6,
            # Expression data file input
            div(
              id = ns("expression_file_container"),
              style = "cursor: help;",
              fileInput(
                inputId = ns("expression_file"),
                label = NULL,
                accept = c(
                  "text/csv",
                  "text/comma-separated-values",
                  "text/tab-separated-values",
                  "text/plain",
                  ".csv",
                  ".tsv",
                  ".xlsx",
                  ".xls"
                ),
                placeholder = "CSV or text",
                width = "100%"
              )
            ),
            tippy::tippy_this(
              ns("expression_file_container"),
              "Upload your expression data. CSV files and tab- or space- delimited text files are recommended. Excel files are also supported.",
              theme = "light"
            )
          ),
          if (has_demo_datasets) {
            column(
              width = 4,
              actionButton(
                inputId = ns("demo_modal_button"),
                label = tags$span("Demo Data", style = "color: red;"),
                class = "btn-default"
              ),
              tippy::tippy_this(
                ns("demo_modal_button"),
                "Load demo datasets.",
                theme = "light"
              )
            )
          },
          column(
            width = 2,
            align = "center",
            actionLink(
              inputId = ns("data_format_help"),
              label = NULL,
              icon = icon("info-circle"),
              class = "load-data-info-icon"
            ),
            tippy::tippy_this(
              ns("data_format_help"),
              "Learn more about the data types iDEP accepts.",
              theme = "light"
            )
          )
        )
      )
    })

    observeEvent(input$demo_modal_button, {
      choices <- demo_choices()
      req(!is.null(choices), length(choices) > 0)

      selected_choice <- isolate({
        current <- selected_demo()
        values <- unlist(choices, use.names = FALSE)
        if (!is.null(current) && current %in% values) {
          current
        } else if (length(values) > 0) {
          values[[1]]
        } else {
          NULL
        }
      })

      shiny::showModal(
        shiny::modalDialog(
          title = "Demo Datasets",
          easyClose = TRUE,
          size = "s",
          footer = NULL,
          tagList(
            selectInput(
              inputId = ns("select_demo"),
              label = NULL,
              choices = choices,
              selected = selected_choice,
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("select_demo"),
              "Pick a demo dataset, then click Load to analyze it.",
              theme = "light"
            ),
            div(
              style = "margin-top: 8px;",
              uiOutput(ns("demo_memo"))
            ),
            br(),
            fluidRow(
              column(
                width = 4,
                align = "left",
                actionButton(
                  inputId = ns("go_button"),
                  label = "Load",
                  class = "btn-primary"
                ),
                tippy::tippy_this(
                  ns("go_button"),
                  "Load the selected dataset.",
                  theme = "light"
                )
              ),
              column(
                width = 8,
                align = "left",
                div(
                  class = "demo-download-actions",
                  style = "display: flex; gap: 10px; flex-wrap: wrap; justify-content: flex-end; align-items: center;",
                  tagList(
                    downloadButton(
                      outputId = ns("download_demo_expression"),
                      label = NULL,
                      class = "btn-default"
                    ),
                    tippy::tippy_this(
                      ns("download_demo_expression"),
                      "Download the selected demo dataset.",
                      theme = "light"
                    )
                  ),
                  uiOutput(ns("demo_design_download_ui"))
                )
              )
            )
          )
        )
      )
    })

    observeEvent(input$go_button, {
      shiny::removeModal()
    })

    # Alternate ui message to reset app once data is loaded ----
    output$load_data_alt <- renderUI({
      if (go_button_count() == 0 && is.null(input$expression_file)) {
        tagList(
          tags$span("Quick Start:", style = "font-size: 18px;"),
          tags$ul(
            tags$li(
              "Watch a 5-min ",
              a("video", href = "https://youtu.be/lqDqrJU-e24", target = "_blank"),
              " tutorial!"
            ),
            tags$li(
              "Select a data type, then click ",
              tags$span("Demo Data.", id = "load-demo", style = "color: #B30000;")
            )
          )
        )
      }
    })


    output$reset_button <- renderUI({
      if (go_button_count() == 0 && is.null(input$expression_file)) {
        NULL
      } else {
        # reset message and action button
        tagList(
          fluidRow(
            column(
              width = 12,
              align = "center",
              actionButton(
                inputId = ns("reset_app_new_data"),
                label = strong(tags$span("Reset", style = "color: red;")),
                align = "left"
              ),
              tippy::tippy_this(
                ns("reset_app_new_data"),
                "Clear the current data and start a fresh analysis.",
                theme = "light"
              )
            )
          ),
          br(),
          br()
        )
      }
    })

    observeEvent(input$reset_app_new_data, {
      session$reload()
    })

    # UI element for design file upload ----
    output$design_file_ui <- renderUI({
      req(go_button_count() == 0)

      tagList(
        tags$style(HTML(paste0(
          "#", ns("experiment_file_container"), " .form-group { margin-bottom: 4px; }",
          "#", ns("experiment_file_progress"), " { display: none; }"
        ))),
        strong("4. Optional: Experiment Design (CSV or text)"),
        fluidRow(
          column(
            width = 10,
            div(
              id = ns("experiment_file_container"),
              style = "cursor: help;",
              fileInput(
                inputId = ns("experiment_file"),
                label = NULL,
                accept = c(
                  "text/csv",
                  "text/comma-separated-values",
                  "text/tab-separated-values",
                  "text/plain",
                  ".csv",
                  ".tsv",
                  ".xlsx",
                  ".xls"
                ),
                placeholder = "CSV or text"
              )
            ),
            tippy::tippy_this(
              ns("experiment_file_container"),
              "Upload a sample information table to define experimental design. CSV files and tab- or space- delimited text files are recommended. Excel files are also supported.",
              theme = "light"
            )
          ),
          column(
            width = 2,
            align = "center",
            actionLink(
              inputId = ns("design_format_help"),
              label = NULL,
              icon = icon("info-circle"),
              class = "load-data-info-icon"
            ),
            tippy::tippy_this(
              ns("design_format_help"),
              "Learn more about experimental design file format.",
              theme = "light"
            )
          )
        )
      )
    })

    # Interactive design builder button (no go_button_count guard — must render for demo data too)
    output$wizard_btn_ui <- renderUI({
      req(!is.null(loaded_data()$data))
      req(is.null(input$experiment_file))
      div(
        style = "margin-top: 4px;margin-bottom: 9px;",
        actionButton(
          inputId = ns("build_design_btn"),
          label = "+ Create design",
          class = "btn-xs btn-default"
        ),
        if (!is.null(gui_design())) {
          tags$span(
            icon("check-circle"),
            " Design applied",
            style = "color: #2e7d32; font-size: 12px; margin-left: 8px;"
          )
        }
      )
    })

    # Helper: build one guided factor row div (used on open and on insertUI)
    make_guided_factor_row <- function(i, default_name = "", default_levels = "") {
      div(
        id = ns(paste0("guided_factor_row_", i)),
        class = "well well-sm",
        style = "margin-bottom: 8px; padding: 10px;",
        fluidRow(
          column(
            4,
            textInput(
              ns(paste0("gf_name_", i)),
              label = paste("Factor", i, "name:"),
              value = default_name,
              placeholder = "e.g. Treatment"
            )
          ),
          column(
            5,
            div(
              id = ns(paste0("gf_levels_container_", i)),
              textInput(
                ns(paste0("gf_levels_", i)),
                label = "Levels (comma-separated):",
                value = default_levels,
                placeholder = "e.g. ctrl, treated"
              )
            ),
            div(
              id = ns(paste0("gf_blocking_note_", i)),
              style = "display: none;",
              tags$small(
                style = "color: #666; margin-top: 5px; display: block;",
                shiny::icon("info-circle"),
                " Pair indices assigned per sample in the next step."
              )
            )
          ),
          column(
            3,
            tags$label("Options:"),
            tags$br(),
            checkboxInput(
              ns(paste0("gf_blocking_", i)),
              label = "Paired / blocking",
              value = FALSE
            )
          )
        )
      )
    }

    # Open design builder modal ----
    observeEvent(input$build_design_btn, {
      req(!is.null(loaded_data()$data))
      sample_names <- colnames(loaded_data()$data)

      # If a design was previously applied, pre-fill the wizard from it so the
      # user can tweak a cell or download the CSV without rebuilding. Treat all
      # factors as non-blocking — blocking metadata is not stored in the
      # applied matrix.
      applied_design <- gui_design()
      has_applied <- !is.null(applied_design) && ncol(applied_design) > 0L

      if (has_applied) {
        applied_factor_defs <- lapply(seq_len(ncol(applied_design)), function(i) {
          col_values <- as.character(applied_design[, i])
          list(
            name = colnames(applied_design)[i],
            levels = unique(col_values),
            is_blocking = FALSE
          )
        })
        applied_assignments <- setNames(
          lapply(seq_along(applied_factor_defs), function(i) {
            as.character(applied_design[, i])
          }),
          as.character(seq_along(applied_factor_defs))
        )
        n_guided_factors(length(applied_factor_defs))
        guided_factor_defs(applied_factor_defs)
        sample_assignments(applied_assignments)

        factor_rows_init <- lapply(
          seq_along(applied_factor_defs),
          function(i) {
            make_guided_factor_row(
              i,
              default_name = applied_factor_defs[[i]]$name,
              default_levels = paste(
                applied_factor_defs[[i]]$levels,
                collapse = ", "
              )
            )
          }
        )
      } else {
        n_guided_factors(2L)
        guided_factor_defs(NULL)
        sample_assignments(list())
        initial_groups <- detect_groups(sample_names)
        level_suggestion <- paste(unique(initial_groups), collapse = ", ")
        factor_rows_init <- list(
          make_guided_factor_row(1L, "group", level_suggestion),
          make_guided_factor_row(2L)
        )
      }

      shiny::showModal(shiny::modalDialog(
        title = "Experiment Design",
        size = "l",
        tabsetPanel(
          id = ns("design_builder_tab"),
          tabPanel(
            "Guided",
            br(),
            # Step 1: define factors and levels
            div(
              id = ns("guided_step1"),
              style = "overflow-x: hidden;",
              tags$p(
                style = "color: #555; font-size: 13px; margin-bottom: 12px;",
                "Name each factor and its levels. Up to 5 factors can be defined. Check 'Paired/blocking' for pair/replicate index factors."
              ),
              div(
                id = ns("guided_factor_rows_container"),
                factor_rows_init
              ),
              hr(),
              fluidRow(
                column(
                  6,
                  actionButton(
                    ns("add_guided_factor"), "+ Factor",
                    class = "btn-xs btn-default"
                  ),
                  actionButton(
                    ns("remove_guided_factor"), "- Factor",
                    class = "btn-xs btn-default",
                    style = "margin-left: 6px;"
                  )
                ),
                column(
                  6, align = "right",
                  actionButton(
                    ns("guided_next_btn"), "Next \u2192",
                    class = "btn-primary btn-sm"
                  )
                )
              )
            ),
            # Step 2: assign samples (hidden until Next clicked)
            shinyjs::hidden(
              div(
                id = ns("guided_step2"),
                style = "max-height: 60vh; overflow-y: auto; overflow-x: hidden;",
                uiOutput(ns("guided_step2_ui")),
                hr(),
                fluidRow(
                  column(
                    6,
                    actionButton(
                      ns("guided_back_btn"), "← Back",
                      class = "btn-default btn-sm"
                    )
                  ),
                  column(
                    6, align = "right",
                    actionButton(
                      ns("apply_guided_design"), "Apply",
                      class = "btn-primary btn-sm"
                    )
                  )
                )
              )
            )
          ),
          tabPanel(
            "Manual",
            br(),
            fluidRow(
              column(
                9,
                tags$p(
                  style = "color: #555; font-size: 13px; margin-bottom: 8px;",
                  "Rows are factors, columns are samples. Type a label per cell. Use consistent labels within each row. Row 1 is automatically filled based on sample names, but can be edited."
                )
              ),
              column(
                3, align = "right",

                tags$div(
                  id = ns("download_design_template_tip"),
                  style = "display: inline-block;",
                  downloadButton(
                    ns("download_design_template"), " ",
                    class = "btn-sm btn-default"
                  )
                ),

                tippy::tippy_this(
                  ns("download_design_template_tip"),
                  "Download CSV",
                  theme = "light-border",
                  placement = "top"
                )
              )
            ),
            tags$style(HTML(
              ".design-grid-table .form-group { margin-bottom: 2px; }
               .design-grid-table input.form-control { font-size: 12px; padding: 2px 4px; height: 26px; }"
            )),
            div(style = "overflow-x: auto;", uiOutput(ns("design_grid_ui"))),
            hr(),
            div(
              style = "text-align: right;",
              shiny::actionButton(
                ns("apply_design"), "Apply",
                class = "btn-primary btn-sm"
              )
            )
          )
        ),
        easyClose = FALSE,
        footer = shiny::modalButton("Cancel")
      ))

      # Inject JS helper so onclick attributes in chip spans can call assignChip()
      shinyjs::runjs(sprintf(
        "window.IDEP_DESIGN_NS = '%s';
         window.assignChip = function(fi, j) {
           var checked = $(\"input[name='\" + IDEP_DESIGN_NS + \"active_level_\" + fi + \"']:checked\");
           if (!checked.length) return;
           Shiny.setInputValue(IDEP_DESIGN_NS + 'chip_clicked',
             {factor: fi, sample: j, level: checked.val()},
             {priority: 'event'}
           );
         };",
        ns("")
      ))

      # If a design is already applied, skip directly to step 2 so the user
      # sees chips colored per the current assignments.
      if (has_applied) {
        shinyjs::hide("guided_step1")
        shinyjs::show("guided_step2")
      }
    })

    # Add a factor row to the guided builder ----
    observeEvent(input$add_guided_factor, {
      current_n <- n_guided_factors()
      if (current_n >= 5L) return()
      new_i <- current_n + 1L
      n_guided_factors(new_i)
      insertUI(
        selector = paste0("#", ns("guided_factor_rows_container")),
        where = "beforeEnd",
        ui = make_guided_factor_row(new_i),
        immediate = TRUE
      )
      if (new_i >= 5L) shinyjs::disable("add_guided_factor")
      shinyjs::enable("remove_guided_factor")
    })

    # Remove the last factor row from the guided builder ----
    observeEvent(input$remove_guided_factor, {
      current_n <- n_guided_factors()
      if (current_n <= 1L) return()
      removeUI(
        selector = paste0("#", ns(paste0("guided_factor_row_", current_n))),
        immediate = TRUE
      )
      n_guided_factors(current_n - 1L)
      if (current_n - 1L <= 1L) shinyjs::disable("remove_guided_factor")
      shinyjs::enable("add_guided_factor")
    })

    # Toggle levels input vs. blocking note for each factor slot ----
    lapply(seq_len(5L), function(i) {
      observeEvent(input[[paste0("gf_blocking_", i)]], {
        if (isTRUE(input[[paste0("gf_blocking_", i)]])) {
          shinyjs::hide(paste0("gf_levels_container_", i))
          shinyjs::show(paste0("gf_blocking_note_", i))
        } else {
          shinyjs::show(paste0("gf_levels_container_", i))
          shinyjs::hide(paste0("gf_blocking_note_", i))
        }
      }, ignoreNULL = TRUE, ignoreInit = TRUE)
    })

    # Validate step 1 and advance to step 2 ----
    observeEvent(input$guided_next_btn, {
      n <- n_guided_factors()
      factor_defs <- list()
      errors <- character(0)
      typo_warnings <- character(0)
      non_ascii_seen <- FALSE

      for (i in seq_len(n)) {
        fname_raw <- input[[paste0("gf_name_", i)]]
        # Strict sanitizer (matches build_sample_info_from_df): keep [A-Za-z0-9_]
        fname <- sanitize_created_names(
          trimws(if (is.null(fname_raw)) "" else fname_raw)
        )

        is_blocking <- isTRUE(input[[paste0("gf_blocking_", i)]])

        if (nchar(fname) == 0L) {
          errors <- c(errors, paste0("Factor ", i, " has no name."))
          next
        }

        if (is_blocking) {
          factor_defs[[length(factor_defs) + 1L]] <- list(
            name = fname, levels = NULL, is_blocking = TRUE
          )
        } else {
          raw_levels_input <- input[[paste0("gf_levels_", i)]]
          raw_levels <- if (is.null(raw_levels_input)) "" else raw_levels_input
          levels_vec <- trimws(strsplit(raw_levels, ",", fixed = TRUE)[[1L]])
          # Strict sanitizer (matches build_sample_info_from_df cell treatment)
          levels_vec <- toupper(sanitize_created_names(levels_vec))
          levels_vec <- levels_vec[nchar(levels_vec) > 0L]
          # Deduplicate after sanitization (e.g. "ctrl" and "Ctrl" both become "CTRL")
          levels_vec <- unique(levels_vec)
          if (length(levels_vec) < 2L) {
            errors <- c(errors, paste0(
              "Factor '", fname, "' needs at least 2 distinct levels ",
              "(after removing spaces and uppercasing)."
            ))
            next
          }
          near_dups <- near_duplicate_pairs(levels_vec)
          if (length(near_dups) > 0L) {
            typo_warnings <- c(typo_warnings, paste0(
              "Factor '", fname, "': possible typo — ",
              paste(near_dups, collapse = ", ")
            ))
          }
          if (has_non_ascii(c(fname, levels_vec))) non_ascii_seen <- TRUE
          factor_defs[[length(factor_defs) + 1L]] <- list(
            name = fname, levels = levels_vec, is_blocking = FALSE
          )
        }
      }

      if (length(errors) > 0L) {
        shiny::showNotification(
          paste(errors, collapse = " "),
          type = "warning", duration = 6L
        )
        return()
      }
      if (length(factor_defs) == 0L) {
        shiny::showNotification(
          "Define at least one factor before continuing.",
          type = "warning", duration = 5L
        )
        return()
      }

      if (length(typo_warnings) > 0L) {
        shiny::showNotification(
          paste(typo_warnings, collapse = " | "),
          type = "warning", duration = 7L
        )
      }
      if (non_ascii_seen) {
        shiny::showNotification(
          "Non-ASCII characters detected — downstream tools may rewrite them.",
          type = "warning", duration = 6L
        )
      }

      guided_factor_defs(factor_defs)
      sample_assignments(list())
      shinyjs::hide("guided_step1")
      shinyjs::show("guided_step2")
    })

    # Return from step 2 to step 1 ----
    observeEvent(input$guided_back_btn, {
      shinyjs::show("guided_step1")
      shinyjs::hide("guided_step2")
    })

    # Colors for chip backgrounds (one per level slot, cycled if >5 levels) ----
    .guided_level_bg  <- c('#d0e9ff', '#d4edda', '#fff3cd', '#f8d7da', '#e8d5f5')
    .guided_level_bdr <- c('#9ecfff', '#a5d6a7', '#ffe082', '#f5c6cb', '#ce93d8')

    # Render sample-assignment UI for step 2 ----
    output$guided_step2_ui <- renderUI({
      fdefs <- guided_factor_defs()
      req(!is.null(fdefs), !is.null(loaded_data()$data))
      sample_names <- colnames(loaded_data()$data)

      tagList(
        fluidRow(
          column(
            9,
            tags$p(
              style = 'color: #555; font-size: 13px; margin-bottom: 12px;',
              'Pick a level, then click samples to assign.'
            )
          ),
          column(
            3, align = "right",

            tags$div(
              id = ns("download_guided_design_tip"),
              style = "display: inline-block;",
              downloadButton(
                ns("download_guided_design"), " ",
                class = "btn-sm btn-default"
              )
            ),

            tippy::tippy_this(
              ns("download_guided_design_tip"),
              "Download CSV",
              theme = "light-border",
              placement = "top"
            )
          )
        ),
        lapply(seq_along(fdefs), function(fi) {
          fdef <- fdefs[[fi]]
          if (fdef$is_blocking) {
            sample_inputs <- lapply(seq_along(sample_names), function(j) {
              div(
                style = paste0(
                  'display: inline-block; margin: 4px; ',
                  'text-align: center; vertical-align: top; min-width: 80px;'
                ),
                tags$div(
                  style = 'font-size: 11px; color: #555; word-break: break-all;',
                  sample_names[j]
                ),
                numericInput(
                  ns(paste0('gs_', fi, '_', j)),
                  label = NULL,
                  value = j,
                  min = 1L, max = length(sample_names), step = 1L,
                  width = '70px'
                )
              )
            })
            div(
              class = 'well well-sm',
              style = 'margin-bottom: 12px; padding: 10px;',
              tags$h5(
                style = 'margin-top: 0; margin-bottom: 6px;',
                tags$strong(fdef$name),
                tags$span(
                  class = 'label label-info',
                  style = 'margin-left: 6px; font-size: 11px;',
                  'Paired / Blocking'
                )
              ),
              tags$p(
                style = 'font-size: 12px; color: #666; margin-bottom: 8px;',
                'Same index = same pair.'
              ),
              div(style = 'overflow-x: auto;', do.call(tagList, sample_inputs))
            )
          } else {
            div(
              class = 'well well-sm',
              style = 'margin-bottom: 12px; padding: 10px;',
              tags$h5(
                style = 'margin-top: 0; margin-bottom: 6px;',
                tags$strong(fdef$name)
              ),
              fluidRow(
                column(
                  4,
                  tags$p(
                    style = 'font-size: 12px; color: #555; margin-bottom: 6px;',
                    'Active level:'
                  ),
                  radioButtons(
                    ns(paste0('active_level_', fi)),
                    label = NULL,
                    choiceValues = fdef$levels,
                    choiceNames  = lapply(seq_along(fdef$levels), function(k) {
                      tags$span(
                        tags$span(
                          style = paste0(
                            'display:inline-block; width:12px; height:12px; ',
                            'border-radius:2px; margin-right:6px; ',
                            'vertical-align:middle; border:1px solid; ',
                            'background:', .guided_level_bg[((k - 1L) %% 5L) + 1L], ';',
                            'border-color:', .guided_level_bdr[((k - 1L) %% 5L) + 1L], ';'
                          )
                        ),
                        fdef$levels[k]
                      )
                    }),
                    selected = fdef$levels[1L]
                  ),
                  if (length(fdef$levels) == 2L) {
                    actionButton(
                      ns(paste0('assign_rest_', fi)),
                      'Assign remaining samples \u2192',
                      class = 'btn-xs btn-info'
                    )
                  }
                ),
                column(
                  8,
                  tags$p(
                    style = 'font-size: 12px; color: #555; margin-bottom: 4px;',
                    paste0('Levels: ', paste(fdef$levels, collapse = ', '))
                  ),
                  div(
                    style = paste0(
                      'max-height: 160px; overflow-y: auto; ',
                      'border: 1px solid #ddd; border-radius: 4px; ',
                      'padding: 8px; background: #fafafa;'
                    ),
                    uiOutput(ns(paste0('guided_chips_', fi)))
                  ),
                  uiOutput(ns(paste0('guided_assign_count_', fi)))
                )
              )
            )
          }
        })
      )
    })

    # Chip renderers - re-render sample chips when assignments change ----
    lapply(seq_len(5L), function(fi) {
      output[[paste0('guided_chips_', fi)]] <- renderUI({
        fdefs <- guided_factor_defs()
        if (is.null(fdefs) || fi > length(fdefs) || fdefs[[fi]]$is_blocking) {
          return(NULL)
        }
        req(!is.null(loaded_data()$data))
        sample_names <- colnames(loaded_data()$data)
        fdef         <- fdefs[[fi]]
        assignments  <- sample_assignments()
        fi_asgn      <- assignments[[as.character(fi)]]
        if (is.null(fi_asgn)) fi_asgn <- rep('', length(sample_names))

        lapply(seq_along(sample_names), function(j) {
          assigned <- fi_asgn[j]
          lvl_idx  <- match(assigned, fdef$levels)
          if (!is.na(lvl_idx) && nchar(assigned) > 0L) {
            k   <- ((lvl_idx - 1L) %% 5L) + 1L
            sty <- paste0(
              'background:', .guided_level_bg[k], ';',
              'border-color:', .guided_level_bdr[k], ';'
            )
            tip <- paste0('Assigned to: ', assigned)
          } else {
            sty <- 'background:#f5f5f5; border-color:#ccc;'
            tip <- 'Unassigned - select a level then click'
          }
          tags$span(
            sample_names[j],
            onclick = sprintf('window.assignChip(%d,%d)', fi, j),
            title   = tip,
            style   = paste0(
              'display:inline-block; margin:2px; padding:4px 10px; ',
              'border-radius:12px; cursor:pointer; font-size:12px; ',
              'border:1px solid; user-select:none; ',
              sty
            )
          )
        })
      })
    })

    # Live assignment counters per factor slot ----
    lapply(seq_len(5L), function(fi) {
      output[[paste0('guided_assign_count_', fi)]] <- renderUI({
        fdefs <- guided_factor_defs()
        if (is.null(fdefs) || fi > length(fdefs) || fdefs[[fi]]$is_blocking) {
          return(NULL)
        }
        req(!is.null(loaded_data()$data))
        sample_names <- colnames(loaded_data()$data)
        assignments  <- sample_assignments()
        fi_asgn      <- assignments[[as.character(fi)]]
        if (is.null(fi_asgn)) fi_asgn <- rep('', length(sample_names))
        assigned <- sum(nchar(fi_asgn) > 0L)
        color    <- if (assigned == length(sample_names)) '#2e7d32' else '#e65100'
        tags$small(
          style = paste0(
            'color: ', color, '; font-weight: bold; ',
            'margin-top: 4px; display: block;'
          ),
          shiny::icon(
            if (assigned == length(sample_names)) 'check-circle' else 'exclamation-circle'
          ),
          paste0(' ', assigned, ' / ', length(sample_names), ' samples assigned')
        )
      })
    })

    # Single chip-click handler (one observer handles all chips) ----
    observeEvent(input$chip_clicked, {
      req(!is.null(loaded_data()$data))
      fi    <- input$chip_clicked$factor
      j     <- input$chip_clicked$sample
      level <- input$chip_clicked$level
      sample_names <- colnames(loaded_data()$data)
      fdefs <- guided_factor_defs()
      # Guard: JS indices could be stale if factors/data changed mid-interaction
      if (!is.numeric(fi) || fi < 1L ||
          is.null(fdefs) || fi > length(fdefs)) return()
      if (!is.numeric(j) || j < 1L || j > length(sample_names)) return()
      asgn  <- sample_assignments()
      key   <- as.character(fi)
      if (is.null(asgn[[key]])) asgn[[key]] <- rep('', length(sample_names))
      asgn[[key]][j] <- level
      sample_assignments(asgn)
    }, ignoreInit = TRUE)

    # 'Assign remaining' buttons (2-level factors only) ----
    lapply(seq_len(5L), function(fi) {
      observeEvent(input[[paste0('assign_rest_', fi)]], {
        fdefs <- guided_factor_defs()
        if (is.null(fdefs) || fi > length(fdefs)) return()
        fdef <- fdefs[[fi]]
        if (fdef$is_blocking || length(fdef$levels) != 2L) return()
        req(!is.null(loaded_data()$data))
        sample_names <- colnames(loaded_data()$data)
        active_level <- input[[paste0('active_level_', fi)]]
        other_level  <- setdiff(fdef$levels, active_level)[1L]
        asgn         <- sample_assignments()
        key          <- as.character(fi)
        current      <- if (is.null(asgn[[key]])) rep('', length(sample_names)) else asgn[[key]]
        current[current == ''] <- other_level
        asgn[[key]]  <- current
        sample_assignments(asgn)
      }, ignoreNULL = TRUE)
    })

    # Apply design built via the guided wizard ----
    observeEvent(input$apply_guided_design, {
      req(!is.null(loaded_data()$data))
      fdefs <- guided_factor_defs()
      req(!is.null(fdefs))
      sample_names <- colnames(loaded_data()$data)
      asgn  <- sample_assignments()

      rows      <- list()
      row_names <- character(0)
      errors    <- character(0)
      all_level_values <- character(0)

      for (fi in seq_along(fdefs)) {
        fdef <- fdefs[[fi]]
        if (fdef$is_blocking) {
          cell_values <- vapply(seq_along(sample_names), function(j) {
            v <- input[[paste0('gs_', fi, '_', j)]]
            if (is.null(v)) as.character(j) else as.character(v)
          }, character(1L))
        } else {
          fi_asgn     <- asgn[[as.character(fi)]]
          if (is.null(fi_asgn)) fi_asgn <- rep('', length(sample_names))
          cell_values  <- fi_asgn
          n_unassigned <- sum(nchar(cell_values) == 0L)
          if (n_unassigned > 0L) {
            errors <- c(errors, paste0(
              "Factor '", fdef$name, "': ", n_unassigned,
              " sample(s) not yet assigned."
            ))
          }
          all_level_values <- c(
            all_level_values, unique(cell_values[nchar(cell_values) > 0L])
          )
        }
        rows[[fi]]  <- cell_values
        row_names   <- c(row_names, fdef$name)
      }

      if (length(errors) > 0L) {
        shiny::showNotification(
          paste(errors, collapse = ' | '),
          type = 'error', duration = 8L
        )
        return()
      }

      level_counts <- table(tolower(all_level_values))
      dup_levels   <- names(level_counts[level_counts > 1L])
      if (length(dup_levels) > 0L) {
        shiny::showNotification(
          paste0(
            'Level names shared across factors (',
            paste(dup_levels, collapse = ', '),
            ') \u2014 iDEP will add factor prefixes to avoid confusion.'
          ),
          type = 'warning', duration = 5L
        )
      }

      if (has_non_ascii(c(row_names, unlist(rows)))) {
        shiny::showNotification(
          'Non-ASCII characters detected \u2014 downstream tools may rewrite them.',
          type = 'warning', duration = 6L
        )
      }

      design_df <- as.data.frame(
        do.call(rbind, rows),
        stringsAsFactors = FALSE
      )
      colnames(design_df) <- sample_names
      rownames(design_df) <- row_names

      si <- build_sample_info_from_df(
        design_df, sample_names, input$max_group_name_length
      )

      if (is.null(si)) {
        shiny::showNotification(
          paste0(
            'Design not applied \u2014 ensure each factor ',
            'has at least two distinct group labels.'
          ),
          type = 'error', duration = 6L
        )
        return()
      }

      gui_design(si)
      shiny::removeModal()
      shiny::showNotification(
        'Experimental design applied.', type = 'message', duration = 4L
      )
    })

    # Download guided design as CSV ----
    output$download_guided_design <- downloadHandler(
      filename = function() 'experimental_design.csv',
      content = function(file) {
        req(!is.null(loaded_data()$data))
        fdefs        <- guided_factor_defs()
        req(!is.null(fdefs))
        sample_names <- colnames(loaded_data()$data)
        asgn         <- sample_assignments()

        rows      <- list()
        row_names <- character(0)
        for (fi in seq_along(fdefs)) {
          fdef <- fdefs[[fi]]
          if (fdef$is_blocking) {
            cell_values <- vapply(seq_along(sample_names), function(j) {
              v <- input[[paste0('gs_', fi, '_', j)]]
              if (is.null(v)) as.character(j) else as.character(v)
            }, character(1L))
          } else {
            fi_asgn     <- asgn[[as.character(fi)]]
            cell_values <- if (is.null(fi_asgn)) rep('', length(sample_names)) else fi_asgn
          }
          rows[[fi]] <- cell_values
          row_names  <- c(row_names,
            if (nchar(fdef$name) > 0L) fdef$name else paste0('factor', fi))
        }
        row_names    <- make.unique(row_names)
        df           <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
        colnames(df) <- sample_names
        rownames(df) <- row_names
        write.csv(df, file)
      }
    )

    # Render the editable design grid inside the modal ----
    output$design_grid_ui <- renderUI({
      req(!is.null(loaded_data()$data))
      sample_names <- colnames(loaded_data()$data)
      n <- length(sample_names)
      n_factors <- 5L

      initial_groups <- detect_groups(sample_names)
      applied_design <- gui_design()
      n_applied <- if (is.null(applied_design)) 0L else ncol(applied_design)

      cell_px <- max(70L, min(110L, as.integer(660L / n)))
      cell_style <- paste0("min-width: ", cell_px, "px; padding: 2px;")

      header <- tags$tr(
        tags$th(
          "Factor Name",
          style = "min-width: 120px; background: #f0f0f0; font-size: 12px; padding: 4px;"
        ),
        lapply(sample_names, function(s) {
          tags$th(
            s,
            style = paste0(
              "background: #f0f0f0; font-size: 11px; padding: 4px; ",
              "word-break: break-all; ", cell_style
            )
          )
        })
      )

      factor_rows <- lapply(seq_len(n_factors), function(i) {
        # Row defaults priority:
        #   1. Pre-fill from applied design if this row is within it.
        #   2. Otherwise row 1 gets detect_groups guess; rest are blank.
        if (i <= n_applied) {
          default_name <- colnames(applied_design)[i]
          default_cells <- as.character(applied_design[, i])
        } else if (is.null(applied_design) && i == 1L) {
          default_name <- "group"
          default_cells <- initial_groups
        } else {
          default_name <- ""
          default_cells <- rep("", n)
        }

        name_td <- tags$td(
          style = "vertical-align: middle; padding: 2px;",
          textInput(ns(paste0("factor_name_", i)), label = NULL,
                    value = default_name, width = "100%")
        )
        cell_tds <- lapply(seq_len(n), function(j) {
          tags$td(
            style = cell_style,
            textInput(ns(paste0("cell_", i, "_", j)), label = NULL,
                      value = default_cells[[j]], width = "100%")
          )
        })
        do.call(tags$tr, c(list(name_td), cell_tds))
      })

      tagList(
        div(
          class = "design-grid-table",
          tags$table(
            class = "table table-bordered table-condensed",
            style = "margin-bottom: 4px; font-size: 12px;",
            tags$thead(header),
            tags$tbody(factor_rows)
          )
        )
      )
    })

    # Apply GUI-built design ----
    observeEvent(input$apply_design, {
      req(!is.null(loaded_data()$data))
      sample_names <- colnames(loaded_data()$data)
      n <- length(sample_names)
      n_factors <- 5L

      rows <- list()
      row_names <- character(0)

      for (i in seq_len(n_factors)) {
        factor_name <- input[[paste0("factor_name_", i)]]
        if (is.null(factor_name) || nchar(trimws(factor_name)) == 0L) next

        cell_values <- vapply(seq_len(n), function(j) {
          v <- input[[paste0("cell_", i, "_", j)]]
          if (is.null(v)) "" else v
        }, character(1L))

        rows[[length(rows) + 1L]] <- cell_values
        row_names <- c(row_names, trimws(factor_name))
      }

      if (length(rows) == 0L) {
        shiny::showNotification(
          "No factors defined. Fill in at least one factor name and group labels.",
          type = "warning", duration = 5
        )
        return()
      }

      if (anyDuplicated(row_names) > 0L) {
        shiny::showNotification(
          "Factor names must be unique. Please rename duplicate factors.",
          type = "warning", duration = 5
        )
        return()
      }

      # Block on near-duplicate cell labels within the same factor row.
      # Compare on the sanitized form (what becomes the design label).
      typo_errors <- character(0)
      for (k in seq_along(rows)) {
        sanitized <- toupper(sanitize_created_names(rows[[k]]))
        near_dups <- near_duplicate_pairs(sanitized)
        if (length(near_dups) > 0L) {
          typo_errors <- c(typo_errors, paste0(
            "Factor '", row_names[k], "': ",
            paste(near_dups, collapse = ", ")
          ))
        }
      }
      if (length(typo_errors) > 0L) {
        shiny::showNotification(
          paste0(
            "Possible typos — fix or confirm before applying: ",
            paste(typo_errors, collapse = " | ")
          ),
          type = "error", duration = 10
        )
        return()
      }

      if (has_non_ascii(c(row_names, unlist(rows)))) {
        shiny::showNotification(
          "Non-ASCII characters detected — downstream tools may rewrite them.",
          type = "warning", duration = 6
        )
      }

      design_df <- as.data.frame(
        do.call(rbind, rows),
        stringsAsFactors = FALSE
      )
      colnames(design_df) <- sample_names
      rownames(design_df) <- row_names

      si <- build_sample_info_from_df(
        design_df, sample_names, input$max_group_name_length
      )

      if (is.null(si)) {
        shiny::showNotification(
          "Design not applied — ensure at least one factor has two or more distinct group labels.",
          type = "error", duration = 6
        )
        return()
      }

      gui_design(si)
      shiny::removeModal()
      shiny::showNotification(
        "Experimental design applied.", type = "message", duration = 4
      )
    })

    # Download current grid state as a CSV template ----
    output$download_design_template <- downloadHandler(
      filename = function() "experimental_design.csv",
      content = function(file) {
        req(!is.null(loaded_data()$data))
        sample_names <- colnames(loaded_data()$data)
        n <- length(sample_names)
        n_factors <- 5L
        initial_groups <- detect_groups(sample_names)

        rows <- list()
        row_names <- character(0)
        for (i in seq_len(n_factors)) {
          factor_name_raw <- input[[paste0("factor_name_", i)]]
          name_provided <- !is.null(factor_name_raw) &&
            nchar(trimws(factor_name_raw)) > 0L

          cell_values <- vapply(seq_len(n), function(j) {
            v <- input[[paste0("cell_", i, "_", j)]]
            if (is.null(v)) {
              if (i == 1L) initial_groups[[j]] else ""
            } else {
              v
            }
          }, character(1L))
          cells_provided <- any(nchar(trimws(cell_values)) > 0L)

          # Skip rows where the user has provided neither a name nor any cell;
          # otherwise the download would carry "factor2", "factor3"... empties.
          if (!name_provided && !cells_provided) next

          factor_name <- if (name_provided) {
            trimws(factor_name_raw)
          } else if (i == 1L) {
            "group"
          } else {
            paste0("factor", i)
          }
          rows[[length(rows) + 1L]] <- cell_values
          row_names <- c(row_names, factor_name)
        }

        if (length(rows) == 0L) {
          # Empty grid: write a header-only template so the user gets a usable file.
          df <- as.data.frame(
            matrix("", nrow = 0L, ncol = n),
            stringsAsFactors = FALSE
          )
          colnames(df) <- sample_names
          write.csv(df, file)
          return()
        }

        row_names <- make.unique(row_names)
        df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
        colnames(df) <- sample_names
        rownames(df) <- row_names
        write.csv(df, file)
      }
    )

    # Disables expression_file input to prevent multiple uploads; reset GUI design
    observeEvent(input$expression_file, {
      shinyjs::disable("expression_file")
      gui_design(NULL)
    })

    show_gene_ids_modal <- function() {
      # Build species_choice on-demand (not loaded at startup for performance)
      species_choice <- setNames(as.list(idep_data$org_info$id), idep_data$org_info$name2)
      species_choice <- append(
        setNames("NEW", "**NEW SPECIES**"),
        species_choice
      )

      shiny::showModal(
        shiny::modalDialog(
          title = "Example Gene IDs",
          tags$style(
            HTML(
              "#DataTables_Table_0_wrapper #DataTables_Table_0_filter label{
                width: 400px;
                float: left;
              }"
            )
          ),
          selectizeInput(
            inputId = ns("gene_id_examples"),
            label = "Select or search for species",
            choices = c("--Select species--", names(species_choice))
          ),
          DT::dataTableOutput(ns("showGeneIDs4Species")),
          size = "l",
          easyClose = FALSE
        )
      )
    }

    observeEvent(input$gene_ids_link, {
      show_gene_ids_modal()
    })

    geneIDs <- reactiveVal(NULL)

    observeEvent(input$gene_id_examples, {
      req(input$gene_id_examples != "--Select species--")
      ix <- which(idep_data$org_info$name2 == input$gene_id_examples)
      dbase <- idep_data$org_info$file[ix]
      geneIDs(
        showGeneIDs(
          species = input$gene_id_examples,
          db = dbase,
          nGenes = 10
        )
      )
    })

    # Render Gene ID example table in gene example Modal
    output$showGeneIDs4Species <- DT::renderDataTable({
      req(!is.null(geneIDs()))
      req(input$gene_id_examples != "--Select species--")
      removeNotification("ExampleIDDataQuery")
      DT::datatable(
        geneIDs(),
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        )
      )
    })

    # Disables experiment_file input to prevent multiple uploads;
    # clears any GUI-built design so the uploaded file is the source of truth.
    observeEvent(input$experiment_file, {
      shinyjs::disable("experiment_file")
      gui_design(NULL)
    })

    # Notification for data type selection (species is now always selected)
    observe({
      req(tab() == "Data")
      req(input$data_file_format == 0)

      showNotification(
        "Choose a species first (default is human). Then select a data type.",
        duration = 30,
        type = "error",
        id = "select_first"
      )
    })

    observe({
      req(input$data_file_format != 0 || tab() != "Data")

      removeNotification("select_first")
    })

    # Show messages when on the Network tab or button is clicked ----
    observe({
      req(is.null(loaded_data()$data) && (
        tab() != "Data" && tab() != "Doc"
      ))

      showNotification(
        ui = paste("Load a demo file or your own data first."),
        id = "load_data_first",
        duration = NULL,
        type = "error"
      )
    })

    # Remove messages if the tab changes --------
    observe({
      req(!is.null(loaded_data()$data) ||
        tab() == "Data" || tab() == "Doc")
      removeNotification("load_data_first")
    })

    observe({
      req(!is.null(loaded_data()$data) && any(apply(loaded_data()$data, 2, function(col) {
        values <- col[!is.na(col)]
        length(values) > 0 && all(values == 0)
      })))

      showNotification(
        ui = paste("A sample has all values as zero. It is recommended to remove that sample."),
        id = "sample_remove_error",
        duration = NULL,
        type = "error"
      )
    })

    # Message for the status of the app ---------
    output$file_format <- renderUI({
      shinyjs::hideElement(id = "load_message")
    })

    # Change demo data based on selected format ----
    # returns a vector with file names  c(data, design)
    demo_data_file <- reactive({
      choice <- selected_demo()
      req(choice)

      files <- idep_data$demo_file_info
      ix <- which(files$ID == choice)
      req(length(ix) > 0)
      return(c(
        files$expression[ix],
        files$design[ix]
      ))
    })

    selected_demo_info <- reactive({
      choice <- selected_demo()
      req(choice)

      files <- idep_data$demo_file_info
      entry <- files[files$ID == choice, , drop = FALSE]
      req(nrow(entry) == 1)

      safe_entry_value <- function(df, column) {
        if (!column %in% names(df)) {
          return(NULL)
        }

        value <- df[[column]][[1]]

        if (length(value) == 0) {
          return(NULL)
        }

        value
      }

      list(
        expression = entry$expression[[1]],
        design = entry$design[[1]],
        name = entry$name[[1]],
        memo = safe_entry_value(entry, "memo"),
        ncbi = safe_entry_value(entry, "ncbi")
      )
    })

    output$demo_memo <- renderUI({
      info <- selected_demo_info()
      memo <- info$memo
      ncbi_id <- info$ncbi

      if (is.null(memo) || is.na(memo) || !nzchar(memo)) {
        return(NULL)
      }

      memo_link <- NULL
      if (!is.null(ncbi_id) && !is.na(ncbi_id)) {
        ncbi_clean <- trimws(ncbi_id)
        if (nzchar(ncbi_clean) && nchar(ncbi_clean) >= 4 && startsWith(tolower(ncbi_clean), "gse")) {
          memo_link <- tags$a(
            href = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", ncbi_clean),
            target = "_blank",
            rel = "noopener noreferrer",
            "NCBI"
          )
        }
      }

      if (is.null(memo_link)) {
        return(tags$span(memo))
      }

      tags$span(memo, " ", memo_link)
    })

    output$demo_design_download_ui <- renderUI({
      info <- selected_demo_info()

      design_path <- info$design
      has_design <- !is.null(design_path) && length(design_path) > 0 &&
        !is.na(design_path) && nzchar(design_path)

      if (has_design) {
        tagList(
          downloadButton(
            outputId = ns("download_demo_design"),
            label = NULL,
            class = "btn-default"
          ),
          tippy::tippy_this(
            ns("download_demo_design"),
            "Download the experimental design file for the selected demo dataset.",
            theme = "light"
          )
        )
      } else {
        NULL
      }
    })

    output$download_demo_expression <- downloadHandler(
      filename = function() {
        info <- selected_demo_info()
        paste0("iDEP demo data ", basename(info$expression))
      },
      content = function(file) {
        info <- selected_demo_info()
        req(!is.null(info$expression), nzchar(info$expression), file.exists(info$expression))

        success <- file.copy(info$expression, file, overwrite = TRUE)
        req(success)
      }
    )

    output$download_demo_design <- downloadHandler(
      filename = function() {
        info <- selected_demo_info()
        design_path <- info$design
        req(!is.null(design_path), length(design_path) > 0, !is.na(design_path), nzchar(design_path))

        paste0("iDEP demo data ", basename(design_path))
      },
      content = function(file) {
        info <- selected_demo_info()
        design_path <- info$design
        req(!is.null(design_path), length(design_path) > 0, !is.na(design_path), nzchar(design_path), file.exists(design_path))

        success <- file.copy(design_path, file, overwrite = TRUE)
        req(success)
      }
    )

    # Reactive element to load the data from the user or demo data ---------
    loaded_data <- reactive({
      base <- input_data(
        expression_file = input$expression_file,
        experiment_file = input$experiment_file,
        go_button = go_button_count(),
        demo_data_file = demo_data_file()[1],
        demo_metadata_file = demo_data_file()[2],
        max_group_name_length = input$max_group_name_length
      )
      # Inject GUI-built design when no experiment file is uploaded.
      # Skip if input_data returned NULL — assigning into a NULL base would
      # silently create list(sample_info = ...) and drop $data, $raw_counts, etc.
      if (!is.null(base) &&
          !is.null(gui_design()) &&
          is.null(input$experiment_file)) {
        base$sample_info <- gui_design()
      }
      base
    })

    # Sample information table -----------
    output$sample_info_table <- renderTable(
      {
        req(!is.null(loaded_data()$sample_info))

        isolate({
          tem <- t(loaded_data()$sample_info)
          tem <- cbind(rownames(tem), tem)
          colnames(tem)[1] <- "Study_design"
          as.data.frame(tem, stringsAsFactors = FALSE)
        })
      },
      striped = TRUE,
      bordered = TRUE,
      spacing = "s",
      align = "l",
      rownames = FALSE
    )

    # First 20 rows of dataset table -----------
    output$sample_20 <- renderTable(
      {
        req(!is.null(conversion_info()$converted_data))

        data_preview <- if (nrow(loaded_data()$data) > 20) {
          loaded_data()$data[1:20, , drop = FALSE]
        } else {
          loaded_data()$data
        }

        as.data.frame(data_preview, stringsAsFactors = FALSE)
      },
      striped = TRUE,
      bordered = TRUE,
      spacing = "s",
      align = "l",
      rownames = TRUE
    )

    observeEvent(input$expression_file, {
      # test data for correct format
      if (
        min(loaded_data()$data, na.rm = TRUE) < 0 && input$data_file_format == 1
      ) {
        showModal(modalDialog(
          title = "Data type does not match data format",
          tags$p("Negative values were detected in this dataset. This is not
                 correct for the selected data type (Read Counts).
                 Please double check data type or the data."),
          tags$br(),
          size = "m",
          easyClose = TRUE
          # footer = actionButton(ns("reset_app"), "Start over")
        ))
      }

      # Check for ratio data in fold-change uploads
      if (input$data_file_format %in% c(3, 4)) {
        data <- loaded_data()$data

        # Determine which columns are fold-changes
        n2 <- ncol(data) %/% 2
        has_p_vals <- input$data_file_format == 3
        fc_cols <- if (has_p_vals) {
          2 * (1:n2) - 1 # Fold-change columns (odd columns)
        } else {
          seq_len(ncol(data)) # All columns are fold-changes
        }

        # Check each fold-change column for ratio characteristics
        ratio_detected <- FALSE
        ratio_columns <- c()

        for (i in fc_cols) {
          col_data <- data[, i]
          col_data <- col_data[!is.na(col_data)]

          if (length(col_data) >= 10) {
            has_no_negatives <- sum(col_data < 0) / length(col_data) < 0.01

            if (has_no_negatives) {
              skewness <- e1071::skewness(col_data)

              if (skewness > 1) {
                ratio_detected <- TRUE
                ratio_columns <- c(ratio_columns, colnames(data)[i])
              }
            }
          }
        }

        # Show warning if ratio data detected
        if (ratio_detected) {
          showModal(modalDialog(
            title = "Ratio Data Detected",
            tags$p(
              "The following fold-change column(s) appear to contain ratio data ",
              "instead of log2 fold-changes:"
            ),
            tags$ul(
              lapply(ratio_columns, function(col) tags$li(tags$strong(col)))
            ),
            tags$p(
              style = "margin-top: 15px;",
              "Ratio data characteristics detected:"
            ),
            tags$ul(
              tags$li("No negative values (or < 1% negative)"),
              tags$li("High right skew (upregulated genes > 1, downregulated 0-1)")
            ),
            tags$p(
              style = "margin-top: 15px; color: #d9534f; font-weight: bold;",
              "iDEP will automatically apply log2 transformation during preprocessing."
            ),
            tags$p(
              "This ensures proper analysis and visualization. If your data is already ",
              "in log2 scale, please check your input file."
            ),
            size = "m",
            easyClose = TRUE,
            footer = tagList(
              modalButton("OK")
            )
          ))
        }
      }
    })

    # Get converted IDs ----------
    conversion_info <- reactive({
      req(!is.null(loaded_data()$data))

      req(!(min(loaded_data()$data, na.rm = TRUE) < 0 && input$data_file_format == 1))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Loading Data",
        color = "#000000"
      )

      converted <- convert_id(
        rownames(loaded_data()$data),
        idep_data = idep_data,
        select_org = select_org(),
        max_sample_ids = 200
      )

      all_gene_info <- gene_info(
        converted = converted,
        select_org = select_org(),
        idep_data = idep_data
      )

      converted_data <- convert_data(
        converted = converted,
        no_id_conversion = input$no_id_conversion,
        data = loaded_data()$data,
        multiple_map = input$multiple_map
      )

      all_gene_names <- get_all_gene_names(
        mapped_ids = converted_data$mapped_ids,
        all_gene_info = all_gene_info
      )

      gmt_choices <- gmt_category(
        converted = converted,
        converted_data = converted_data$data,
        select_org = select_org(),
        gmt_file = input$gmt_file,
        idep_data = idep_data
      )

      shinybusy::remove_modal_spinner()

      return(list(
        converted = converted,
        all_gene_info = all_gene_info,
        converted_data = converted_data$data,
        all_gene_names = all_gene_names,
        gmt_choices = gmt_choices
      ))
    })

    # download database for selected species
    observeEvent(select_org(), {
      req(!is.null(select_org()))
      if (identical(select_org(), "NEW")) {
        return()
      }

      ix <- which(idep_data$org_info$id == select_org())
      if (length(ix) == 0) {
        return()
      }

      db_file <- idep_data$org_info[ix, "file"]
      if (length(db_file) == 0) {
        return()
      }

      db_file <- db_file[[1]]
      if (is.na(db_file) || !nzchar(db_file)) {
        return()
      }

      dbname <- file.path(DATAPATH, "db", db_file)
      if (!file.exists(dbname)) {
        download_message <- paste(
          "Downloading database for",
          idep_data$org_info[ix, "name2"],
          "(~5 minutes)"
        )
        shinybusy::show_modal_spinner(
          spin = "orbit",
          text = download_message,
          color = "#000000"
        )
        on.exit(shinybusy::remove_modal_spinner(), add = TRUE)

        # download org_info and demo files to current folder
        options(timeout = 3000)
        download.file(
          url = paste0(db_url, db_ver, "/db/", db_file, ".gz"),
          destfile = paste0(dbname, ".gz"),
          mode = "wb",
          quiet = FALSE
        )
        R.utils::gunzip(
          paste0(dbname, ".gz"),
          remove = TRUE
        ) # untar and unzip the files

        shinybusy::remove_modal_spinner()
        on.exit(NULL)
      }
    })


    species_match_data <- reactive({
      req(select_org())
      if (is.null(input$expression_file) && go_button_count() == 0) {
        return(NULL)
      }
      req(input$data_file_format)
      isolate({
        if (is.null(conversion_info()$converted)) {
          return(as.data.frame("ID not recognized."))
        }
        tem <- conversion_info()$converted$species_match
        if (nrow(tem) > 50) { # show only 50
          tem <- tem[1:50, , drop = FALSE]
        }
        if (is.null(tem)) {
          as.data.frame("ID not recognized.")
        } else {
          data.frame(
            "Species(genes matched)" = tem[, 1],
            check.names = FALSE
          )
        }
      })
    })

    # Track when we've shown the modal to prevent multiple displays
    modal_shown <- reactiveVal(FALSE)

    # Reset modal flag when new data is loaded
    observeEvent(list(input$expression_file, input$go_button),
      {
        modal_shown(FALSE)
      },
      ignoreInit = TRUE
    )

    # Show gene ID error modal when IDs are not recognized
    observe({
      req(tab() == "Data")
      req(select_org() != "NEW")
      req(!modal_shown()) # Only show once per data load

      # Check if we have species match data indicating ID not recognized
      match_data <- species_match_data()
      req(!is.null(match_data))
      req(nrow(match_data) > 0)
      req(match_data[1, 1] == "ID not recognized.")

      # Mark that we've shown the modal
      modal_shown(TRUE)

      showModal(modalDialog(
        title = "Gene IDs not recognized",
        tags$p("Possible causes:"),
        tags$ul(
          tags$li("Wrong species is selected."),
          tags$li(
            "Correct species is selected but we cannot map your ",
            "gene IDs to Ensembl gene IDs or STRING protein IDs."
          ),
          tags$li("Your species is not included in our database.")
        ),
        tags$p(
          "You can still run many analyses except pathway and ",
          "enrichment."
        ),
        size = "m",
        easyClose = FALSE,
        footer = modalButton("OK")
      ))
    })


    output$welcome_ui <- renderUI({
      req(go_button_count() == 0)
      req(input$data_format_help == 0 && input$design_format_help == 0)

      tagList(
        fluidRow(
          column(
            width = 9,
            h4(paste0("iDEP   v", as.character(packageVersion("idepGolem")))),
            h4("Integrated Differential Expression & Pathway analysis"),
          ),
          column(
            width = 3,
            img(
              src = "www/idep_logo.png",
              width = "43",
              height = "50"
            )
          )
        ),
        htmlOutput(ns("file_format")),
        # alternative UI output message for once expression data is loaded
        uiOutput(ns("load_data_alt")),
        includeHTML(app_sys("app/www/messages.html")),
        div(
          style = "position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%;",
          tags$iframe(
            src = "https://www.youtube.com/embed/lqDqrJU-e24?rel=0",
            style = "position: absolute; top: 0; left: 0; width: 100%; height: 100%; border: 0;",
            allow = "accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture",
            allowfullscreen = NA
          )
        ),
        img(
          src = "www/flowchart.png",
          align = "center",
          width = "562",
          height = "383"
        ),
        br(),
        img(
          src = "www/figs.gif",
          align = "center",
          width = "640",
          height = "480"
        )
      )
    })

    # Track which help button was clicked last
    last_help_clicked <- reactiveVal("none")

    observeEvent(input$data_format_help, {
      last_help_clicked("data")
    })

    observeEvent(input$design_format_help, {
      last_help_clicked("design")
    })

    output$format_help_ui <- renderUI({
      req(input$data_format_help != 0 || input$design_format_help != 0)

      # Scroll to section 2 only if design_format_help was clicked last
      should_scroll_to_section2 <- last_help_clicked() == "design"

      tagList(
        includeHTML(app_sys("app/www/format.html")),
        # Scroll to section 2 only if design_format_help was clicked
        if (should_scroll_to_section2) {
          tags$script(HTML("
            setTimeout(function() {
              var anchor = document.querySelector('a[name=\"LFCs\"]');
              if (anchor) {
                anchor.scrollIntoView({ behavior: 'smooth', block: 'start' });
              }
            }, 100);
          "))
        }
      )
    })

    list(
      data_file_format = reactive(input$data_file_format),
      no_fdr = reactive(input$data_file_format == 4),
      select_org = select_org,
      gmt_file = reactive(input$gmt_file),
      sample_info = reactive(loaded_data()$sample_info),
      all_gene_info = reactive(conversion_info()$all_gene_info),
      converted_data = reactive(conversion_info()$converted_data),
      all_gene_names = reactive(conversion_info()$all_gene_names),
      matched_ids = reactive(conversion_info()$converted$ids),
      gmt_choices = reactive(conversion_info()$gmt_choices),
      converted = reactive(conversion_info()$converted),
      no_id_conversion = reactive(input$no_id_conversion),
      plots_color_select = reactive(input$plots_color_select),
      heatmap_color_select = reactive(input$heatmap_color_select),
      select_gene_id = reactive(input$select_gene_id),
      multiple_map = reactive(input$multiple_map),
      plot_grid_lines = reactive(input$plot_grid_lines),
      ggplot2_theme = reactive(input$ggplot2_theme),
      max_groups = reactive(input$max_groups),
      max_group_name_length = reactive(input$max_group_name_length)
    )
  })
}
