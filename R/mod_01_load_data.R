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
    # move notifications and progress bar to the center of screen
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
              "
        )
      )
    ),
    sidebarLayout(

      ##################################################################
      #       Load Data sidebar panel ----
      ##################################################################
      sidebarPanel(
        uiOutput(ns("reset_button")),
        # Species Match Drop Down ------------
        selectInput(
          inputId = ns("select_org"),
          label = NULL,
          choices = list(),
          multiple = FALSE,
          selectize = TRUE,
          selected = NULL
        ),

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
                                 }"
        )
        ),

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
        strong("2. Data type"),
        
        selectInput(
          inputId = ns("data_file_format"),
          label = NULL,
          choices = list(
            "..." = 0,
            "Read counts data" = 1,
            "Normalized Expression data" = 2,
            "Fold-changes & adjusted P-vals" = 3
          ),
          selected = 0,
          selectize = FALSE
        ),
        tippy::tippy_this(
          ns("data_file_format"),
          "We recommend raw read counts so iDEP can run DESeq2. Choose normalized expression if you have TPM/FPKM, microarray, or proteomics values. Select Fold-change plus P-values when statistical analysis was done elsewhere.",
          theme = "light"
        ),
        
        
        # Conditional panel for fold changes data file ----------
        conditionalPanel(
          condition = "input.data_file_format == 3",
          checkboxInput(
            inputId = ns("no_fdr"),
            label = "Fold-changes only",
            value = FALSE
          ),
          tippy::tippy_this(
            ns("no_fdr"),
            "Your data contains fold changes but not adjusted p-values.",
            theme = "light"
          ),
          ns = ns
        ),
        # Load expression data options ----------
        # Includes load demo action button, demo data dropdown, and expression
        # file upload box
        conditionalPanel(
          condition = "input.data_file_format != 0",
          uiOutput(ns("load_data_ui")),
          ns = ns
        ),

        # Experiment design file input ----------
        conditionalPanel(
          condition = "input.data_file_format != 0",
          uiOutput(ns("design_file_ui")),
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
              href = "http://bioinformatics.sdstate.edu/reads/"
            ),
          ),
          column(
            width = 4,
            align = "center",
            a(
              "Cite iDEP",
              href = "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6#citeas",
              target = "_blank"
            ),
          ),
          column(
            width = 4,
            align = "right",
            # Action link to reveal Gene ID examples -----------
            actionLink(
              inputId = ns("gene_ids_link"),
              label = "Gene IDs"
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

        DT::dataTableOutput(ns("sample_info_table")),

        # Display first 20 rows of the data ----------
        DT::dataTableOutput(ns("sample_20")),
        div(
          id = ns("load_message"),
          h3("From data to discoveries", style = "color: #d9534f; font-weight: 700;"),
          br(),
          h4("Visualize, analyze, & unveil pathways — in minutes!", style = "color: green;"),
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

    selected_demo <- reactiveVal(NULL)
    demo_preview_content <- reactiveVal(NULL)

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
      )
    )

    get_data_type_details <- function(type) {
      switch(
        as.character(type),
        `1` = list(
          title = "Read Counts",
          body = tagList(
            p("Raw gene-by-sample count matrix. Values indicate the number of sequencing reads assigned to each gene. Counts are typically integers, but estimated counts (e.g., from kallisto or Salmon) may be non-integer and should still be treated as count data."),
            tags$ul(
              tags$li("First column: gene IDs, such as Ensembl, Entrez, symbols, etc."),
              tags$li("Column headers: sample names; avoid spaces and '-' characters.")
            )
          )
        ),
        `2` = list(
          title = "Normalized Expression Matrix",
          body = tagList(
            p("Provide a gene-by-sample matrix with normalized values such as TPM/FPKM, microarray intensities, proteomics, etc."),
            tags$ul(
              tags$li("First column: gene IDs such as Ensembl, Entrez, symbols, etc."),
              tags$li("Column headers: sample names; avoid spaces and '-' characters.")
            )
          )
        ),
        `3` = list(
          title = "Fold Change & Adjusted P-values",
          body = tagList(
            p("Upload summary statistics for one or more contrasts."),
            tags$ul(
              tags$li("First column: gene IDs (Ensembl, symbols, ...)"),
              tags$li("Log2 fold-change and its matching adjusted P-value/FDR"),
              tags$li("Pair each contrast with an adjusted P-value/FDR column.")
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
      if(!is.null(idep_data$org_info) && nrow(idep_data$org_info) > 0) {
        # Get first species from org_info (ordered by 'top' field)
        first_species_id <- idep_data$org_info$id[1]

        # Create all species choices including NEW for custom species
        all_choices <- setNames(as.list(idep_data$org_info$id), idep_data$org_info$name2)
        all_choices <- c(all_choices, setNames("NEW", "**NEW SPECIES**"))

        updateSelectInput(
          session = session,
          inputId = "select_org",
          choices = all_choices,
          selected = first_species_id
        )

        # Set initial selected species name
        selected_species_name(idep_data$org_info$name2[1])
      }
    })

    # Hide the species dropdown
    shinyjs::hideElement(id = "select_org")

    # Pop-up modal for gene assembl information ----
    observeEvent(input$genome_assembl_button, {
      # Create species count summary by source
      org_info_df <- idep_data$org_info
      source_counts <- table(org_info_df$group)

      # Categorize sources according to requirements
      ensembl_sources <- source_counts[!grepl("^newSpecies_|^STRING", names(source_counts))]
      custom_count <- sum(source_counts[grepl("^newSpecies_", names(source_counts))])
      string_count <- sum(source_counts[grepl("^STRING", names(source_counts))])

      # Create Ensembl breakdown
      ensembl_parts <- character(0)
      if("ENSEMBL" %in% names(ensembl_sources)) {
        ensembl_parts <- c(ensembl_parts, paste0("main (", ensembl_sources["ENSEMBL"], ")"))
      }
      if("plants" %in% names(ensembl_sources)) {
        ensembl_parts <- c(ensembl_parts, paste0("plants (", ensembl_sources["plants"], ")"))
      }
      if("bacteria" %in% names(ensembl_sources)) {
        ensembl_parts <- c(ensembl_parts, paste0("bacteria (", ensembl_sources["bacteria"], ")"))
      }
      if("fungi" %in% names(ensembl_sources)) {
        ensembl_parts <- c(ensembl_parts, paste0("fungi (", ensembl_sources["fungi"], ")"))
      }
      if("protists" %in% names(ensembl_sources)) {
        ensembl_parts <- c(ensembl_parts, paste0("protists (", ensembl_sources["protists"], ")"))
      }
      if("metazoa" %in% names(ensembl_sources)) {
        ensembl_parts <- c(ensembl_parts, paste0("metazoa (", ensembl_sources["metazoa"], ")"))
      }

      # Create count summary text
      all_parts <- character(0)
      if(length(ensembl_parts) > 0) {
        all_parts <- c(all_parts, paste("Ensembl:", paste(ensembl_parts, collapse = ", ")))
      }
      if(string_count > 0) {
        all_parts <- c(all_parts, paste0("STRING-db (", string_count, ")"))
      }
      if(custom_count > 0) {
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
            df <- idep_data$org_info[,
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
        updateSelectInput(
          session = session,
          inputId = "select_org",
          selected = "NEW"
        )
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

        updateSelectInput(
          session = session,
          inputId = "select_org",
          selected = first_species_id
        )
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
      if(length(lines) > 3){
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

      # Update species selection
      updateSelectInput(
        session = session,
        inputId = "select_org",
        selected = selected_id
      )

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
      if(input$select_org == first_species_id && input$select_org != "NEW") {
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
      req(input$select_org != "NEW")
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
        style = paste0("background-color: ", bg_color, "; border: 1px solid ",
                       border_color, "; border-radius: 4px; padding: 10px; margin: 10px 0;
                       color: ", text_color, ";"),
        h5(
          icon(icon_name),
          paste0(" ", message_type, " ", converted_count, " out of ",
                 original_count, " genes (", conversion_rate,
                 "%) converted to Ensembl/STRING IDs.")
        )
      )
    })

    observeEvent(input$data_file_format, {
      req(input$data_file_format != 0)

      # Update dropdown text
      updateSelectInput(
        session = session,
        inputId = "data_file_format",
        choices = list("Read counts data" = 1,
                       "Normalized Expression data" = 2,
                       "Fold-change & adjusted P-val" = 3),
        selected = input$data_file_format
      )

      # Show notification only when on Data tab
      req(tab() == "Data")

      details <- get_data_type_details(input$data_file_format)
      preview <- demo_preview_tables[[as.character(input$data_file_format)]]

      demo_preview_content(preview)

      # Build notification UI with details and preview
      notification_ui <- if (!is.null(preview)) {
        # Convert preview to character to preserve exact display
        df <- preview
        df[] <- lapply(df, as.character)

        # Create HTML table manually
        table_html <- tags$table(
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
    }, ignoreNULL = TRUE)

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

      if (input$data_file_format > 0) {
        files <- idep_data$demo_file_info
        files <- files[files$type == input$data_file_format, ]
        setNames(as.list(files$ID), files$name)
      } else {
        NULL
      }
    })

    observeEvent(demo_choices(), {
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
    }, ignoreNULL = FALSE)

    observeEvent(input$select_demo, {
      req(!is.null(input$select_demo))
      selected_demo(input$select_demo)
    }, ignoreNULL = TRUE)

    # UI elements for load demo action button, demo data drop down, and -----
    # expression file upload
    output$load_data_ui <- renderUI({
      req(go_button_count() == 0)
      req(input$data_file_format)

      choices <- demo_choices()
      tagList(
        strong("3. Expression matrix (CSV, text, or xlsx)"),
        fluidRow(
          column(
            width = 6,
            # Expression data file input
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
              placeholder = "",
              width = "100%"
            ),
            tippy::tippy_this(
              ns("expression_file"),
              "Upload your expression matrix in CSV, TSV, or Excel format.",
              theme = "light"
            )
          ),

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
          ),
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
              "Watch a ",
              a("video", href = "https://youtu.be/Hs5SamHHG9s", target = "_blank"),
              " tutorial!"
            ),
            tags$li(
              "Try it with demo data. After selecting a data type, just click ",
              tags$span("Load Demo.", id = "load-demo", style = "color: red;")
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
              align = "right",
              actionButton(
                inputId = ns("reset_app_new_data"),
                label = strong("Reset"),
                align = "right"
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

        strong("4. Optional: Exp. Design (CSV or text)"),
        fluidRow(
          column(
            width = 10,
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
            ),
            tippy::tippy_this(
              ns("experiment_file"),
              "Upload a sample information table to define experimental design.",
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

    # Disables expression_file input to prevent multiple uploads
    observeEvent(input$expression_file, {
      shinyjs::disable("expression_file")
    })

    show_gene_ids_modal <- function() {
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
            choices = c("--Select species--", names(idep_data$species_choice))
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

    # Disables experiment_file input to prevent multiple uploads
    observeEvent(input$experiment_file, {
      shinyjs::disable("experiment_file")
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
      input_data(
        expression_file = input$expression_file,
        experiment_file = input$experiment_file,
        go_button = go_button_count(),
        demo_data_file = demo_data_file()[1],
        demo_metadata_file = demo_data_file()[2]
      )
    })

    # Sample information table -----------
    output$sample_info_table <- DT::renderDataTable({
      req(!is.null(loaded_data()$sample_info))

      DT::datatable(
        isolate({
          tem <- t(loaded_data()$sample_info)
          tem <- cbind(rownames(tem), tem)
          colnames(tem)[1] <- "Study_design"
          tem
        }),
        options = list(
          pageLength = 10,
          scrollX = "400px",
          dom = "t",
          ordering = F
        ),
        rownames = FALSE
      )
    })

    # First 20 rows of dataset table -----------
    output$sample_20 <- DT::renderDataTable({
      req(!is.null(conversion_info()$converted_data))

      DT::datatable(
        #conversion_info()$converted_data[1:20, ],
        if(nrow(loaded_data()$data) > 20){
          loaded_data()$data[1:20, ]
        } else (loaded_data()$data),
        options = list(
          pageLength = 10,
          scrollX = "400px",
          dom = "t"
        ),
        rownames = TRUE
      )
    })

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
          #footer = actionButton(ns("reset_app"), "Start over")
        ))
      }

      # Check for ratio data in fold-change uploads
      if (input$data_file_format == 3) {
        data <- loaded_data()$data

        # Determine which columns are fold-changes
        n2 <- ncol(data) %/% 2
        fc_cols <- if (!input$no_fdr) {
          2 * (1:n2) - 1  # Fold-change columns (odd columns)
        } else {
          1:ncol(data)    # All columns are fold-changes
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
        select_org = input$select_org,
        max_sample_ids = 200
      )

      all_gene_info <- gene_info(
        converted = converted,
        select_org = input$select_org,
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
        select_org = input$select_org,
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
    observeEvent(input$select_org, {
      req(!is.null(input$select_org))
      if (identical(input$select_org, "NEW")) {
        return()
      }

      ix <- which(idep_data$org_info$id == input$select_org)
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
        withProgress(
          message = paste(
            "Download database for",
            idep_data$org_info[ix, "name2"],
           "(~5 minutes)"
           ), {
          incProgress(0.2)
          # download org_info and demo files to current folder
          options(timeout = 3000)
          download.file(
            url = paste0(db_url, db_ver, "/db/", db_file, ".gz"),
            destfile = paste0(dbname, ".gz"),
            mode = "wb",
            quiet = FALSE
          )
          incProgress(0.7)
          R.utils::gunzip(
            paste0(dbname, ".gz"), 
            remove = TRUE
          ) # untar and unzip the files
        })

      }
    })


    species_match_data <- reactive({
      req(input$select_org)
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
    observeEvent(list(input$expression_file, input$go_button), {
      modal_shown(FALSE)
    }, ignoreInit = TRUE)

    # Show gene ID error modal when IDs are not recognized
    observe({
      req(tab() == "Data")
      req(input$select_org != "NEW")
      req(!modal_shown())  # Only show once per data load

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
          tags$li("Correct species is selected but we cannot map your ",
                  "gene IDs to Ensembl gene IDs or STRING protein IDs."),
          tags$li("Your species is not included in our database.")
        ),
        tags$p("You can still run many analyses except pathway and ",
               "enrichment."),
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
            h4("iDEP: integrated Differential Expression & Pathway analysis (v2.20)"),
            h5("The power of 100s of R packages and annotation databases, at your fingertips!")
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
        br(),
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
      no_fdr = reactive(input$no_fdr),
      select_org = reactive(input$select_org),
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
      ggplot2_theme = reactive(input$ggplot2_theme)
    )
  })
}
