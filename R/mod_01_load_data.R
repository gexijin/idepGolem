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
    title = "Load Data",
    # move notifications and progress bar to the center of screen
    tags$head(
      tags$style(
        HTML(".shiny-notification {
              width: 300px;
              position:fixed;
              top: calc(85%);
              left: calc(5%);
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
              "Check this to analyze a new/custom species not in our database.",
              theme = "light-border"
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
            "Upload a custom pathway .GMT file to perform pathway analysis for your new species. This is optional - you can proceed without it.",
            theme = "light-border"
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
            "Read counts data (recommended)" = 1,
            "Normalized Expression data" = 2,
            "Fold-changes & adjusted P-values" = 3
          ),
          selected = 0,
          selectize = FALSE
        ),
        tippy::tippy_this(
          ns("data_file_format"),
          "We recommend uploading Read Counts Data, which can be analyzed using DESeq2.
          Choose Normalized Expression Data if your data is derived from 
          RNA-Seq(FPKM, RPKM, TPM), DNA microarray, proteomics data, etc. 
          Select Fold Change and Adjusted P-value, if you have already conducted D.E.G. analysis.
            ",
          theme = "light-border"
        ),
        
        
        # Conditional panel for fold changes data file ----------
        conditionalPanel(
          condition = "input.data_file_format == 3",
          checkboxInput(
            inputId = ns("no_fdr"),
            label = "Fold-changes only",
            value = FALSE
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
        # tags$style(
        #   HTML("
        #     #load_data-ui {
        #       display: block !important;
        #     }
        #   ")
        # ),
        
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
            "Reveal appearance and ID-conversion settings shared across the app.",
            theme = "light-border"
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
                width = "100%"
              )
            ),
            tippy::tippy_this(
              ns("heatmap_color_select_container"),
              "Choose the color palette used for heatmaps, sample trees, and networks.",
              theme = "light-border"
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
              "When multiple IDs map to the same gene, we can summarize
              data in a certain way (sum, mean, median, max),
              or just keep the rows with the most variation (max SD).
              When uploading transcript level counts, choose \"sum\" to aggregate gene level counts. ",
              theme = "light-border"
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
                width = "100%"
              )
            ),
            tippy::tippy_this(
              ns("plots_color_select_container"),
              "Palette applied to PCA and QC plots so sample groups stand out.",
              theme = "light-border"
            ),
            selectInput(
              inputId = ns("select_gene_id"),
              label = "Gene ID type for plots:",
              choices = c("symbol", "ensembl_ID", "User_ID"),
              selected = "symbol"
            ),
            tippy::tippy_this(
              ns("select_gene_id"),
              "Pick which gene identifier appears on plots and tables.",
              theme = "light-border"
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
              "Changes the styles for all 20 ggplot2 plots.",
              theme = "light-border"
            ),
            checkboxInput(
              inputId = ns("plot_grid_lines"),
              label = "Add grid lines to plots",
              value = FALSE
            ),
            tippy::tippy_this(
              ns("plot_grid_lines"),
              "Overlay light grid lines to help compare values across samples.",
              theme = "light-border"
            ),
            checkboxInput(
              inputId = ns("no_id_conversion"),
              label = "Do not convert gene IDs",
              value = FALSE
            ),
            tippy::tippy_this(
              ns("no_id_conversion"),
              "If selected, uploaded gene IDs will not be converted to ENSEMBL gene IDs,
              which is used as a central id type in pathway databases.",
              theme = "light-border"
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
          h4("Loading R packages, please wait ... ... ...")
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
      updateSelectInput(
        session = session,
        inputId = "data_file_format",
        choices = list("Read counts data (recommended)" = 1,
                       "Normalized Expression data" = 2,
                       "Fold-changes & adjusted P-values" = 3),
        selected = input$data_file_format
        )
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
              "Load demo data",
              theme = "light-border"
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
              "Additional info on accepted data types",
              theme = "light-border"
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
          title = "Demo Data",
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
              "Select a demo file then click the \"Load Demo\" button",
              theme = "light-border"
            ),
            br(),
            actionButton(
              inputId = ns("go_button"),
              label = "Load",
              class = "btn-primary"
            ),
            tippy::tippy_this(
              ns("go_button"),
              "Load the selected demo file",
              theme = "light-border"
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
                a(
                  "video",
                  href = "https://youtu.be/Hs5SamHHG9s",
                  target = "_blank"
                ),
                "tutorial!"
              ),
              tags$li("Select a Species & Data Type"),
              tags$li("Upload data or click ",
                      tags$span("Demo", id = "load-demo"),
                      " to try a sample data set!"
              )
          ),
          tags$script("
            document.getElementById('load-demo').style.color = 'red';
          ")
        )
      } else {
        NULL
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
              "Additional info on accepted data types",
              theme = "light-border"
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
      req(tab() == "Load Data")
      req(input$data_file_format == 0)

        showNotification(
          "Select a data type before uploading data.",
          duration = 30,
          type = "error",
          id = "select_first"
          )
    })

    observe({
      req(input$data_file_format != 0 || tab() != "Load Data")

      removeNotification("select_first")
    })

    # Show messages when on the Network tab or button is clicked ----
    observe({
      req(is.null(loaded_data()$data) && (
        tab() != "Load Data" && tab() != "About"
      ))

      showNotification(
        ui = paste("Please load a demo file or your own data first."),
        id = "load_data_first",
        duration = NULL,
        type = "error"
      )
    })

    # Remove messages if the tab changes --------
    observe({
      req(!is.null(loaded_data()$data) ||
        tab() == "Load Data" || tab() == "About")
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
    # Species match message ----------
    observe({
      req(tab() == "Load Data")

      match_data <- species_match_data()
      req(!is.null(match_data))
      req(nrow(match_data) > 0)
      req(match_data[1, 1] == "ID not recognized.")
      req(input$select_org != "NEW")

      showModal(modalDialog(
        title = "Please double check the selected species",
        tags$p("None of the gene IDs are recognized. Possible causes: 1. Wrong species is selected. 
        2. Correct species is selected but we cannot map your gene IDs to Ensembl gene IDs. 
        3. Your species is not included in our database.  
        You can still run many analyses except pathway and enrichment."),
        size = "s",
        easyClose = TRUE
      ))
    })


    output$welcome_ui <- renderUI({
      req(go_button_count() == 0)
      req(input$data_format_help == 0)

      tagList(
        fluidRow(
          column(
            width = 9,
            h4("iDEP: integrated Differential Expression & Pathway analysis")
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

    output$format_help_ui <- renderUI({
      req(input$data_format_help != 0)

      includeHTML(app_sys("app/www/format.html"))
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
