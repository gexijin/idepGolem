#' 06_pathway UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_06_pathway_ui <- function(id) {
  ns <- shiny::NS(id)
  tabPanel(
    title = "Pathway",
    sidebarLayout(
      sidebarPanel(
        htmlOutput(
          outputId = ns("list_comparisons_pathway")
        ),
        tags$style(
          type = "text/css",
          "#pathway-list_comparisons_pathway { width:100%;   margin-top:-12px}"
        ),
        selectInput(
          inputId = ns("pathway_method"),
          label = "Select method:",
          choices = list(
            "GAGE" = 1,
            "GSEA (preranked fgsea)" = 3,
            "PGSEA" = 2,
            "PGSEA w/ all samples" = 4,
            "ReactomePA" = 5
          ),
          selected = 1
        ),
        tags$style(
          type = "text/css",
          "#pathway-pathway_method { width:100%;   margin-top:-12px}"
        ),
        htmlOutput(
          outputId = ns("select_go_selector")
        ),
        tags$style(
          type = "text/css",
          "#pathway-select_go { width:100%;   margin-top:-12px}"
        ),
        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = ns("min_set_size"),
              label = "Geneset size: Min.",
              min = 5,
              max = 30,
              value = 15,
              step = 1
            )
          ),
          column(
            width = 6,
            numericInput(
              inputId = ns("max_set_size"),
              label = "Max.",
              min = 1000,
              max = 2000,
              value = 2000,
              step = 100
            )
          )
        ),
        tags$style(
          type = "text/css",
          "#pathway-min_set_size { width:100%;   margin-top:-12px}"
        ),
        tags$style(
          type = "text/css",
          "#pathway-max_set_size { width:100%;   margin-top:-12px}"
        ),
        numericInput(
          inputId = ns("pathway_p_val_cutoff"),
          label = "Pathway signifiance cutoff (FDR)",
          value = 0.2,
          min = 1e-20,
          max = 1,
          step = .05
        ),
        tags$style(
          type = "text/css",
          "#pathway-pathway_p_val_cutoff { width:100%;   margin-top:-12px}"
        ),
        numericInput(
          inputId = ns("n_pathway_show"),
          label = "Number of top pathways to show",
          value = 30,
          min = 5,
          max = 100,
          step = 5
        ),
        tags$style(
          type = "text/css",
          "#pathway-n_pathway_show { width:100%;   margin-top:-12px}"
        ),
        checkboxInput(
          inputId = ns("absolute_fold"),
          label = "Use absolute values of fold changes for GSEA and GAGE",
          value = FALSE
        ),
        numericInput(
          inputId = ns("gene_p_val_cutoff"),
          label = "Remove genes with big FDR before pathway analysis:",
          value = 1,
          min = 1e-20,
          max = 1,
          step = .05
        ),
        tags$style(
          type = "text/css",
          "#pathway-gene_p_val_cutoff { width:100%;   margin-top:-12px}"
        ),
        h6("* Warning! The many combinations can lead to false positives in pathway analyses."),
        # Download report button
        downloadButton(
          outputId = ns("report"),
          label = "Generate Report"
        ),
        tippy::tippy_this(
          ns("report"),
          "Generate HTML report of pathway tab",
          theme = "light-border"
        ),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/pathways/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          id = ns("pathway_tabs"),
          tabPanel(
            title = "Significant pathways",
            br(),
            conditionalPanel(
              condition = "input.pathway_method == 1",
              tableOutput(ns("gage_pathway_table")),
              ns = ns
            ),
            conditionalPanel(
              condition = "input.pathway_method == 2",
              h5("Red and blue indicates relatively activated
               and suppressed pathways, respectively.
              GS just indicates a color scale."),
              plotOutput(
                outputId = ns("pgsea_plot"),
                inline = TRUE
              ),
              ns = ns
            ),
            conditionalPanel(
              condition = "input.pathway_method == 3",
              tableOutput(outputId = ns("fgsea_pathway")),
              ns = ns
            ),
            conditionalPanel(
              condition = "input.pathway_method == 4",
              h5("Red and blue indicates activated and suppressed pathways, respectively."),
              plotOutput(
                outputId = ns("pgsea_plot_all_samples"),
                inline = TRUE
              ),
              ns = ns
            ),
            conditionalPanel(
              condition = "input.pathway_method == 5",
              DT::dataTableOutput(outputId = ns("reactome_pa_pathway")),
              ns = ns
            )
          ),
          tabPanel(
            title = "Tree",
            plotOutput(
              outputId = ns("enrichment_tree"),
              width = "100%"
            ),
            br(),
            p("Adjusting the width of the browser
            window can render figure differently and
            resolve the \"Figure margin too wide\" error. ")
          ),
          tabPanel(
            title = "Network",
            p("Connected gene sets share more genes. Color of node correspond to adjuested Pvalues."),
            fluidRow(
              column(
                width = 2,
                actionButton(
                  inputId = ns("layout_vis_deg"),
                  label = "Change layout"
                )
              ),
              column(
                width = 1,
                h5("Cutoff:"),
                align = "right"
              ),
              column(
                width = 2,
                numericInput(
                  inputId = ns("edge_cutoff_deg"),
                  label = NULL,
                  value = 0.30,
                  min = 0,
                  max = 1,
                  step = .1
                ),
                align = "left"
              ),
              column(
                width = 2,
                checkboxInput(
                  inputId = ns("wrap_text_network_deg"),
                  label = "Wrap text",
                  value = TRUE
                )
              ),
              column(
                width = 3,
                selectInput(
                  inputId = ns("up_down_reg_deg"),
                  label = NULL,
                  choices = c(
                    "Both Up & Down" = "All Groups",
                    "Up regulated" = "Up",
                    "Down regulated" = "Down"
                  )
                )
              )
            ),
            h6(
              "Two pathways (nodes) are connected if they share 30% (default, adjustable) or more genes.
              Green and red represents down- and up-regulated pathways. You can move the nodes by
              dragging them, zoom in and out by scrolling, and shift the entire network by click on an
              empty point and drag. Darker nodes are more significantly enriched gene sets. Bigger nodes
              represent larger gene sets. Thicker edges represent more overlapped genes."
            ),
            visNetwork::visNetworkOutput(
              outputId = ns("vis_network_path"),
              height = "800px",
              width = "100%"
            )
          ),
          tabPanel(
            title = "Heatmap",
            htmlOutput(outputId = ns("list_sig_pathways")),
            mod_12_heatmap_ui(ns("12_heatmap_1"))
          ),
          tabPanel(
            title = "KEGG",
            conditionalPanel(
              condition = "(input.pathway_method == 1 | input.pathway_method == 2 |
                            input.pathway_method == 3 | input.pathway_method == 4) &
                            input.select_go != 'KEGG'",
              h5("Please select KEGG database, if available, from left and perform pathway analysis first."),
              ns = ns
            ),
            conditionalPanel(
              condition = "(input.pathway_method == 1 | input.pathway_method == 2 |
                            input.pathway_method == 3 | input.pathway_method == 4) &
                            input.select_go == 'KEGG'",
              fluidRow(
                column(
                  width = 6,
                  htmlOutput(outputId = ns("list_sig_pathways_kegg"))
                ),
                column(
                  width = 3,
                  checkboxInput(
                    inputId = ns("kegg_sig_only"),
                    label = "All KEGG pathways",
                    value = FALSE
                  )
                ),
                column(
                  width = 3,
                  selectInput(
                    inputId = ns("kegg_color_select"),
                    label = "Colors (low-high)",
                    choices = "green-red",
                    width = "100%"
                  )
                )
              ),
              imageOutput(
                outputId = ns("kegg_image"),
                width = "100%",
                height = "100%"
              ),
              h5("Red and green represent up- and down-regulated genes, respectively."),
              ns = ns
            )
          ),
          tabPanel(
            title = "Info",
            includeHTML("inst/app/www/pathway.html")
          )
        )
      )
    )
  )
}

#' mod_06_pathway Server Functions
#'
#' @noRd
mod_06_pathway_server <- function(id, pre_process, deg, idep_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Interactive heatmap environment
    path_env <- new.env()

    # GMT choices for enrichment ----------
    output$select_go_selector <- renderUI({
      req(!is.null(pre_process$gmt_choices()))
      req(input$pathway_method != 5)
      # if there is KEGG, use KEGG as default
      selected <- "GOBP"
      if ("KEGG" %in% pre_process$gmt_choices()) {
        selected <- "KEGG"
      }
      selectInput(
        inputId = ns("select_go"),
        label = "Select Geneset:",
        choices = pre_process$gmt_choices(),
        selected = selected
      )
    })

    # Dynamic Barplot Tab ----------
    observe({
      if (input$pathway_method == 5) {
        hideTab(inputId = "pathway_tabs", target = "Tree")
        hideTab(inputId = "pathway_tabs", target = "Network")
        hideTab(inputId = "pathway_tabs", target = "Heatmap")
        hideTab(inputId = "pathway_tabs", target = "KEGG")
      } else {
        showTab(inputId = "pathway_tabs", target = "Tree")
        showTab(inputId = "pathway_tabs", target = "Network")
        showTab(inputId = "pathway_tabs", target = "Heatmap")
        showTab(inputId = "pathway_tabs", target = "KEGG")
      }
    })

    # If data is uploaded, but DEG1 is not run
    observe({
      req(!is.null(pre_process$data()) && is.null(deg$limma()) && (
        tab() == "DEG1" || tab() == "DEG2" ||
          tab() == "Pathway" || tab() == "Genome"
      ))

      showNotification(
        ui = paste("Differentially expressed genes need to
        be identified first. Please select factors and comparisons and
        click Submit on the DEG1 tab."),
        id = "click_submit_DEG1",
        duration = NULL,
        type = "error"
      )
    })

    # Remove messages if the tab changes --------
    observe({
      req(!is.null(deg$limma()) || (
        tab() != "DEG1" && tab() != "DEG2" &&
          tab() != "Pathway" && tab() != "Genome"
      ))
      removeNotification("click_submit_DEG1")
    })

    output$list_comparisons_pathway <- renderUI({
      if (is.null(deg$limma()$comparisons)) {
        selectInput(
          inputId = ns("select_contrast"),
          label = NULL,
          choices = list("All" = "All"),
          selected = "All"
        )
      } else {
        selectInput(
          inputId = ns("select_contrast"),
          label =
            "Select a comparison to examine. \"A-B\" means A vs. B (See heatmap).
            Interaction terms start with \"I:\"",
          choices = deg$limma()$comparisons
        )
      }
    })

    output$list_sig_pathways <- renderUI({
      req(!is.null(input$pathway_method))
      # Default, sometimes these methods returns "No significant pathway found"
      choices <- "All"
      if (input$pathway_method == 1) {
        if (!is.null(gage_pathway_data())) {
          if (dim(gage_pathway_data())[2] > 1) {
            choices <- gage_pathway_data()[, 2]
          }
        }
      } else if (input$pathway_method == 2) {
        if (!is.null(pgsea_plot_data())) {
          if (dim(pgsea_plot_data())[2] > 1) {
            pathways <- as.data.frame(pgsea_plot_data())
            choices <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
          }
        }
      } else if (input$pathway_method == 3) {
        if (!is.null(fgsea_pathway_data())) {
          if (dim(fgsea_pathway_data())[2] > 1) {
            choices <- fgsea_pathway_data()[, 2]
          }
        }
      } else if (input$pathway_method == 4) {
        if (!is.null(pgsea_plot_all_samples_data())) {
          if (dim(pgsea_plot_all_samples_data())[2] > 1) {
            pathways <- as.data.frame(pgsea_plot_all_samples_data())
            choices <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
          }
        }
      } else if (input$pathway_method == 5) {
        if (!is.null(reactome_pa_pathway_data())) {
          if (dim(reactome_pa_pathway_data())[2] > 1) {
            choices <- reactome_pa_pathway_data()[, 2]
          }
        }
      }

      selectInput(
        inputId = ns("sig_pathways"),
        label = "Select a significant pathway:",
        choices = choices
      )
    })

    output$list_sig_pathways_kegg <- renderUI({
      req(!is.null(input$pathway_method))
      # Default, sometimes these methods returns "No significant pathway found"
      choices <- "All"
      if (input$pathway_method == 1) {
        if (!is.null(gage_pathway_data())) {
          if (dim(gage_pathway_data())[2] > 1) {
            choices <- gage_pathway_data()[, 2]
          }
        }
      } else if (input$pathway_method == 2) {
        if (!is.null(pgsea_plot_data())) {
          if (dim(pgsea_plot_data())[2] > 1) {
            pathways <- as.data.frame(pgsea_plot_data())
            choices <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
          }
        }
      } else if (input$pathway_method == 3) {
        if (!is.null(fgsea_pathway_data())) {
          if (dim(fgsea_pathway_data())[2] > 1) {
            choices <- fgsea_pathway_data()[, 2]
          }
        }
      } else if (input$pathway_method == 4) {
        if (!is.null(pgsea_plot_all_samples_data())) {
          if (dim(pgsea_plot_all_samples_data())[2] > 1) {
            pathways <- as.data.frame(pgsea_plot_all_samples_data())
            choices <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
          }
        }
      } else if (input$pathway_method == 5) {
        if (!is.null(reactome_pa_pathway_data())) {
          if (dim(reactome_pa_pathway_data())[2] > 1) {
            choices <- reactome_pa_pathway_data()[, 2]
          }
        }
      }

      # show all pathways if selected
      if (input$kegg_sig_only && !is.null(gene_sets())) {
        choices <- names(gene_sets()$gene_lists)
      }
      selectInput(
        inputId = ns("sig_pathways_kegg"),
        label = "Select a KEGG pathway:",
        choices = choices
      )
    })

    gene_sets <- reactive({
      #      req(tab() == "Pathway")
      req(!is.null(input$select_go))

      withProgress(message = "Reading pathway database", {
        incProgress(0.2)
        if (pre_process$select_org() == "NEW" && !is.null(pre_process$gmt_file())) {
          in_file <- pre_process$gmt_file()
          in_file <- in_file$datapath
          gene_sets <- read_gmt_robust(in_file)
        } else {
          gene_sets <- read_gene_sets(
            converted = pre_process$converted(),
            all_gene_names = pre_process$all_gene_names(),
            go = input$select_go,
            select_org = pre_process$select_org(),
            idep_data = idep_data,
            my_range = c(input$min_set_size, input$max_set_size)
          )
        }
      })

      return(gene_sets)
    })

    gage_pathway_data <- reactive({
      req(input$pathway_method == 1)
      req(!is.null(deg$limma()))
      req(!is.null(gene_sets()))

      withProgress(message = "Running GAGE", {
        incProgress(0.2)
        gage_data(
          select_go = input$select_go,
          select_contrast = input$select_contrast,
          min_set_size = input$min_set_size,
          max_set_size = input$max_set_size,
          limma = deg$limma(),
          gene_p_val_cutoff = input$gene_p_val_cutoff,
          gene_sets = gene_sets()$gene_lists,
          absolute_fold = input$absolute_fold,
          pathway_p_val_cutoff = input$pathway_p_val_cutoff,
          n_pathway_show = input$n_pathway_show
        )
      })
    })

    output$gage_pathway_table <- renderTable(
      {
        req(!is.null(gage_pathway_data()))

        res <- gage_pathway_data()
        if (ncol(res) > 1) {
          # add URL
          ix <- match(res[, 2], gene_sets()$pathway_info$description)
          res[, 2] <- hyperText(
            res[, 2],
            gene_sets()$pathway_info$memo[ix]
          )
          res$Genes <- as.character(res$Genes)
        }
        return(res)
      },
      digits = -1,
      spacing = "s",
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = TRUE,
      sanitize.text.function = function(x) x
    )

    contrast_samples <- reactive({
      req(!is.null(input$select_contrast))
      req(!is.null(pre_process$data()))

      find_contrast_samples(
        select_contrast = input$select_contrast,
        all_sample_names = colnames(pre_process$data()),
        sample_info = pre_process$sample_info(),
        select_factors_model = deg$select_factors_model(),
        select_model_comprions = deg$select_model_comprions(),
        reference_levels = deg$reference_levels(),
        counts_deg_method = deg$counts_deg_method(),
        data_file_format = pre_process$data_file_format()
      )
    })

    output$pgsea_plot <- renderPlot(
      {
        req(input$pathway_method == 2)
        withProgress(message = "Running PGSEA...", {
          incProgress(0.2)
          plot_pgsea(
            my_range = c(input$min_set_size, input$max_set_size),
            processed_data = pre_process$data(),
            contrast_samples = contrast_samples(),
            gene_sets = gene_sets()$gene_lists,
            pathway_p_val_cutoff = input$pathway_p_val_cutoff,
            n_pathway_show = input$n_pathway_show
          )
        })
      },
      height = 800,
      width = 800
    )

    pgsea_plot_data <- reactive({
      req(input$pathway_method == 2)
      req(!is.null(gene_sets()))
      withProgress(message = "Running PGSEA...", {
        incProgress(0.2)
        get_pgsea_plot_data(
          my_range = c(input$min_set_size, input$max_set_size),
          data = pre_process$data(),
          select_contrast = input$select_contrast,
          gene_sets = gene_sets()$gene_lists,
          sample_info = pre_process$sample_info(),
          select_factors_model = deg$select_factors_model(),
          select_model_comprions = deg$select_model_comprions(),
          pathway_p_val_cutoff = input$pathway_p_val_cutoff,
          n_pathway_show = input$n_pathway_show
        )
      })
    })

    fgsea_pathway_data <- reactive({
      req(input$pathway_method == 3)
      req(!is.null(deg$limma()))
      req(!is.null(gene_sets()))
      withProgress(message = "Permutations in fgsea may take several minutes ...", {
        incProgress(0.2)
        fgsea_data(
          select_contrast = input$select_contrast,
          my_range = c(input$min_set_size, input$max_set_size),
          limma = deg$limma(),
          gene_p_val_cutoff = input$gene_p_val_cutoff,
          gene_sets = gene_sets()$gene_lists,
          absolute_fold = input$absolute_fold,
          pathway_p_val_cutoff = input$pathway_p_val_cutoff,
          n_pathway_show = input$n_pathway_show
        )
      })
    })

    output$fgsea_pathway <- renderTable(
      {
        req(input$pathway_method == 3)
        req(!is.null(fgsea_pathway_data()))

        res <- fgsea_pathway_data()
        if (ncol(res) > 1) {
          # add URL
          ix <- match(res[, 2], gene_sets()$pathway_info$description)
          res[, 2] <- hyperText(
            res[, 2],
            gene_sets()$pathway_info$memo[ix]
          )
          res$Genes <- as.character(res$Genes)
        }

        return(res)
      },
      digits = -1,
      spacing = "s",
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = TRUE,
      sanitize.text.function = function(x) x
    )



    reactome_pa_pathway_data <- reactive({
      req(input$pathway_method == 5)
      req(!is.null(deg$limma()))
      withProgress(message = "ReactomePA may take 5 minutes...", {
        incProgress(0.2)
        reactome_data(
          select_contrast = input$select_contrast,
          my_range = c(input$min_set_size, input$max_set_size),
          limma = deg$limma(),
          gene_p_val_cutoff = input$gene_p_val_cutoff,
          converted = pre_process$converted(),
          idep_data = idep_data,
          pathway_p_val_cutoff = input$pathway_p_val_cutoff,
          n_pathway_show = input$n_pathway_show,
          absolute_fold = input$absolute_fold
        )
      })
    })

    output$pgsea_plot_all_samples <- renderPlot(
      {
        req(input$pathway_method == 4)
        withProgress(message = "Running PGSEA on all samples ...", {
          incProgress(0.2)
          pgsea_plot_all(
            go = input$select_go,
            my_range = c(input$min_set_size, input$max_set_size),
            data = pre_process$data(),
            select_contrast = input$select_contrast,
            gene_sets = gene_sets()$gene_lists,
            pathway_p_val_cutoff = input$pathway_p_val_cutoff,
            n_pathway_show = input$n_pathway_show
          )
        })
      },
      height = 800,
      width = 800
    )

    pgsea_plot_all_samples_data <- reactive({
      req(input$pathway_method == 4)
      req(!is.null(gene_sets()))
      withProgress(message = "Running PGSEA on all samples ...", {
        incProgress(0.2)
        get_pgsea_plot_all_samples_data(
          data = pre_process$data(),
          select_contrast = input$select_contrast,
          gene_sets = gene_sets()$gene_lists,
          my_range = c(input$min_set_size, input$max_set_size),
          pathway_p_val_cutoff = input$pathway_p_val_cutoff,
          n_pathway_show = input$n_pathway_show
        )
      })
    })

    output$reactome_pa_pathway <- DT::renderDataTable({
      req(input$pathway_method == 5)
      req(!is.null(reactome_pa_pathway_data()))

      DT::datatable(
        reactome_pa_pathway_data(),
        options = list(
          pageLength = 15,
          scrollX = "500px"
        ),
        rownames = FALSE
      )
    })

    selected_pathway_data <- reactive({
      req(!is.null(input$sig_pathways))
      req(!is.null(gene_sets()$gene_lists))

      pathway_select_data(
        sig_pathways = input$sig_pathways,
        gene_sets = gene_sets()$gene_lists,
        contrast_samples = contrast_samples(),
        data = pre_process$data(),
        select_org = pre_process$select_org(),
        all_gene_names = pre_process$all_gene_names()
      )
    })

    # Kegg Colors --------
    kegg_colors <- list(
      "Green-Red" = c("green", "red"),
      "Blue-Orange" = c("blue", "orange")
    )
    kegg_choices <- c(
      "Green-Red",
      "Blue-Orange"
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "kegg_color_select",
        choices = kegg_choices
      )
    })

    heatmap_module <- mod_12_heatmap_server(
      id = "12_heatmap_1",
      data = reactive({
        selected_pathway_data()
      }),
      bar = function() {
        return(NULL)
      },
      all_gene_names = reactive({
        pre_process$all_gene_names()
      }),
      cluster_rows = TRUE,
      heatmap_color = reactive({
        pre_process$heatmap_color_select()
      }),
      select_gene_id = reactive({
        pre_process$select_gene_id()
      })
    )


    output$kegg_image <- renderImage(
      {
        req(!is.null(input$sig_pathways_kegg))
        withProgress(message = "Downloading KEGG pathway", {
          incProgress(0.2)
          kegg_pathway(
            go = input$select_go,
            gage_pathway_data = pathway_list_data()[, 1:5],
            sig_pathways = input$sig_pathways_kegg,
            select_contrast = input$select_contrast,
            limma = deg$limma(),
            converted = pre_process$converted(),
            idep_data = idep_data,
            select_org = pre_process$select_org(),
            low_color = kegg_colors[[input$kegg_color_select]][1],
            high_color = kegg_colors[[input$kegg_color_select]][2]
          )
        })
      },
      deleteFile = TRUE
    )

    # List of pathways with details
    pathway_list_data <- reactive({
      get_pathway_list_data(
        pathway_method = input$pathway_method,
        gage_pathway_data = gage_pathway_data(),
        fgsea_pathway_data = fgsea_pathway_data(),
        pgsea_plot_data = pgsea_plot_data(),
        pgsea_plot_all_samples_data = pgsea_plot_all_samples_data(),
        go = input$select_go,
        select_org = pre_process$select_org(),
        gene_info = pre_process$all_gene_info(),
        gene_sets = gene_sets()$gene_lists
      )
    })

    # Enrichment Tree -----------
    output$enrichment_tree <- renderPlot({
      req(!is.null(pathway_list_data()))

      enrichment_tree_plot(
        go_table = pathway_list_data(),
        group = "All Groups",
        right_margin = 45
      )
    })

    # Define a Network
    network_data_path <- reactive({
      req(!is.null(pathway_list_data()))
      req(input$up_down_reg_deg)
      network_data(
        network = pathway_list_data(),
        up_down_reg_deg = input$up_down_reg_deg,
        wrap_text_network_deg = input$wrap_text_network_deg,
        layout_vis_deg = input$layout_vis_deg,
        edge_cutoff_deg = input$edge_cutoff_deg
      )
    })

    # Interactive vis network plot
    output$vis_network_path <- visNetwork::renderVisNetwork({
      req(!is.null(network_data_path()))

      vis_network_plot(
        network_data = network_data_path()
      )
    })


    # Markdown report------------
    output$report <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "pathway_report.html",
      content = function(file) {
        # Set up parameters to pass to Rmd document
        params <- list(
          pre_processed = pre_process$data(),
          sample_info = pre_process$sample_info(),
          all_gene_info = pre_process$all_gene_info(),
          deg = deg,
          idep_data = idep_data,
          converted = pre_process$converted(),
          all_gene_names = pre_process$all_gene_names(),
          go = input$select_go,
          select_org = pre_process$select_org(),
          my_range = c(input$min_set_size, input$max_set_size),
          select_contrast = input$select_contrast,
          min_set_size = input$min_set_size,
          max_set_size = input$max_set_size,
          limma = deg$limma(),
          gene_p_val_cutoff = input$gene_p_val_cutoff,
          gene_sets = gene_sets()$gene_lists,
          absolute_fold = input$absolute_fold,
          pathway_p_val_cutoff = input$pathway_p_val_cutoff,
          n_pathway_show = input$n_pathway_show,
          contrast_samples = contrast_samples(),
          sig_pathways = input$sig_pathways,
          pathway_method = input$pathway_method,
          pathway_list_data = pathway_list_data(),
          date = Sys.Date()
        )

        req(params)
        withProgress(message = "Generating Report", {
          incProgress(0.2)

          # Copy the report file to a temporary directory before processing it, in
          # case we don't have write permissions to the current working dir (which
          # can happen when deployed).
          tempReport <- file.path(tempdir(), "pathway_workflow.Rmd")
          # tempReport
          tempReport <- gsub("\\", "/", tempReport, fixed = TRUE)

          # This should retrieve the project location on your device:
          # "C:/Users/bdere/Documents/GitHub/idepGolem"
          wd <- getwd()

          markdown_location <- paste0(wd, "/vignettes/Reports/pathway_workflow.Rmd")
          file.copy(from = markdown_location, to = tempReport, overwrite = TRUE)

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
# mod_08_pathway_ui("08_pathway_ui_1")

## To be copied in the server
# mod_08_pathway_server("08_pathway_ui_1")
