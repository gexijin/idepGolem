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
        actionButton(
          inputId = ns("submit_pathway_button"),
          label = "Submit",
          style = "float:right"
        ),
        tippy::tippy_this(
          ns("submit_pathway_button"),
         "Run pathway analysis",
          theme = "light-border"
        ),
        tags$head(tags$style(
          "#pathway-submit_pathway_button{font-size: 16px;color: red}"
        )),
        br(),
        br(),
        htmlOutput(
          outputId = ns("list_comparisons_pathway")
        ),
        selectInput(
          inputId = ns("pathway_method"),
          label = "Pathway analysis method:",
          choices = list(
            "GAGE" = 1,
            "GSEA (preranked fgsea)" = 3,
            "PGSEA" = 2,
            "PGSEA w/ all samples" = 4,
            "ReactomePA" = 5,
            "GSVA" = 6,
            "ssGSEA" = 7,
            "PLAGE" = 8
          ),
          selected = 3
        ),

        htmlOutput(
          outputId = ns("select_go_selector")
        ),

        numericInput(
          inputId = ns("pathway_p_val_cutoff"),
          label = "Pathway signifiance cutoff (FDR):",
          value = 0.1,
          min = 1e-20,
          max = 1,
          step = .05
        ),
        tags$style(
          type = "text/css",
          "#pathway-pathway_p_val_cutoff { width:100%;}"
        ),

        checkboxInput(ns("customize_button"), strong("More options")),

        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = ns("min_set_size"),
              label = "Geneset size: Min.",
              min = 2,
              max = 30,
              value = 5,
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
          "#pathway-min_set_size { width:100%;   margin-top:-10px}"
        ),
        tags$style(
          type = "text/css",
          "#pathway-max_set_size { width:100%;   margin-top:-10px}"
        ),
        numericInput(
          inputId = ns("n_pathway_show"),
          label = "Number of top pathways to show",
          value = 20,
          min = 5,
          max = 100,
          step = 5
        ),
        tags$style(
          type = "text/css",
          "#pathway-n_pathway_show { width:100%;   margin-top:-12px}"
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
        checkboxInput(
          inputId = ns("absolute_fold"),
          label = "Use absolute values of fold changes for GSEA and GAGE",
          value = FALSE
        ),
        checkboxInput(
          inputId = ns("show_pathway_id"),
          label = "Show pathway IDs in results",
          value = FALSE
        ),
        tippy::tippy_this(
          ns("show_pathway_id"),
          "If selected, pathway IDs, such as Path:mmu04115 and GO:0042770,  will be appended to pathway name.",
          theme = "light-border"
        ),
        # Download report button
        downloadButton(
          outputId = ns("report"),
          label = "Report"
        ),
        tippy::tippy_this(
          ns("report"),
          "Generate HTML report of pathway tab",
          theme = "light-border"
        ),
        h6("Beware of P-hacking! If you try all the combinations, you can find evidence for anything."),
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
              condition = "input.submit_pathway_button == 0",
              br(),
              br(),
              h3("Adjust parameters and click the Submit button to perform analysis."),
              ns = ns
            ),
            htmlOutput(
              outputId = ns("main_pathway_result")
            ),
            downloadButton(
              outputId = ns("download_sig_paths"),
              label = "CSV file"
            ),
            tippy::tippy_this(
              ns("download_sig_paths"),
              "Download Significant Pathways",
              theme = "light-border"
            ),
          ),

          tabPanel(
            title = "Tree",
            br(),
            selectInput(
              inputId = ns("leaf_colors"),
              label = "Select leaf colors (low-high)",
              choices = "green-red"
            ),
            plotOutput(
              outputId = ns("enrichment_tree"),
              width = "100%"
            ),
            br(),
            p("Adjusting the width of the browser
            window can render figure differently and
            resolve the \"Figure margin too wide\" error. "),
            br(),
            ottoPlots::mod_download_figure_ui(ns("download_pathway_tree"))
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
                conditionalPanel( #up and down only applies to GSEA and GAGE
                                                  #GAGE                 GSEA
                  condition = "input.pathway_method == 1 | input.pathway_method == 3",
                  selectInput(
                    inputId = ns("up_down_reg_deg"),
                    label = NULL,
                    choices = c(
                      "Both Up & Down" = "All Groups",
                      "Up regulated" = "Up",
                      "Down regulated" = "Down"
                    )
                  ),
                  ns = ns
                )
              )
            ),
            conditionalPanel(
              condition = "input.up_down_reg_deg == 'All Groups'",
              h6(
                "Two pathways (nodes) are connected if they share 30% (default, adjustable) or more genes.
                Green and red represents up- and down-regulated pathways, respectively. You can move the nodes by
                dragging them, zoom in and out by scrolling, and shift the entire network by click on an
                empty point and drag. Darker nodes are more significantly enriched gene sets. Bigger nodes
                represent larger gene sets. Thicker edges represent more overlapped genes."
              ),
              ns = ns
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
            br(),
            conditionalPanel(
              condition = "(input.pathway_method == 1 | input.pathway_method == 2 |
                            input.pathway_method == 3 | input.pathway_method == 4
                            | input.pathway_method == 6 | input.pathway_method == 7
                            | input.pathway_method == 8) 
                            & input.select_go != 'KEGG'",
              h5("Please select KEGG database, if available, from left and perform pathway analysis first."),
              ns = ns
            ),
            conditionalPanel(
              condition = "(input.pathway_method == 1 | input.pathway_method == 2 |
                            input.pathway_method == 3 | input.pathway_method == 4 
                            | input.pathway_method == 6 | input.pathway_method == 7
                            | input.pathway_method == 8) &
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

              br(),
              h5("Red and green (or orange and blue) represent up- and down-
                 regulated genes, respectively."),
              downloadButton(ns('download_kegg'),''),
              tippy::tippy_this(
                ns("download_kegg"),
                "Download KEGG plot",
                theme = "light-border"
              ),
              br(),

              ns = ns
            )
          ),
          tabPanel(
            title = "Info",
            includeHTML(app_sys("app/www/pathway.html"))
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

      #This make it reactive: every time method changes, it reset to KEGG.
      #req(input$pathway_method != 5)

      # if there is KEGG, use KEGG as default
      selected <- "GOBP"
      if ("KEGG" %in% pre_process$gmt_choices()) {
        selected <- "KEGG"
      }
      selectInput(
        inputId = ns("select_go"),
        label = "Pathway database:",
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
        tab() == "DEG2" ||
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

    observe({
      shinyjs::toggle(id = "max_set_size", condition = input$customize_button)
      shinyjs::toggle(id = "min_set_size", condition = input$customize_button)
      shinyjs::toggle(id = "n_pathway_show", condition = input$customize_button)
      shinyjs::toggle(id = "gene_p_val_cutoff", condition = input$customize_button)
      shinyjs::toggle(id = "absolute_fold", condition = input$customize_button)
      shinyjs::toggle(id = "show_pathway_id", condition = input$customize_button)
      shinyjs::toggle(id = "absolute_fold", condition = input$customize_button)
    })

    output$main_pathway_result <- renderUI({
      req(input$submit_pathway_button)

      isolate({

        # use the switch function to deliver results for different methods
        switch(
          as.integer(input$pathway_method),    
          
          #1 GAGE
          tableOutput(ns("gage_pathway_table")),

          #2 PGSEA
          tagList(
            h5("Red and blue indicates relatively activated
              and suppressed pathways, respectively.
              GS just indicates a color scale."),
            plotOutput(
              outputId = ns("pgsea_plot"),
              inline = TRUE
            )
          ),

          #3 FGSEA
          tableOutput(outputId = ns("fgsea_pathway")),

          #4 PGSEA-All
          tagList(
            h5("Red and blue indicates relatively activated
              and suppressed pathways, respectively.
              GS just indicates a color scale."),
            plotOutput(
              outputId = ns("pgsea_plot_all_samples"),
              inline = TRUE
            )
          ),

          #5 ReactomePA
          DT::dataTableOutput(outputId = ns("reactome_pa_pathway")),

          #6 GSVA
          tagList(
            h5("Red and blue indicates relatively activated
              and suppressed pathways, respectively.
              GS just indicates a color scale."),
            plotOutput(
              outputId = ns("gsva_plot"),
              inline = TRUE
            )
          ),

          #7 ssGSEA   same as above
          tagList(
            h5("Red and blue indicates relatively activated
              and suppressed pathways, respectively.
              GS just indicates a color scale."),
            plotOutput(
              outputId = ns("gsva_plot"),
              inline = TRUE
            )
          ),

          #8 PLAGE   same as above
          tagList(
            h5("Red and blue indicates relatively activated
              and suppressed pathways, respectively.
              GS just indicates a color scale."),
            plotOutput(
              outputId = ns("gsva_plot"),
              inline = TRUE
            )
          ),
        )
      })
    })
    
    output$download_sig_paths <- downloadHandler(
      filename = function() {
        "sig_pathways.csv"
      },
      content = function(file) {
        write.csv(res_pathway()[,c(1,6,7,3:5)], file)
      }
    )

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
            "Select a comparison:",
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
      } else if (input$pathway_method >= 6 && input$pathway_method <= 8 ) {
        if (!is.null(gsva_plot_data())) {
          if (dim(gsva_plot_data())[2] > 1) {
            pathways <- as.data.frame(gsva_plot_data())
            choices <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
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
      } else if (input$pathway_method >= 6 && input$pathway_method <= 8) {
        if (!is.null(gsva_plot_data())) {
          if (dim(gsva_plot_data())[2] > 1) {
            pathways <- as.data.frame(gsva_plot_data())
            choices <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
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
    gene_sets <- eventReactive(input$submit_pathway_button, {
#    gene_sets <- reactive({
      #      req(tab() == "Pathway")
      req(!is.null(input$select_go))

      withProgress(message = "Reading pathway database", {
        incProgress(0.2)

        if (
          !is.null(pre_process$gmt_file())
        ) {
          in_file <- pre_process$gmt_file()
          in_file <- in_file$datapath
          gene_lists <- read_gmt_robust(in_file)
          pathway_info <- data.frame(
            id = 1:length(gene_lists),
            description = names(gene_lists),
            memo = rep("", length(gene_lists))
          )
          gene_sets <- list(
            gene_lists = gene_lists,
            pathway_info = pathway_info
          )     
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

    gage_pathway_data <- eventReactive(input$submit_pathway_button, {
#    gage_pathway_data <- reactive({
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
        input$submit_pathway_button
        isolate({
          req(!is.null(gage_pathway_data()))

          res <- gage_pathway_data()
          if (ncol(res) > 1) {
            # add URL
            ix <- match(res[, 2], gene_sets()$pathway_info$description)

            # remove pathway ID  only in Ensembl species
            if (!input$show_pathway_id && pre_process$select_org() > 0) {
              res[, 2] <- remove_pathway_id(res[, 2], input$select_go)
            }
            res[, 2] <- hyperText(
              res[, 2],
              gene_sets()$pathway_info$memo[ix]
            )
            res$Genes <- as.character(res$Genes)
          }
          return(res)
        })
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
        input$submit_pathway_button
        isolate({
          req(input$pathway_method == 2)
          withProgress(message = "Running PGSEA...", {
            incProgress(0.2)

            # only remove pathway ID for Ensembl species
            show_pathway_id <- input$show_pathway_id
            # always show pathway ID for STRING species
            if (pre_process$select_org() < 0) {
              show_pathway_id <- TRUE
            }

            plot_pgsea(
              my_range = c(input$min_set_size, input$max_set_size),
              processed_data = pre_process$data(),
              contrast_samples = contrast_samples(),
              gene_sets = gene_sets()$gene_lists,
              pathway_p_val_cutoff = input$pathway_p_val_cutoff,
              n_pathway_show = input$n_pathway_show,
              select_go = input$select_go,
              show_pathway_id = show_pathway_id
            )
          })
        })
      },
      height = 800,
      width = 800
    )

    pgsea_plot_data <- eventReactive(input$submit_pathway_button, {
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


    output$gsva_plot <- renderPlot(
      {
        input$submit_pathway_button
        isolate({
          req(input$pathway_method >= 6 && input$pathway_method <= 8)
          gsva_algorithm <- switch(
            as.numeric(input$pathway_method) - 5, 
            "gsva",   #6
            "ssgsea", #7
            "plage"   #8
          )
          withProgress(message = paste("Running", toupper(gsva_algorithm), "..."), {
            incProgress(0.2)

            # only remove pathway ID for Ensembl species
            show_pathway_id <- input$show_pathway_id
            # always show pathway ID for STRING species
            if (pre_process$select_org() < 0) {
              show_pathway_id <- TRUE
            }


            plot_gsva(
              my_range = c(input$min_set_size, input$max_set_size),
              processed_data = pre_process$data(),
              contrast_samples = contrast_samples(),
              gene_sets = gene_sets()$gene_lists,
              pathway_p_val_cutoff = input$pathway_p_val_cutoff,
              n_pathway_show = input$n_pathway_show,
              select_go = input$select_go,
              show_pathway_id = show_pathway_id,
              algorithm = gsva_algorithm
            )
          })
        })
      },
      height = 800,
      width = 800
    )

    gsva_plot_data <- eventReactive(input$submit_pathway_button, {
      req(input$pathway_method >= 6 && input$pathway_method <= 8)
      req(!is.null(gene_sets()))
      gsva_algorithm <- switch(
        as.numeric(input$pathway_method) - 5, 
        "gsva",   #6
        "ssgsea", #7
        "plage"   #8
      )
      withProgress(message = paste("Running", toupper(gsva_algorithm), "..."), {
        incProgress(0.2)
        get_gsva_plot_data(
          my_range = c(input$min_set_size, input$max_set_size),
          data = pre_process$data(),
          select_contrast = input$select_contrast,
          gene_sets = gene_sets()$gene_lists,
          sample_info = pre_process$sample_info(),
          select_factors_model = deg$select_factors_model(),
          select_model_comprions = deg$select_model_comprions(),
          pathway_p_val_cutoff = input$pathway_p_val_cutoff,
          n_pathway_show = input$n_pathway_show,
          algorithm = gsva_algorithm
        )
      })
    })


    fgsea_pathway_data <- eventReactive(input$submit_pathway_button, {
#    fgsea_pathway_data <- reactive({
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
    
    res_pathway <- reactive({
      req(input$pathway_method == 3)
      req(!is.null(fgsea_pathway_data()))
      res <- fgsea_pathway_data()
      # copy name from res[,2]
      pathway_csv_name <- colnames(res)[2]
      non_hypertext_name <- paste(pathway_csv_name, "Pathways")
      
      if (ncol(res) > 1) {
        # add URL
        ix <- match(res[, 2], gene_sets()$pathway_info$description)
        # remove pathway ID, but only in Ensembl species
        if (!input$show_pathway_id && pre_process$select_org() > 0) {
          res[, 2] <- remove_pathway_id(res[, 2], input$select_go)
        }
        # copy res[,2] to new column before hypertext
        res[non_hypertext_name] <- res[,2]
        res[, 2] <- hyperText(
          res[, 2],
          gene_sets()$pathway_info$memo[ix]
        )
        # create separate URL column for download
        res$URL <- NULL
        res$URL <- gene_sets()$pathway_info$memo[ix]
        res$Genes <- as.character(res$Genes)
      }
      return(res)
    })

    output$fgsea_pathway <- renderTable(
      {
        res_pathway()[,1:5]
      },
      digits = -1,
      spacing = "s",
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = TRUE,
      sanitize.text.function = function(x) x
    )



    reactome_pa_pathway_data <- eventReactive(input$submit_pathway_button, {
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
        input$submit_pathway_button
        isolate({
          req(input$pathway_method == 4)
          withProgress(message = "Running PGSEA on all samples ...", {
            incProgress(0.2)

            # only remove pathway ID for Ensembl species
            show_pathway_id <- input$show_pathway_id
            # always show pathway ID for STRING species
            if (pre_process$select_org() < 0) {
              show_pathway_id <- TRUE
            }

            pgsea_plot_all(
              go = input$select_go,
              my_range = c(input$min_set_size, input$max_set_size),
              data = pre_process$data(),
              select_contrast = input$select_contrast,
              gene_sets = gene_sets()$gene_lists,
              pathway_p_val_cutoff = input$pathway_p_val_cutoff,
              n_pathway_show = input$n_pathway_show,
              select_go = input$select_go,
              show_pathway_id = show_pathway_id
            )
          })
        })
      },
      height = 800,
      width = 800
    )

    pgsea_plot_all_samples_data <- eventReactive(input$submit_pathway_button, {
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

      req(!is.null(gene_sets()$gene_lists))
      if(is.null(input$sig_pathways)) {
        return(NULL)
      } else {
        pathway_select_data(
          sig_pathways = input$sig_pathways,
          gene_sets = gene_sets()$gene_lists,
          contrast_samples = contrast_samples(),
          data = pre_process$data(),
          select_org = pre_process$select_org(),
          all_gene_names = pre_process$all_gene_names()
        )
      }
    })

    # Kegg Colors --------
    kegg_colors <- list(
      "Green-Red" = c("green", "red"),
      "Red-Green" = c("red", "green"),
      "Blue-Red" = c("blue", "red"),
      "Green-Magenta" = c("green", "magenta"),
      "Blue-Orange" = c("blue", "orange")
    )
    kegg_choices <- c(
      "Green-Red",
      "Red-Green",
      "Blue-Red",
      "Green-Magenta",
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
      list(
        src = kegg_image(), 
        height = "100%", 
        width = "100%", 
        contentType = "image/png"
      )
      },
      deleteFile = FALSE
    )

    kegg_image <- reactive({
      req(!is.null(input$sig_pathways_kegg))
      
      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Generating KEGG...",
        color = "#000000"
      )
      tmpfile <- kegg_pathway(
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
      shinybusy::remove_modal_spinner()
     tmpfile$src
    })
    
    output$download_kegg <- downloadHandler(
      filename = function() {
        "KEGGplot.png"
      },
      content = function(file) {
        file.copy(kegg_image(), file)
      },
      contentType = "image/png"
    ) 

    # List of pathways with details
    pathway_list_data <- reactive({
      # only remove pathway ID for Ensembl species
      show_pathway_id <- input$show_pathway_id
      # always show pathway ID for STRING species
      if (pre_process$select_org() < 0) {
        show_pathway_id <- TRUE
      }

      get_pathway_list_data(
        pathway_method = input$pathway_method,
        gage_pathway_data = gage_pathway_data(),
        fgsea_pathway_data = fgsea_pathway_data(),
        pgsea_plot_data = pgsea_plot_data(),
        pgsea_plot_all_samples_data = pgsea_plot_all_samples_data(),
        gsva_plot_data = gsva_plot_data(),
        go = input$select_go,
        select_org = pre_process$select_org(),
        gene_info = pre_process$all_gene_info(),
        gene_sets = gene_sets()$gene_lists,
        show_pathway_id = show_pathway_id
      )
    })
    
    observe({
      updateSelectInput(
        session = session,
        inputId = "leaf_colors",
        choices = kegg_choices
      )
    })

    enrichment_tree_p <- reactive({
      req(!is.null(pathway_list_data()))
      enrichment_tree_plot(
        go_table = pathway_list_data(),
        group = "All Groups",
        right_margin = 30,
        leaf_color_choices = kegg_colors[[input$leaf_colors]]
      )
      p <- recordPlot()
      return(p)
    })
    output$enrichment_tree <- renderPlot({
      req(!is.null(enrichment_tree_p()))
      print(enrichment_tree_p())
    })
    download_pathway_tree <- ottoPlots::mod_download_figure_server(
      id = "download_pathway_tree",
      filename = "pathway_tree",
      figure = reactive({
        enrichment_tree_p()
      }),
      label = ""
    )

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
      filename = paste0(
              "pathway_workflow_",
              format(Sys.time(), "%Y-%m-%d_%H-%M"),
              ".html"
            ),
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
          up_down_reg_deg = input$up_down_reg_deg,
          wrap_text_network_deg = input$wrap_text_network_deg,
          layout_vis_deg = input$layout_vis_deg,
          edge_cutoff_deg = input$edge_cutoff_deg,
          selected_pathway_data = selected_pathway_data(),
          heatmap_color_select = pre_process$heatmap_color_select(),
          sig_pathways_kegg = input$sig_pathways_kegg, 
          kegg_color_select = input$kegg_color_select,
          kegg_colors = kegg_colors,
          descr = deg$limma()[["description"]],
          show_pathway_id = input$show_pathway_id
        )

        req(params)
        withProgress(message = "Generating Report", {
          incProgress(0.2)

          # Copy the report file to a temporary directory before processing it, in
          # case we don't have write permissions to the current working dir (which
          # can happen when deployed).

          tempReport <- file.path(
            tempdir(), 
            "pathway_workflow.Rmd"
          )
          # tempReport
          tempReport <- gsub("\\", "/", tempReport, fixed = TRUE)

          wd <- getwd()

          markdown_location <- app_sys("app/www/RMD/pathway_workflow.Rmd")
          #markdown_location <- "C:/work/idepGolem/vignettes/Reports/pathway_workflow.Rmd"
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
