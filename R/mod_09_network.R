#' 09_network UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_09_network_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "Network",
    sidebarLayout(
      sidebarPanel(
        h5(
          "Identify co-expression networks and sub-modules using",
          a(
            "WGCNA.",
            href = "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559",
            target = "_blank"
          ),
          "Only useful when  sample size is large(> 15)."
        ),
        numericInput(
          inputId = ns("n_genes_network"),
          label = "Most variable genes to include (< 3001)",
          min = 10,
          max = 3000,
          value = 1000
        ),
        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = ns("soft_power"),
              label = "Soft Threshold",
              min = 1,
              max = 20,
              value = 5
            )
          ),
          column(
            width = 6,
            numericInput(
              inputId = ns("min_module_size"),
              label = "Min. Module Size",
              min = 10,
              max = 100,
              value = 20
            )
          )
        ),
        tags$style(
          type = "text/css",
          "#network-soft_power{ width:100%;   margin-top:-12px}"
        ),
        tags$style(
          type = "text/css",
          "#network-min_module_size{ width:100%;   margin-top:-12px}"
        ),
        uiOutput(outputId = ns("heatmap_color_ui")),
        HTML(
          "<hr style='height:1px;border:none;color:#333;background-color:#333;' />"
        ),
        htmlOutput(outputId = ns("list_wgcna_modules")),
        downloadButton(outputId = ns("download_all_WGCNA_module"),"Download all modules"),
        br(),
        downloadButton(outputId = ns("download_selected_WGCNA_module"),"Download network for selected module"),
        textOutput(ns("module_statistic")),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/network/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          id = ns("network_tabs"),
          tabPanel(
            "Network Plot",
            br(),
            fluidRow(
              column(
                width = 4,
                numericInput(
                  inputId = ns("edge_threshold"),
                  label = "Edge Threshold",
                  min = 0,
                  max = 1,
                  value = .4,
                  step = .1
                )
              ),
              column(
                width = 4,
                numericInput(
                  inputId = ns("top_genes_network"),
                  label = "Top genes",
                  min = 10,
                  max = 2000,
                  value = 10,
                  step = 10
                )
              ),
              column(
                width = 4,
                style = "margin-top: 25px;",
                actionButton(
                  inputId = ns("network_layout"),
                  label = "Change network layout",
                  style = "float:center"
                )
              )
            ),
            tags$style(
              type = "text/css",
              "#network-edge_threshold{ width:100%;   margin-top:-12px}"
            ),
            tags$style(
              type = "text/css",
              "#network-top_genes_network{ width:100%;   margin-top:-12px}"
            ),
            br(),
            plotOutput(outputId = ns("module_network"))
            #            ,h5(
            #              "The network file can be imported to",
            #              a("VisANT", href = "http://visant.bu.edu/", target = "_blank"),
            #              " or ",
            #              a("Cytoscape.", href = "http://www.cytoscape.org/", target = "_blank")
            #            ),
          ),
          tabPanel(
            "Module Plot",
            plotOutput(
              outputId = ns("module_plot"),
              width = "100%",
              height = "500px"
            )
          ),
          tabPanel(
            title = "Heatmap",
            mod_12_heatmap_ui(ns("12_heatmap_1"))
          ),
          tabPanel(
            "Enrichment",
            p("Enriched pathways in the selected module"),
            mod_11_enrichment_ui(ns("enrichment_table_cluster"))
          ),
          tabPanel(
            "Scale Independence",
            br(),
            plotOutput(
              outputId = ns("scale_independence_plot"),
              width = "100%",
              height = "500px"
            )
          ),
          tabPanel(
            "Mean Connectivity",
            br(),
            plotOutput(
              outputId = ns("mean_connectivity_plot"),
              width = "100%",
              height = "500px"
            )
          )
        )
      )
    )
  )
}

#' 09_network Server Functions
#'
#' @noRd
mod_09_network_server <- function(id, pre_process, idep_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Interactive heatmap environment
    network_env <- new.env()

    output$list_wgcna_modules <- renderUI({
      req(!is.null(wgcna()))
      module_list <- get_wgcna_modules(wgcna = wgcna())
      req(!is.null(module_list))
      selectInput(
        inputId = ns("select_wgcna_module"),
        label = "Select a module",
        choices = module_list
      )
    })

    wgcna <- reactive({
      req(!is.null(pre_process$data()))
      withProgress(message = "Runing WGCNA ...", {
        incProgress(0.2)
        get_wgcna(
          data = pre_process$data(),
          n_genes = input$n_genes_network,
          soft_power = input$soft_power,
          min_module_size = input$min_module_size
        )
      })
    })
    module_csv_data <- reactive({
      req(!is.null(wgcna()))

      prepare_module_csv(
        wgcna = wgcna(), select_org = pre_process$select_org(),
        all_gene_info = pre_process$all_gene_info()
      )
    })
    output$download_all_WGCNA_module <- downloadHandler(
      filename = function() {
        "WGCNA_modules.csv"
      },
      content = function(file) {
        write.csv(module_csv_data(), file)
      }
    )
    module_csv_data_filter <- reactive({
      req(!is.null(module_csv_data()))

      prepare_module_csv_filter(
        module_data = module_csv_data(), select_org = pre_process$select_org()
      )
    })
    output$download_selected_WGCNA_module <- downloadHandler(
      filename <- function() {
        paste0("module_", input$select_wgcna_module, ".csv")
      },
      content <- function(file) {
        write.csv(module_csv_data(), file)
      },
    )
    
    output$module_plot <- renderPlot({
      req(!is.null(wgcna()))
      get_module_plot(wgcna())
    })

    network <- reactiveValues(network_plot = NULL)

    observe({
      req(!is.null(input$select_wgcna_module))
      req(!is.null(wgcna()))

      network$network_plot <- get_network_plot(
        select_wgcna_module = input$select_wgcna_module,
        wgcna = wgcna(),
        top_genes_network = input$top_genes_network,
        select_org = pre_process$select_org(),
        all_gene_info = pre_process$all_gene_info(),
        edge_threshold = input$edge_threshold
      )
    })

    observeEvent(
      input$network_layout,
      {
        req(!is.null(input$select_wgcna_module))
        req(!is.null(wgcna()))

        network$network_plot <- get_network_plot(
          select_wgcna_module = input$select_wgcna_module,
          wgcna = wgcna(),
          top_genes_network = input$top_genes_network,
          select_org = pre_process$select_org(),
          all_gene_info = pre_process$all_gene_info(),
          edge_threshold = input$edge_threshold
        )

        
      }
    )

    output$module_network <- renderPlot({
      req(!is.null(input$select_wgcna_module))
      req(!is.null(wgcna()))
      network$network_plot()
    })

    network_query <- reactive({
      req(!is.null(input$select_wgcna_module))
      req(!is.null(wgcna()))

      network_query <- network_enrich_data(
        select_wgcna_module = input$select_wgcna_module,
        wgcna = wgcna()
      )
    })

    gene_lists <- reactive({
      req(!is.null(network_query()))
      withProgress(message = "Running Analysis", {
        incProgress(0.2)
        gene_lists <- list()
        gene_lists[["Cluster"]] <- dplyr::filter(
          pre_process$all_gene_names(),
          ensembl_ID %in% network_query()
        )
      })
      return(gene_lists)
    })


    enrichment_table_cluster <- mod_11_enrichment_server(
      id = "enrichment_table_cluster",
      gmt_choices = reactive({
        pre_process$gmt_choices()
      }),
      gene_lists = reactive({
        gene_lists()
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

    output$scale_independence_plot <- renderPlot({
      req(!is.null(wgcna()))

      p <- plot_scale_independence(
        wgcna = wgcna()
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })

    output$mean_connectivity_plot <- renderPlot({
      req(!is.null(wgcna()))

      p <- plot_mean_connectivity(
        wgcna = wgcna()
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })

    output$module_statistic <- renderText({
      req(!is.null(wgcna()))
      paste(
        "A network of",
        wgcna()$n_genes,
        "genes was divided into ",
        wgcna()$n_modules,
        "modules."
      )
    })


    # Remove messages if the tab changes --------
    observe({
      req(tab() != "Network")
      removeNotification("network_summary")
    })

    network_data <- reactive({
      req(!is.null(network_query()))

      data <- pre_process$data()[rownames(pre_process$data()) %in% network_query(), ]
    })

    heatmap_module <- mod_12_heatmap_server(
      id = "12_heatmap_1",
      data = reactive({
        network_data()
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
  })
}

## To be copied in the UI
# mod_09_network_ui("09_network_ui_1")

## To be copied in the server
# mod_09_network_server("09_network_ui_1")
