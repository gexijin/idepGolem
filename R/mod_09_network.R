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
          "WGCNA requires at least 15 samples to be meaningful."
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
        uiOutput(outputId = ns("heatmap_color_ui")),
        HTML(
          "<hr style='height:1px;border:none;color:#333;background-color:#333;' />"
        ),
        htmlOutput(outputId = ns("list_wgcna_modules")),
        div(
          style = "display: flex; gap: 10px;",
          downloadButton(
            outputId = ns("download_all_WGCNA_module"),
            "All modules"
          ),
          downloadButton(
            outputId = ns("download_selected_WGCNA_module"),
            "Selected module"
          ),
          conditionalPanel(
            condition = "input.network_tabs == 'Heatmap'",
            downloadButton(
              outputId = ns("download_heat_data"),
              "Heatmap Data"
            ),
            ns = ns
          )
        ),
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
            br(),
            plotOutput(outputId = ns("module_network")),
            downloadButton(outputId = ns("download_module_network"), "Network file"),
            tippy::tippy_this(
              ns("download_module_network"),
              "This file can be imported to CytoScape or VisANT for further analysis.",
              theme = "light-border"
            ),
            #ottoPlots::mod_download_figure_ui(
            #  id = ns("dl_network_plot")
            #)
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
      req(!is.na(input$n_genes_network))
      req(!is.na(input$min_module_size))
      req(!is.na(input$soft_power))
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
        write.csv(module_csv_data(), file, row.names = FALSE)
      }
    )
    module_csv_data_filter <- reactive({
      req(!is.null(module_csv_data()))

      prepare_module_csv_filter(
        module_data = module_csv_data(), 
        module = input$select_wgcna_module
      )
    })
    output$download_selected_WGCNA_module <- downloadHandler(
      filename <- function() {
        paste0("module_", input$select_wgcna_module, ".csv")
      },
      content <- function(file) {
        write.csv(module_csv_data_filter(), file, row.names = FALSE)
      }
    )
    
    output$module_plot <- renderPlot({
      req(!is.null(wgcna()))
      get_module_plot(wgcna())
    })

    network <- reactiveValues(network_plot = NULL)

    adj_matrix <- reactive({
      req(!is.null(input$select_wgcna_module))
      req(!is.na(input$network_layout))
      req(!is.na(input$edge_threshold))
      req(!is.na(input$top_genes_network))
      req(!is.na(input$n_genes_network))
      req(!is.na(input$min_module_size))
      req(!is.na(input$soft_power))
      req(!is.null(wgcna()))
      
      tem <- input$network_layout
      tem <- input$edge_threshold
      get_network(
        select_wgcna_module = input$select_wgcna_module,
        wgcna = wgcna(),
        top_genes_network = input$top_genes_network,
        select_org = pre_process$select_org(),
        all_gene_info = pre_process$all_gene_info(),
        edge_threshold = input$edge_threshold
      )
    })

    observe({
      req(!is.null(input$select_wgcna_module))
      req(!is.na(input$network_layout))
      req(!is.na(input$edge_threshold))
      req(!is.na(input$top_genes_network))
      req(!is.na(input$n_genes_network))
      req(!is.na(input$min_module_size))
      req(!is.na(input$soft_power))
      req(!is.null(wgcna()))

      tem = input$network_layout
      network$network_plot <- get_network_plot(
        adj_matrix(),
        edge_threshold = input$edge_threshold
      )
    })



    output$download_module_network <- downloadHandler(
      filename = function() {
        paste0("module_network_", input$select_wgcna_module, ".csv")
      },
      content = function(file) {
        # convert adjacency matrix to edge list, i.e. from wide to long format
        network <- reshape2::melt(adj_matrix(), id.vars = "gene")
        network <- dplyr::rename(network, gene1 = Var1, gene2 = Var2, weight = value)
        network <- dplyr::filter(network, weight > 0 & gene1 != gene2)
        write.csv(network, file, row.names = FALSE)
      }
    )

    output$module_network <- renderPlot({
      req(!is.null(input$select_wgcna_module))
      req(!is.null(wgcna()))
      network$network_plot()
    })

    # not working
    #dl_network_plot <- ottoPlots::mod_download_figure_server(
    #  id = "dl_network_plot",
    #  filename = "module_network",
    #  figure = reactive({
    #    network$network_plot
    #  })
    #  ,
    #  label = ""
    #)

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
    
    output$download_heat_data <- downloadHandler(
      filename = function() {
        req(!is.null(input$select_wgcna_module))
        # trim module name to the color
        paste0(
          gsub(
            "\\d|\\s|genes|\\.|\\(|\\)", 
            "", 
            input$select_wgcna_module, 
            perl = TRUE
          ),
          "_Heatmap_Data.csv"
        )
      },
      content = function(file) {
        req(!is.null(network_data()))
        
        df <- network_data()
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
  })
}

## To be copied in the UI
# mod_09_network_ui("09_network_ui_1")

## To be copied in the server
# mod_09_network_server("09_network_ui_1")
