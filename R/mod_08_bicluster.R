#' 08_bicluster UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_08_bicluster_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "Bicluster",
    sidebarLayout(
      sidebarPanel(
        h5(
          "Biclustering can discover genes correlated on subset of samples.
           Only useful when sample size is large(N>15)
           and more than 2 sample groups.
           Based on the",
          a("biclust",
            href = "https://cran.r-project.org/web/packages/biclust/index.html"
          ),
          " and ",
          a("QUBIC",
            href = "https://www.bioconductor.org/packages/release/bioc/html/QUBIC.html"
          ),
          " R packages."
        ),
        numericInput(
          inputId = ns("n_genes"),
          label = "Most variable genes to include: ",
          min = 10,
          max = 2000,
          value = 1000
        ),
        selectInput(
          inputId = ns("biclust_method"),
          label = "Method:",
          choices = list(
            "BCCC" = "biclust::BCCC()",
            "QUBIC" = "QUBIC::BCQU()",
            "runibic" = "runibic::BCUnibic()",
            "BCXmotifs" = "biclust::BCXmotifs()",
            "BCPlaid" = "biclust::BCPlaid()",
            "BCSpectral" = "biclust::BCSpectral()",
            "BCBimax" = "biclust::BCBimax()",
            "BCQuest" = "biclust::BCQuest()"
          ),
          selected = "BCCC()"
        ),
        htmlOutput(outputId = ns("list_biclusters")),
        textOutput(ns("bicluster_info")),
        a(
          h5(
            "Questions?",
            align = "right"
          ),
          href = "https://idepsite.wordpress.com/biclustering/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            title = "Heatmap",
            mod_12_heatmap_ui(ns("12_heatmap_1"))
          ),
          tabPanel(
            "Enrichment",
            p("Enriched pathways in the selected cluster:"),
            mod_11_enrichment_ui(ns("enrichment_table_cluster"))
          ),
          tabPanel(
            "Genes",
            DT::dataTableOutput(
              outputId = ns("gene_list_bicluster"),
            ),
            br(),
            uiOutput(ns("download_biclust_button"))
          )
        )
      )
    )
  )
}





#------------------------------------------------------------------------



#' 08_bicluster Server Functions
#'
#' @noRd
mod_08_bicluster_server <- function(id, pre_process, idep_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Interactive heatmap environment
    biclust_env <- new.env()

    # all clusters
    biclustering <- reactive({
      req(!is.null(pre_process$data()))

      get_biclustering(
        data = pre_process$data(),
        n_genes = input$n_genes,
        biclust_method = input$biclust_method
      )
    })

    output$list_biclusters <- renderUI({
      req(!is.null(biclustering()))
      req(biclustering()$res@Number != 0)

      selectInput(
        inputId = ns("select_bicluster"),
        label = "Select a cluster",
        selected = 1,
        choices = 1:biclustering()$res@Number
      )
    })

    # information for a specific cluster
    biclust_data <- reactive({
      req(!is.null(biclustering()) && !is.null(input$select_bicluster))
      req(biclustering()$res@Number != 0)
      withProgress(message = "Runing biclustering", {
        incProgress(0.2)
        return(
          biclust::bicluster(
            biclustering()$data,
            biclustering()$res,
            as.numeric(input$select_bicluster)
          )[[1]]
        )
      })
    })

    # split "green-white-red" to c("green", "white", "red")
    heatmap_color_select <- reactive({
      req(pre_process$heatmap_color_select())
      unlist(strsplit(pre_process$heatmap_color_select(), "-"))
    })

    heatmap_module <- mod_12_heatmap_server(
      id = "12_heatmap_1",
      data = reactive({
        biclust_data()
      }),
      bar = function() {
        return(NULL)
      },
      all_gene_names = reactive({
        pre_process$all_gene_names()
      }),
      cluster_rows = TRUE,
      heatmap_color = reactive({
        heatmap_color_select()
      }),
      select_gene_id = reactive({
        pre_process$select_gene_id()
      })
    )

    # Biclustering summary message -----------
    output$bicluster_info <- renderText({
      req(!is.null(biclustering()) && !is.null(input$select_bicluster))

      bicluster_summary_message(
        biclustering = biclustering(),
        select_bicluster = input$select_bicluster
      )
    })

    gene_lists <- reactive({
      req(!is.null(biclust_data()))

      withProgress(message = "Generating gene lists", {
        incProgress(0.1)
        gene_names <- merge_data(
          all_gene_names = pre_process$all_gene_names(),
          data = biclust_data(),
          merge_ID = "ensembl_ID"
        )
        incProgress(0.3)
        gene_lists <- list()
        # Only keep the gene names and scrap the data
        gene_lists[["Cluster"]] <- dplyr::select_if(gene_names, is.character)
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
      })
    )

    # list of genes and symbols in the currently selected cluster
    genes_in_selected_cluster <- reactive({
      req(!is.null(biclust_data()))
      req(!is.null(biclustering()) && !is.null(input$select_bicluster))

      get_biclust_table_data(
        res = biclustering()$res,
        biclust_data = biclust_data(),
        select_org = pre_process$select_org(),
        all_gene_info = pre_process$all_gene_info()
      )
    })

    output$gene_list_bicluster <- DT::renderDataTable({
      req(!is.null(biclust_data()) && !is.null(genes_in_selected_cluster()))
      req(!is.null(biclustering()) && !is.null(input$select_bicluster))

      DT::datatable(
        genes_in_selected_cluster(),
        options = list(
          pageLength = 20,
          scrollX = "400px"
          # dom = 'ftg' # hide "Show 20 entries"
        ),
        class = "cell-border stripe",
        rownames = FALSE
      )
    })

    output$download_biclust <- downloadHandler(
      filename = "bicluster.csv",
      content = function(file) {
        write.csv(genes_in_selected_cluster(), file, row.names = FALSE)
      }
    )

    output$download_biclust_button <- renderUI({
      req(!is.null(biclust_data()) && !is.null(genes_in_selected_cluster()))
      req(!is.null(biclustering()) && !is.null(input$select_bicluster))
      downloadButton(
        outputId = ns("download_biclust"),
        "All genes"
      )
    })
  })
}
