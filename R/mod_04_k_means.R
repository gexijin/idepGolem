#' 04_k_means UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_04_k_means_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    "k-Means",
    sidebarLayout(
      sidebarPanel(
        sliderInput(
          inputId = ns("n_genes_kmn"),
          label = h4("Most variable genes to include:"),
          min = 0,
          max = 12000,
          value = c(0, 100),
          step = 50
        ),
        sliderInput(
          inputId = ns("n_clusters"),
          label = h4("Number of Clusters"),
          min   = 2,
          max   = 20,
          value = 4,
          step  = 1
        ),
        actionButton(
          inputId = ns("k_means_re_run"),
          label = "Re-Run"
        ),
        actionButton(
          inputId = ns("n_clusters"),
          label = "How many clusters?"
        ),
        actionButton(
          inputId = ns("show_gene_sd"),
          label = "Gene SD distribution"
        ),
        actionButton(
          inputId = ns("gene_tsne"),
          label = "t-SNE map"
        ),
        selectInput(
          inputId = ns("k_means_normalization"),
          label = h5("Normalize by gene:"),
          choices = list(
            "Mean center" = "geneMean",
            "Standardization" = "geneStandardization",
            "L1 Norm" = "L1Norm"
          ),
          selected = "Standardization"
        ),
        tags$style(
          type = "text/css",
          "#k_means-k_means_normalization {width:100%; margin-top:-9px}"
        ),
        actionButton(
          inputId = ns("show_motif_k_means"),
          label = "Enriched TF binding motifs"
        ),
        br(),
        downloadButton(
          outputId = ns("download_data_k_means"),
          label = "K-means data"
        ),
        HTML(
          '<hr style="height:1px;border:none;color:#333;
            background-color:#333;" />'
        ),
        h5("Pathway database"),
        htmlOutput(
          outputId = ns("select_go_3")
        ),
        tags$style(
          type = "text/css",
          "#k_means-select_go_3 { width:100%;   margin-top:-9px}"
        ),
        checkboxInput(
          inputId = ns("remove_redudant_sets"),
          label = "Remove redudant genesets",
          value = TRUE
        ),
        actionButton(
          inputId = ns("modal_enrichment_plot_k_means"),
          label = "Visualize enrichment"
        ),
        downloadButton(
          outputId = ns("download_k_means_go"),
          label = "Enrichment details"
        ),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/k-means/",
          target = "_blank"
        )
      ),
      mainPanel(
        NULL
      )
    )
  )
}

#' 04_k_means Server Functions
#'
#' @noRd
mod_04_k_means_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Update Slider Input ---------
    #observe({
      #req(tab() == "k-Means")
      #req(!is.null(heatmap$data()))

      #if (nrow(heatmap$data()) > 12000) {
        #max_genes <- 12000
      #} else {
        #max_genes <- round(nrow(heatmap$data()) + 50, -2)
      #}
      #updateSliderInput(
        #inputId = "n_genes_kmn",
        #value = c(0, 100),
        #max = max_genes
      #)
    #})
  })
}

## To be copied in the UI
# mod_04_k_means_ui("04_k_means_ui_1")

## To be copied in the server
# mod_04_k_means_server("04_k_means_ui_1")
