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
        tabsetPanel(
          id = ns("k_means"),

          tabPanel(
            title = "Cluster Map",
            plotOutput(ns("k_means_heatmap")),
            br(),
            h4("Enriched pathways for each cluster"),
            tableOutput(outputId = ns("k_means_go"))
          ),
          tabPanel(
            title = "N-Clusters",
            h4("Determining the number of clusters (k)"),
            h5(
              "Following the elbow method, one should choose k so that adding another 
               cluster does not substantially reduce the within groups sum of squares.",
              a(
                "Wikipedia",
                href = "https://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set",
                target = "_blank"
              )
            ),
            plotOutput(outputId = ns("n_clusters"))
          ),
          tabPanel(
            title = "t-SNE",
            h5(
              "We use the dimension reduction algorith ",
              a(
                "t-SNE",
                href="https://lvdmaaten.github.io/tsne/",
                target="_blank"
              ), 
              "to map the top genes. Examine the distribution
               can help choose the nubmer of clusters in k-Means."
            ),
            fluidRow(
              column(
                width = 4,
                checkboxInput(
                  inputId = ns("color_genes"),
                  label = "Color genes by the results of k-Means",
                  value = TRUE
                ),
              ),
              column(
                width = 4,
                actionButton(
                  inputId = ns("seed_tsne"),
                  label = "Re-calculate using different random numbers"
                )
              )
            ),
            plotOutput(outputId = ns("t_sne_gene_plot"))
          ),
          tabPanel(
            title = "Enrichment",
            h5(
              "Gene sets closer on the tree share more genes. Sizes
              of dot correspond to adjusted p-values"
            ),
            plotOutput(outputId = ns("enrichment_plot_k_means"))
          ),
          tabPanel(
            title = "Enriched Motifs",
            h5("Enriched TF binding motifs in promoters of k-means clusters"),
            br(),
            radioButtons(
              inputId = ns("radio_promoter_k_means"), 
              label    = NULL, 
              choices  = list(
                "Upstream 300bp as promoter" = 300, 
                "Upstream 600bp as promoter" = 600
              ),
              selected = 300
            ),
            tableOutput(outputId = ns("k_means_promoter"))
          )
        )
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
