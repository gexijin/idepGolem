#' 03_heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_03_heatmap_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Heatmap",
    sidebarLayout(

      # Heatmap Panel Sidebar ----------
      sidebarPanel(
        sliderInput(
          inputId = ns("n_genes"),
          label = h4("Most variable genes to include:"),
          min = 0,
          max = 12000,
          value = c(0, 100),
          step = 50
        ),
        HTML(
          '<hr style="height:1px;border:none;
           color:#333;background-color:#333;" />'
        ),

        # Gene ID Selection -----------
        selectInput(
          inputId = ns("select_gene_id"),
          label = "Select Gene ID Label (<= 50 genes)",
          choices = NULL,
          selected = NULL
        ),

        # Heatmap customizing features ----------
        strong("Customize hierarchical clustering (Default values work well):"),
        fluidRow(
          column(width = 3, h5("Color")),
          column(
            width = 9,
            selectInput(
              inputId = ns("heatmap_color_select"),
              label = NULL,
              choices = "green-black-red",
              width = "100%"
            )
          )
        ),
        fluidRow(
          column(width = 4, h5("Distance")),
          column(
            width = 8,
            selectInput(
              inputId = ns("dist_function"),
              label = NULL,
              choices = NULL,
              width = "100%"
            )
          )
        ),
        fluidRow(
          column(width = 4, h5("Linkage")),
          column(
            width = 8,
            selectInput(
              inputId = ns("hclust_function"),
              label = NULL,
              choices = c(
                "average", "complete", "single",
                "median", "centroid", "mcquitty"
              ),
              width = "100%"
            )
          )
        ),
        fluidRow(
          column(width = 8, h5("Cut-off Z score")),
          column(
            width = 4,
            numericInput(
              inputId = ns("heatmap_cutoff"),
              label = NULL,
              value = 4,
              min = 2,
              step = 1
            )
          )
        ),

        # Checkbox features ------------
        checkboxInput(
          inputId = ns("gene_centering"),
          label = "Center genes (substract mean)",
          value = TRUE
        ),
        checkboxInput(
          inputId = ns("gene_normalize"),
          label = "Normalize genes (divide by SD)",
          value = FALSE
        ),
        checkboxInput(
          inputId = ns("sample_centering"),
          label = "Center samples (substract mean)",
          value = FALSE
        ),
        checkboxInput(
          inputId = ns("sample_normalize"),
          label = "Normalize samples(divide by SD)",
          value = FALSE
        ),
        checkboxInput(
          inputId = ns("no_sample_clustering"),
          label = "Do not re-order or cluster samples",
          value = FALSE
        ),
        checkboxInput(
          inputId = ns("show_row_dend"),
          label = "Show Row Dendogram",
          value = TRUE
        ),

        # Sample coloring bar -----------
        htmlOutput(ns("list_factors_heatmap")),
        br(),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/heatmap/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          id = ns("heatmap_panels"),

          # Main heatmap ---------
          tabPanel(
            title = "Heatmap",
            br(),
            plotOutput(
              outputId = ns("heatmap_main"),
              width = "100%",
              height = "700px"
            )
          ),

          # Interactive heatmap panel ----------
          tabPanel(
            title = "Interactive Heatmap",
            br(),
            # INSERT PLOTLY OUTPUT
          ),

          # Gene Standard Deviation Distribution ----------
          tabPanel(
            title = "Gene SD Distribution",
            br(),
            plotOutput(
              outputId = ns("sd_density_plot"),
              width = "100%",
              height = "500px"
            )
          ),

          # Correlation matrix panel ----------
          tabPanel(
            title = "Correlation Matrix",
            br(),
            # Insert Correlation Matrix
          ),

          # Sample Tree Plot ---------
          tabPanel(
            title = "Sample Tree",
            br(),
            # INSERT SAMPLE TREE PLOT
          )
        )
      )
    )
  )
}

#' 03_heatmap Server Functions
#'
#' @noRd
mod_03_heatmap_server <- function(id, pre_process, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Update Slider Input ---------
    observe({
      req(tab() == "Heatmap")
      req(!is.null(pre_process$data()))
      if(nrow(pre_process$data()) > 12000) {
        max_genes <- 12000
      } else {
        max_genes <- round(nrow(pre_process$data()) + 50, -2)
      }
      updateSliderInput(
        inputId = "n_genes",
        value = c(0, 100),
        max = max_genes
      )
    })

    # Heatmap Colors ----------
    heatmap_colors <- list(
      "Green-Black-Red" = c("green", "black", "red"),
      "Blue-White-Red" = c("blue", "white", "red"),
      "Green-Black-Magenta" = c("green", "black", "magenta"),
      "Blue-Yellow-Red" = c("blue", "yellow", "red"),
      "Blue-White-Brown" = c("blue", "white", "brown")
    )
    heatmap_choices <- c(
      "Green-Black-Red",
      "Blue-White-Red",
      "Green-Black-Magenta",
      "Blue-Yellow-Red",
      "Blue-White-Brown"
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "heatmap_color_select",
        choices = heatmap_choices
      )
    })

    # Distance functions -----------
    dist_funs <- dist_functions()
    dist_choices <- setNames(
      1:length(dist_funs),
      names(dist_funs)
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "dist_function",
        choices = dist_choices
      )
    })

    # Hclust Functions ----------
    hclust_funs <- hcluster_functions()
    hclust_choices <- setNames(
      1:length(hclust_funs),
      names(hclust_funs)
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "hclust_functions",
        choices = hclust_choices
      )
    })

    # Gene ID Name Choices ----------
    observe({
      req(!is.null(pre_process$all_gene_names()))

      updateSelectInput(
        session = session,
        inputId = "select_gene_id",
        choices = colnames(pre_process$all_gene_names())
      )
    })

    # Sample color bar render ----------
    output$list_factors_heatmap <- renderUI({
      selectInput(
        inputId = ns("select_factors_heatmap"),
        label = "Sample color bar:",
        choices = c("Sample_Name", colnames(pre_process$sample_info()))
      )
    })

    # Standard Deviation Density Plot ----------
    output$sd_density_plot <- renderPlot({
      req(!is.null(pre_process$data()))

      sd_density(
        data = pre_process$data(),
        n_genes_max = input$n_genes[2],
        n_genes_min = input$n_genes[1]
      )
    })

    # Heatmap Data -----------
    heatmap_data <- reactive({
      req(!is.null(pre_process$data()))

      process_heatmap_data(
        data = pre_process$data(),
        n_genes_max = input$n_genes[2],
        n_genes_min = input$n_genes[1],
        gene_centering = input$gene_centering,
        gene_normalize = input$gene_normalize,
        sample_centering = input$sample_centering,
        sample_normalize = input$sample_normalize,
        all_gene_names = pre_process$all_gene_names(),
        select_gene_id = input$select_gene_id
      )
    })

    # Main heatmap ----------
    output$heatmap_main <- renderPlot({
      req(!is.null(heatmap_data()))
      req(!is.null(input$select_factors_heatmap))

      heatmap_main(
        data = heatmap_data(),
        n_genes = input$n_genes[2] - input$n_genes[1],
        heatmap_cutoff = input$heatmap_cutoff,
        sample_info = pre_process$sample_info(),
        select_factors_heatmap = input$select_factors_heatmap,
        dist_funs = dist_funs,
        dist_function = input$dist_function,
        hclust_function = input$hclust_function,
        no_sample_clustering = input$no_sample_clustering,
        heatmap_color_select = heatmap_colors[[input$heatmap_color_select]],
        row_dend = input$show_row_dend
      )
    })
  })
}

## To be copied in the UI
# mod_03_heatmap_ui("03_heatmap_ui_1")

## To be copied in the server
# mod_03_heatmap_server("03_heatmap_ui_1")
