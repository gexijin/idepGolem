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
        downloadButton(
          outputId = ns("download_heatmap_data"),
          label = "Heatmap data"
        ),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/heatmap/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          id = ns("heatmap_panels"),

          # Heatmap panel ----------
          tabPanel(
            title = "Heatmap",
            h5("Brush for sub-heatmap, click for value. (Shown Below)"),
            br(),
            plotOutput(
              outputId = ns("heatmap_main"),
              height = "700px",
              width = "100%",
              brush = ns("ht_brush"),
              click = ns("ht_click")
            ),
            verbatimTextOutput(ns("ht_click_content")),
            plotOutput(
              outputId = ns("sub_heatmap"),
              height = "700px",
              width = "100%"
            ),
            h4("Sub-heatmap Data Table", align = "center"),
            DT::dataTableOutput(outputId = ns("subheat_data"))
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
            h5(
              "Using genes with maximum expression level at the top 75%.
               Data is transformed and clustered as specified in the sidebar."
            ),
            br(),
            plotOutput(
              outputId = ns("sample_tree"),
              width = "100%",
              height = "400px"
            )
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

    # Interactive heatmap environment
    shiny_env <- new.env()

    # Update Slider Input ---------
    observe({
      req(tab() == "Heatmap")
      req(!is.null(pre_process$data()))
      if (nrow(pre_process$data()) > 12000) {
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

    # HEATMAP -----------
    # Information on interactivity
    # https://jokergoo.github.io/2020/05/15/interactive-complexheatmap/
    output$heatmap_main <- renderPlot({
      req(!is.null(heatmap_data()))
      req(!is.null(input$select_factors_heatmap))

      # Assign heatmap to be used in multiple components
      shiny_env$ht <- heatmap_main(
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

      # Use heatmap position in multiple components
      shiny_env$ht_pos <- ComplexHeatmap::ht_pos_on_device(shiny_env$ht)
    })

    # Heatmap Click Value ---------
    output$ht_click_content <- renderText({

      if (is.null(input$ht_click)) {
        "Click for Info."
      } else {

        pos1 <- ComplexHeatmap:::get_pos_from_click(input$ht_click)

        ht <- shiny_env$ht
        pos <- ComplexHeatmap::selectPosition(
          ht,
          mark = FALSE,
          pos = pos1,
          verbose = FALSE,
          ht_pos = shiny_env$ht_pos
        )

        row_index <- pos[1, "row_index"]
        column_index <- pos[1, "column_index"]

        m <- ht@ht_list[[1]]@matrix
        value <- m[row_index, column_index]
        sample <- colnames(m)[column_index]
        gene <- rownames(m)[row_index]


        glue::glue(
          "Sample Name: {sample}",
          "Gene: {gene}",
          "Value: {value}",
          .sep = "\n"
        )
      }
    })

    # Subheatmap creation ---------
    output$sub_heatmap <- renderPlot({

      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("No region is selected.", 0.5, 0.5)
      } else {
        lt <- ComplexHeatmap:::get_pos_from_brush(input$ht_brush)
        pos1 <- lt[[1]]
        pos2 <- lt[[2]]

        ht <- shiny_env$ht
        pos <- ComplexHeatmap::selectArea(
          ht,
          mark = FALSE,
          pos1 = pos1,
          pos2 = pos2,
          verbose = FALSE,
          ht_pos = shiny_env$ht_pos
        )

        row_index <- unlist(pos[1, "row_index"])
        column_index <- unlist(pos[1, "column_index"])
        m <- ht@ht_list[[1]]@matrix
        if (length(row_index) > 50) {
          show_rows <- FALSE
        } else {
          show_rows <- TRUE
        }
        shiny_env$submap <- m[row_index, column_index, drop = FALSE]

        ht_select <- ComplexHeatmap::Heatmap(
          shiny_env$submap,
          col = ht@ht_list[[1]]@matrix_color_mapping@col_fun,
          show_heatmap_legend = FALSE,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = show_rows
        )
        ComplexHeatmap::draw(ht_select)
      }
    })

    # Subheatmap Data Table ----------
    output$subheat_data <- DT::renderDataTable({
      req(!is.null(input$ht_brush))

      DT::datatable(
        shiny_env$submap,
        options = list(
          pageLength = 10,
          scrollX = "400px"
        ),
        rownames = TRUE)
    })

    # Sample Tree ----------
    output$sample_tree <- renderPlot({
      req(!is.null(pre_process$data()))

      draw_sample_tree(
        tree_data = pre_process$data(),
        gene_centering = input$gene_centering,
        gene_normalize = input$gene_normalize,
        sample_centering = input$sample_centering,
        sample_normalize = input$sample_normalize,
        hclust_funs = hclust_funs,
        hclust_function = input$hclust_function,
        dist_funs = dist_funs,
        dist_function = input$dist_function
      )
    })

     # Heatmap Download Data -----------
    heatmap_data_download <- reactive({
      req(!is.null(pre_process$data()))

      down_data <- process_heatmap_data(
        data = pre_process$data(),
        n_genes_max = input$n_genes[2],
        n_genes_min = input$n_genes[1],
        gene_centering = input$gene_centering,
        gene_normalize = input$gene_normalize,
        sample_centering = input$sample_centering,
        sample_normalize = input$sample_normalize,
        all_gene_names = pre_process$all_gene_names(),
        select_gene_id = "User_ID"
      )

      merged_data <- merge_data(
        pre_process$all_gene_names(),
        down_data,
        merge_ID = "User_ID"
      )
    })

    output$download_heatmap_data <- downloadHandler(
      filename = function() {
        "heatmap_data.csv"
      },
      content = function(file) {
        write.csv(heatmap_data_download(), file)
      }
    )

  })
}

## To be copied in the UI
# mod_03_heatmap_ui("03_heatmap_ui_1")

## To be copied in the server
# mod_03_heatmap_server("03_heatmap_ui_1")
