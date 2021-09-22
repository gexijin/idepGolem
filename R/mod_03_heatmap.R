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

        # Select Clustering Method ----------
        selectInput(
          inputId = ns("cluster_meth"),
          label = "Select Clustering Method:",
          choices = list(
            "Hierarchical" = 1,
            "k-Means" = 2
          ),
          selected = 1
        ),

        # Gene ID Selection -----------
        selectInput(
          inputId = ns("select_gene_id"),
          label = "Select Gene ID Label (<= 50 genes):",
          choices = NULL,
          selected = NULL
        ),

        # Sample coloring bar -----------
        htmlOutput(ns("list_factors_heatmap")),

        # Heatmap customizing features ----------
        strong("Customize heatmap (Default values work well):"),
        fluidRow(
          br(),
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

        # Clustering methods for hierarchical ----------
        conditionalPanel(
          condition = "input.cluster_meth == 1",
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
          ns = ns
        ),

        # k- means slidebar -----------
        conditionalPanel(
          condition = "input.cluster_meth == 2",
          sliderInput(
            inputId = ns("k_clusters"),
            label = "Number of Clusters",
            min   = 2,
            max   = 20,
            value = 4,
            step  = 1
          ),
          actionButton(
            inputId = ns("k_means_re_run"),
            label = "Re-Run"
          ),
          ns = ns
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
            fluidRow(
              column(
                width = 3,
                plotOutput(
                  outputId = ns("heatmap_main"),
                  height = "450px",
                  width = "100%",
                  brush = ns("ht_brush")
                ),
                br(),
                h5("Selected Cell (Submap):"),
                uiOutput(
                  outputId = ns("ht_click_content"),
                  placeholder = TRUE
                )
              ),
              column(
                width = 9,
                plotOutput(
                  outputId = ns("sub_heatmap"),
                  height = "650px",
                  width = "100%",
                  click = ns("ht_click")
                )
              )
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
            checkboxInput(
              ns("label_pcc"),
              label = "Label w/ Pearson's correlation coefficients",
              value = TRUE
            ),
            plotOutput(
              outputId = ns("correlationMatrix")
            )
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
          ),

          # K-means elbow plot ----------
          tabPanel(
            title = "k-Cluster Plot",
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
            plotOutput(
              outputId = ns("k_clusters"),
              width = "100%",
              height = "500px"
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
        label = "Sample Color Bar:",
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

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Creating Heatmap",
        color = "#000000"
      )

      # Assign heatmap to be used in multiple components
      shiny_env$ht <- heatmap_main(
        data = heatmap_data(),
        cluster_meth = input$cluster_meth,
        heatmap_cutoff = input$heatmap_cutoff,
        sample_info = pre_process$sample_info(),
        select_factors_heatmap = input$select_factors_heatmap,
        dist_funs = dist_funs,
        dist_function = input$dist_function,
        hclust_function = input$hclust_function,
        no_sample_clustering = input$no_sample_clustering,
        heatmap_color_select = heatmap_colors[[input$heatmap_color_select]],
        row_dend = input$show_row_dend,
        k_clusters = input$k_clusters,
        re_run = input$k_means_re_run
      )

      # Use heatmap position in multiple components
      shiny_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht)

      shinybusy::remove_modal_spinner()

      return(shiny_env$ht)
    })

    # Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      if (is.null(input$ht_click)) { 
        "Click for Info."
      } else {
        heat_click_info(
          click = input$ht_click,
          ht_sub = shiny_env$ht_sub,
          ht_sub_obj = shiny_env$ht_sub_obj,
          ht_pos_sub = shiny_env$ht_pos_sub,
          sub_groups = shiny_env$sub_groups,
          group_colors = shiny_env$group_colors,
          cluster_meth = input$cluster_meth,
          click_data = shiny_env$click_data
        )
      }
    })

    # Subheatmap creation ---------
    output$sub_heatmap <- renderPlot({

      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("No region is selected.", 0.5, 0.5)
      } else {
        submap_return <- heat_sub(
          ht_brush = input$ht_brush,
          ht = shiny_env$ht,
          ht_pos_main = shiny_env$ht_pos_main,
          heatmap_data = heatmap_data(),
          sample_info = pre_process$sample_info(),
          select_factors_heatmap = input$select_factors_heatmap,
          cluster_meth = input$cluster_meth
        )

        # Objects used in other components ----------
        shiny_env$ht_sub_obj <- submap_return$ht_select
        shiny_env$submap_data <- submap_return$submap_data
        shiny_env$sub_groups <- submap_return$sub_groups
        shiny_env$group_colors <- submap_return$group_colors
        shiny_env$click_data <- submap_return$click_data
        
        shiny_env$ht_sub <- ComplexHeatmap::draw(
          shiny_env$ht_sub_obj,
          annotation_legend_list = submap_return$lgd,
          annotation_legend_side = "top"
        )

        shiny_env$ht_pos_sub <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht_sub)

        return(shiny_env$ht_sub)
      }
    })

    # Subheatmap Data Table ----------
    output$subheat_data <- DT::renderDataTable({
      req(!is.null(input$ht_brush))

      DT::datatable(
        shiny_env$submap_data,
        options = list(
          pageLength = 10,
          scrollX = "400px"
        ),
        rownames = TRUE
      )
    })

    # Correlation Matrix ----------
    output$correlationMatrix <- renderPlot({
		  cor_plot(
        data = pre_process$data(),
        label_pcc = input$label_pcc
      )
    }, height = 600, width = 700)

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

    # k-Cluster elbow plot ----------
    output$k_clusters <- renderPlot({
      req(!is.null(heatmap_data()))

      k_means_elbow(
        heatmap_data = heatmap_data()
      )
    })

     # Heatmap Download Data -----------
    heatmap_data_download <- reactive({
      req(!is.null(pre_process$all_gene_names()))
      req(!is.null(heatmap_data()))

      merged_data <- merge_data(
        pre_process$all_gene_names(),
        heatmap_data(),
        merge_ID = input$select_gene_id
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
