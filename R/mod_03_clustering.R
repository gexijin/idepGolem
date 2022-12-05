#' 03_heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_03_clustering_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "Clustering",
    sidebarLayout(

      # Heatmap Panel Sidebar ----------
      sidebarPanel(
        width = 3,
        # Select Clustering Method ----------
        conditionalPanel(
          condition = "input.cluster_panels == 'Heatmap' |
            input.cluster_panels == 'sample_tab'",
          radioButtons(
            inputId = ns("cluster_meth"),
            label = NULL,
            choices = list(
              "Hierarchical Clustering" = 1,
              "k-Means Clustering" = 2
            ),
            selected = 1
          ),
          ns = ns
        ),
        HTML('<hr style="height:1px;border:none;color:#333;background-color:#333;" />'),
        conditionalPanel(
          condition = "input.cluster_panels == 'Heatmap' |
          input.cluster_panels == 'Gene SD Distribution' ",
          fluidRow(
            column(width = 6, p("Top Genes:")),
            column(
              width = 6,
              numericInput(
                inputId = ns("n_genes"),
                label = NULL,
                min = 10,
                max = 12000,
                value = 1000,
                step = 10
              ),
              tippy::tippy_this(
                ns("n_genes"),
                "Genes are ranked by standard deviations across samples
                based on the transformed data.
                By showing the patterns of the genes with most varability,
                this give us a big picture view of the general pattern of
                gene expression.",
                theme = "light-border"
              )
            )
          ),
          ns = ns
        ),
        conditionalPanel(
          condition = "(input.cluster_panels == 'Heatmap' |
            input.cluster_panels == 'sample_tab') &&  input.cluster_meth == 2",

          # k- means slidebar -----------

          sliderInput(
            inputId = ns("k_clusters"),
            label = "Number of Clusters:",
            min = 2,
            max = 20,
            value = 6,
            step = 1
          ),

          # Re-run k-means with a different seed
          actionButton(
            inputId = ns("k_means_re_run"),
            label = "Re-Run"
          ),
          tippy::tippy_this(
            ns("k_means_re_run"),
            "Re-run the k-Means algorithm using different seeds for random number generator.",
            theme = "light-border"
          ),
          # Elbow plot pop-up
          actionButton(
            inputId = ns("elbow_pop_up"),
            label = "How many clusters?"
          ),
          tippy::tippy_this(
            ns("elbow_pop_up"),
            "k-Means elbow plot",
            theme = "light-border"
          ),
          # Line break ---------
          HTML(
            '<hr style="height:1px;border:none;
           color:#333;background-color:#333;" />'
          ),
          ns = ns
        ),

        # Clustering methods for Heatmap ----------
        conditionalPanel(
          condition = "input.cluster_meth == 1 &&
            (input.cluster_panels == 'Heatmap' |
            input.cluster_panels == 'sample_tab')",
          fluidRow(
            column(width = 4, p("Distance")),
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
            column(width = 4, p("Linkage")),
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
          ns = ns
        ),
        conditionalPanel(
          condition = "input.cluster_panels == 'Heatmap' ",
          fluidRow(
            column(width = 4, p("Samples color")),
            column(
              width = 8,
              htmlOutput(ns("list_factors_heatmap"))
            )
          ),
          fluidRow(
            column(width = 4, p("Label Genes:")),
            column(
              width = 8,
              htmlOutput(ns("selected_genes_ui"))
            )
          ),
          checkboxInput(ns("customize_button"), "More options"),
          checkboxInput(
            inputId = ns("sample_clustering"),
            label = "Cluster samples",
            value = FALSE
          ),
          checkboxInput(
            inputId = ns("show_row_dend"),
            label = "Show Row Dendogram",
            value = TRUE
          ),
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
          numericInput(
            inputId = ns("heatmap_cutoff"),
            label = "Max Z score:",
            value = 3,
            min = 2,
            step = 1
          ),
          fluidRow(
            column(
              width = 6,
              downloadButton(
                outputId = ns("download_heatmap_data"),
                label = "Heatmap data"
              )
            )
          ),
          ns = ns
        ),
        br(),
        downloadButton(
          outputId = ns("report"),
          label = "Report"
        ),
        tippy::tippy_this(
          ns("report"),
          "Generate HTML report of clustering tab",
          theme = "light-border"
        ),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/heatmap/",
          target = "_blank"
        )
      ),





      #########################################################################
      # Main Panel
      #########################################################################

      mainPanel(
        tabsetPanel(
          id = ns("cluster_panels"),
          # Heatmap panel ----------
          tabPanel(
            title = "Heatmap",
            br(),
            fluidRow(
              column(
                width = 4,
                plotOutput(
                  outputId = ns("heatmap_main"),
                  height = "450px",
                  width = "100%",
                  brush = ns("ht_brush")
                ),
                fluidRow(
                  column(
                    width = 6,
                    ottoPlots::mod_download_figure_ui(
                      ns("dl_heatmap_main")
                    )
                  ),
                  column(
                    width = 6,
                    ottoPlots::mod_download_figure_ui(
                      ns("dl_heatmap_sub")
                    )
                  )
                ),
                br(),
                uiOutput(
                  outputId = ns("ht_click_content")
                )
              ),
              column(
                width = 8,
                # align = "right",
                p("Broaden your browser window if there is overlap.-->"),
                checkboxInput(
                  inputId = ns("cluster_enrichment"),
                  label = strong("Show enrichment"),
                  value = FALSE
                ),
                tippy::tippy_this(
                  ns("cluster_enrichment"),
                  "Conducts GO enrichment on the selected genes.
                  For hierarchical clustering, users need to select a
                  region to zoom in first.
                  When k-means is used, enrichment analyses are
                  conducted on all clusters, regardless of your selection.",
                  theme = "light-border"
                ),
                conditionalPanel(
                  condition = "input.cluster_enrichment == 1 ",
                  mod_11_enrichment_ui(ns("enrichment_table_cluster")),
                  ns = ns
                ),
                plotOutput(
                  outputId = ns("sub_heatmap"),
                  height = "100%",
                  width = "100%",
                  click = ns("ht_click")
                )
              )
            )
          ),

          # Gene Standard Deviation Distribution ----------
          tabPanel(
            title = "Gene SD Distribution",
            br(),
            plotOutput(
              outputId = ns("sd_density_plot"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(ns("dl_gene_dist"))
          ),
          # Sample Tree -----------------
          tabPanel(
            title = "Sample Tree",
            value = "sample_tab",
            h5(
              "Using genes with maximum expression level at the top 75%.
               Data is transformed and clustered as specified in the sidebar."
            ),
            br(),
            plotOutput(
              outputId = ns("sample_tree"),
              width = "100%",
              height = "400px"
            ),
            ottoPlots::mod_download_figure_ui(ns("dl_sample_tree"))
          )
        )
      )
    )
  )
}








#########################################################################
# Server function
#########################################################################

#' 03_heatmap Server Functions
#'
#' @noRd
mod_03_clustering_server <- function(id, pre_process, load_data, idep_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Interactive heatmap environment
    shiny_env <- new.env()

    # Update Slider Input ---------
    observe({
      req(tab() == "Clustering")
      req(!is.null(pre_process$data()))
      if (nrow(pre_process$data()) > 12000) {
        max_genes <- 12000
      } else {
        max_genes <- round(nrow(pre_process$data()), -2)
      }
      updateNumericInput(
        inputId = "n_genes",
        max = max_genes
      )
    })

    observe({
      shinyjs::toggle(id = "sample_clustering", condition = input$customize_button)
      shinyjs::toggle(id = "show_row_dend", condition = input$customize_button)
      shinyjs::toggle(id = "heatmap_cutoff", condition = input$customize_button)
      shinyjs::toggle(id = "gene_normalize", condition = input$customize_button)
      shinyjs::toggle(id = "gene_centering", condition = input$customize_button)
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

    # Sample color bar selector ----------
    output$list_factors_heatmap <- renderUI({
      choices <- "Names"
      selected <- choices
      if (!is.null(colnames(pre_process$sample_info()))) {
        factors <- colnames(pre_process$sample_info())
        choices <- c(
          choices,
          factors,
          "All factors"
        )
        selected <- choices[length(choices)]
      }
      selectInput(
        inputId = ns("select_factors_heatmap"),
        label = NULL,
        choices = choices,
        selected = selected
      )
    })

    # Standard Deviation Density Plot ----------
    sd_density_plot <- reactive({
      req(!is.null(pre_process$data()))

      p <- sd_density(
        data = pre_process$data(),
        n_genes_max = input$n_genes
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })

    output$sd_density_plot <- renderPlot({
      print(sd_density_plot())
    })

    dl_gene_dist <- ottoPlots::mod_download_figure_server(
      id = "dl_gene_dist",
      filename = "sd_density_plot",
      figure = reactive({
        sd_density_plot()
      }),
      label = ""
    )



    # Heatmap Data -----------
    heatmap_data <- reactive({
      req(!is.null(pre_process$data()))

      process_heatmap_data(
        data = pre_process$data(),
        n_genes_max = input$n_genes,
        gene_centering = input$gene_centering,
        gene_normalize = input$gene_normalize,
        sample_centering = FALSE,
        sample_normalize = FALSE,
        all_gene_names = pre_process$all_gene_names(),
        select_gene_id = pre_process$select_gene_id()
      )
    })


    # Heatmap Click Value ---------
    output$selected_genes_ui <- renderUI({
      req(!is.null(pre_process$all_gene_names()))
      req(!is.null(pre_process$data()))
      req(!is.null(heatmap_data()))

      selectizeInput(
        inputId = ns("selected_genes"),
        label = NULL,
        choices = row.names(heatmap_data()),
        selected = NULL,
        multiple = TRUE
      )
    })

    # split "green-white-red" to c("green", "white", "red")
    heatmap_color_select <- reactive({
      req(pre_process$heatmap_color_select())
      unlist(strsplit(pre_process$heatmap_color_select(), "-"))
    })

    # HEATMAP -----------
    # Information on interactivity
    # https://jokergoo.github.io/2020/05/15/interactive-complexheatmap/
    output$heatmap_main <- renderPlot(
      {
        req(!is.null(heatmap_data()))
        req(input$select_factors_heatmap)
        #      req(input$selected_genes)

        shinybusy::show_modal_spinner(
          spin = "orbit",
          text = "Creating Heatmap",
          color = "#000000"
        )

        shiny_env$ht <- heatmap_main_object()

        # Use heatmap position in multiple components
        shiny_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht)
        shinybusy::remove_modal_spinner()
        return(shiny_env$ht)
      },
      width = 240 # , # this avoids the heatmap being redraw
      # height = 600
    )

    heatmap_main_object <- reactive({
      req(!is.null(heatmap_data()))
      req(input$select_factors_heatmap)

      # Assign heatmap to be used in multiple components
      try(
        obj <- heatmap_main(
          data = heatmap_data(),
          cluster_meth = input$cluster_meth,
          heatmap_cutoff = input$heatmap_cutoff,
          sample_info = pre_process$sample_info(),
          select_factors_heatmap = input$select_factors_heatmap,
          dist_funs = dist_funs,
          dist_function = input$dist_function,
          hclust_function = input$hclust_function,
          sample_clustering = input$sample_clustering,
          heatmap_color_select = heatmap_color_select(),
          row_dend = input$show_row_dend,
          k_clusters = input$k_clusters,
          re_run = input$k_means_re_run,
          selected_genes = input$selected_genes
        )
      )

      return(obj)
    })


    dl_heatmap_main <- ottoPlots::mod_download_figure_server(
      id = "dl_heatmap_main",
      filename = "heatmap_main",
      figure = reactive({
        heatmap_main_object()
      }),
      width = 6,
      height = 16,
      label = "Above"
    )


    # Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      # zoomed in, but not clicked
      if (is.null(input$ht_click) &&
        !is.null(shiny_env$ht_sub) &&
        !is.null(input$ht_brush)
      ) {
        p <- '<br><p style="color:red;text-align:left;">Click on the sub-heatmap for more info</p>'
        html <- GetoptLong::qq(p)
        return(HTML(html))
      }
      # if not zoomed in, show nothing
      if (is.null(input$ht_click) ||
        is.null(shiny_env$ht_sub) ||
        is.null(input$ht_brush)
      ) {
        return(NULL)
      }

      cluster_heat_click_info(
        click = input$ht_click,
        ht_sub = shiny_env$ht_sub,
        ht_sub_obj = shiny_env$ht_sub_obj,
        ht_pos_sub = shiny_env$ht_pos_sub,
        sub_groups = shiny_env$sub_groups,
        group_colors = shiny_env$group_colors,
        cluster_meth = input$cluster_meth,
        click_data = shiny_env$click_data
      )
    })

    # depending on the number of genes selected
    # change the height of the sub heatmap
    height_sub_heatmap <- reactive({
      if (is.null(input$ht_brush)) {
        return(400)
      }

      # Get the row ids of selected genes
      lt <- InteractiveComplexHeatmap::getPositionFromBrush(input$ht_brush)
      pos1 <- lt[[1]]
      pos2 <- lt[[2]]
      pos <- InteractiveComplexHeatmap::selectArea(
        shiny_env$ht,
        mark = FALSE,
        pos1 = pos1,
        pos2 = pos2,
        verbose = FALSE,
        ht_pos = shiny_env$ht_pos_main
      )
      row_index <- unlist(pos[1, "row_index"])
      # convert to height, pxiels
      height1 <- max(
        400, # minimum
        min(
          30000, # maximum
          12 * length(row_index)
        )
      )
      return(height1) # max width is 1000
    })

    # Subheatmap creation ---------
    output$sub_heatmap <- renderPlot(
      {
        if (is.null(input$ht_brush)) {
          grid::grid.newpage()
          grid::grid.text("Select a region on the heatmap to zoom in.
        Selection can be adjusted from the sides.
        It can also be dragged around.
        ", 0.5, 0.5)
        } else {
          shinybusy::show_modal_spinner(
            spin = "orbit",
            text = "Creating sub-heatmap",
            color = "#000000"
          )
          try(
            submap_return <- heatmap_sub_object_calc()
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

          shinybusy::remove_modal_spinner()
          return(shiny_env$ht_sub)
        }
      },
      # adjust height of the zoomed in heatmap dynamically based on selection
      height = reactive(height_sub_heatmap())
      # width = 500 # this avoids the heatmap being redraw
    )



    heatmap_sub_object_calc <- reactive({
      try( # tolerates error; otherwise stuck with spinner
        submap_return <- heat_sub(
          ht_brush = input$ht_brush,
          ht = shiny_env$ht,
          ht_pos_main = shiny_env$ht_pos_main,
          heatmap_data = heatmap_data(),
          sample_info = pre_process$sample_info(),
          select_factors_heatmap = input$select_factors_heatmap,
          cluster_meth = input$cluster_meth
        )
      )

      return(submap_return)
    })
    # Subheatmap creation ---------
    heatmap_sub_object <- reactive({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("Select a region on the heatmap to zoom in.", 0.5, 0.5)
      } else {
        try(
          submap_return <- heatmap_sub_object_calc()
        )

        # Objects used in other components ----------
        shiny_env$ht_sub_obj <- submap_return$ht_select
        shiny_env$submap_data <- submap_return$submap_data
        shiny_env$sub_groups <- submap_return$sub_groups
        shiny_env$group_colors <- submap_return$group_colors
        shiny_env$click_data <- submap_return$click_data

        return(ComplexHeatmap::draw(
          shiny_env$ht_sub_obj,
          annotation_legend_list = submap_return$lgd,
          annotation_legend_side = "top"
        ))
      }
    })

    dl_heatmap_sub <- ottoPlots::mod_download_figure_server(
      id = "dl_heatmap_sub",
      filename = "heatmap_zoom",
      figure = reactive({
        heatmap_sub_object()
      }),
      width = 8,
      height = 12,
      label = "Right"
    )

    # gene lists for enrichment analysis
    gene_lists <- reactive({
      req(!is.null(pre_process$select_gene_id()))
      req(!is.null(input$ht_brush) || input$cluster_meth == 2)

      gene_lists <- list()

      if (input$cluster_meth == 1) {
        gene_names <- merge_data(
          all_gene_names = pre_process$all_gene_names(),
          data = shiny_env$submap_data,
          merge_ID = pre_process$select_gene_id()
        )

        # Only keep the gene names and scrap the data
        gene_lists[["Selection"]] <- dplyr::select_if(gene_names, is.character)

        # k-means-----------------------------------------------------
      } else if (input$cluster_meth == 2) {
        # Get the cluster number and Gene

        req(heatmap_data())
        req(input$k_clusters)
        req(pre_process$select_gene_id())
        req(shiny_env$ht)

        row_ord <- ComplexHeatmap::row_order(shiny_env$ht)

        req(!is.null(names(row_ord)))

        for (i in 1:length(row_ord)) {
          if (i == 1) {
            clusts <- data.frame(
              "cluster" = rep(names(row_ord[i]), length(row_ord[[i]])),
              "row_order" = row_ord[[i]]
            )
          } else {
            tem <- data.frame(
              "cluster" = rep(names(row_ord[i]), length(row_ord[[i]])),
              "row_order" = row_ord[[i]]
            )
            clusts <- rbind(clusts, tem)
          }
        }
        clusts$id <- rownames(heatmap_data()[clusts$row_order, ])

        # disregard user selection use clusters for enrichment
        for (i in 1:input$k_clusters) {
          cluster_data <- subset(clusts, cluster == i)
          row.names(cluster_data) <- cluster_data$id

          gene_names <- merge_data(
            all_gene_names = pre_process$all_gene_names(),
            data = cluster_data,
            merge_ID = pre_process$select_gene_id()
          )

          # Only keep the gene names and scrap the data
          gene_lists[[paste0("", i)]] <-
            dplyr::select_if(gene_names, is.character)
        }
      }

      return(gene_lists)
    })

    # Sample Tree ----------
    sample_tree <- reactive({
      req(!is.null(pre_process$data()), input$cluster_meth == 1)

      draw_sample_tree(
        tree_data = pre_process$data(),
        gene_centering = input$gene_centering,
        gene_normalize = input$gene_normalize,
        sample_centering = FALSE,
        sample_normalize = FALSE,
        hclust_funs = hclust_funs,
        hclust_function = input$hclust_function,
        dist_funs = dist_funs,
        dist_function = input$dist_function
      )
      p <- recordPlot()
      return(p)
    })

    output$sample_tree <- renderPlot({
      print(sample_tree())
    })

    dl_sample_tree <- ottoPlots::mod_download_figure_server(
      id = "dl_sample_tree",
      filename = "sample_tree",
      figure = reactive({
        sample_tree()
      }),
      label = ""
    )

    observeEvent(input$cluster_meth, {
      if (input$cluster_meth == 1) {
        showTab(
          inputId = "cluster_panels",
          target = "sample_tab"
        )
      }
    })

    observeEvent(input$cluster_meth, {
      if (input$cluster_meth == 2) {
        hideTab(
          inputId = "cluster_panels",
          target = "sample_tab"
        )
      }
    })

    # k-Cluster elbow plot ----------
    output$k_clusters <- renderPlot({
      req(!is.null(heatmap_data()))

      k_means_elbow(
        heatmap_data = heatmap_data()
      )
    })
    # pop-up modal
    observeEvent(input$elbow_pop_up, {
      showModal(modalDialog(
        plotOutput(ns("k_clusters")),
        footer = NULL,
        easyClose = TRUE,
        title = tags$h5(
          "Following the elbow method, one should choose k so that adding
          another cluster does not substantially reduce the within groups sum of squares.",
          tags$a(
            "Wikipedia",
            href = "https://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set",
            target = "_blank"
          )
        ),
      ))
    })

    # Heatmap Download Data -----------
    heatmap_data_download <- reactive({
      req(!is.null(pre_process$all_gene_names()))
      req(!is.null(heatmap_data()))

      data <- prep_download(
        heatmap = heatmap_main_object(),
        heatmap_data = heatmap_data(),
        cluster_meth = input$cluster_meth
      )

      merged_data <- merge_data(
        pre_process$all_gene_names(),
        data,
        merge_ID = pre_process$select_gene_id()
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

    # Markdown report------------
    output$report <- downloadHandler(

      # For PDF output, change this to "report.pdf"
      filename = "clustering_report.html",
      content = function(file) {
        withProgress(message = "Generating report", {
          incProgress(0.2)
          # Copy the report file to a temporary directory before processing it, in
          # case we don't have write permissions to the current working dir (which
          # can happen when deployed).
          tempReport <- file.path(tempdir(), "clustering_workflow.Rmd")
          # tempReport
          tempReport <- gsub("\\", "/", tempReport, fixed = TRUE)

          # This should retrieve the project location on your device:
          # "C:/Users/bdere/Documents/GitHub/idepGolem"
          wd <- getwd()

          markdown_location <- paste0(wd, "/vignettes/Reports/clustering_workflow.Rmd")
          file.copy(from = markdown_location, to = tempReport, overwrite = TRUE)

          # Set up parameters to pass to Rmd document
          params <- list(
            pre_processed_data = pre_process$data(),
            sample_info = pre_process$sample_info(),
            descr = pre_process$descr(),
            all_gene_names = pre_process$all_gene_names(),
            n_genes = input$n_genes,
            k_clusters = input$k_clusters,
            cluster_meth = input$cluster_meth,
            select_gene_id = pre_process$select_gene_id(),
            list_factors_heatmap = input$list_factors_heatmap,
            heatmap_color_select = heatmap_color_select(),
            dist_function = input$dist_function,
            hclust_function = input$hclust_function,
            heatmap_cutoff = input$heatmap_cutoff,
            gene_centering = input$gene_centering,
            gene_normalize = input$gene_normalize,
            sample_clustering = input$sample_clustering,
            show_row_dend = input$show_row_dend,
            selected_genes = input$selected_genes
          )

          req(params)

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
# mod_03_heatmap_ui("03_heatmap_ui_1")

## To be copied in the server
# mod_03_heatmap_server("03_heatmap_ui_1")
