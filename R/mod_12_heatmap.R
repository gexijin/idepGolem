#' 12_heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_12_heatmap_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        width = 3,
        plotOutput(
          outputId = ns("main_heatmap"),
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
        uiOutput(
          outputId = ns("ht_click_content")
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
    )
  )
}

#' 12_heatmap Server Functions
#'
#' @noRd
mod_12_heatmap_server <- function(id,
                                  data,
                                  bar,
                                  all_gene_names,
                                  cluster_rows = FALSE,
                                  heatmap_color) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Interactive heatmap environment
    shiny_env <- new.env()

    output$main_heatmap <- renderPlot({
      req(!is.null(data()))
      withProgress(message = "Creating Heatmap", {
        incProgress(0.2)
        # Assign heatmap to be used in multiple components
        shiny_env$ht <- deg_heatmap(
          data = data(),
          bar = bar(),
          heatmap_color_select = heatmap_color(),
          cluster_rows = cluster_rows
        )

        # Use heatmap position in multiple components
        shiny_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht)
      })
      return(shiny_env$ht)
    })

    main_heatmap_object <- reactive({
      req(!is.null(data()))
      withProgress(message = "Creating Heatmap for export", {
        incProgress(0.2)
        deg_heatmap(
          data = data(),
          bar = bar(),
          heatmap_color_select = heatmap_color(),
          cluster_rows = cluster_rows
        )
      })
    })

    dl_heatmap_main <- ottoPlots::mod_download_figure_server(
      id = "dl_heatmap_main",
      filename = "heatmap_main",
      figure = reactive({
        main_heatmap_object()
      }),
      width = 5,
      height = 8,
      label = "Above"
    )


    output$sub_heatmap <- renderPlot({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("Select a region on the heatmap to zoom in.

        Selection can be adjusted from the sides.
        It can also be dragged around.", 0.5, 0.5)
      } else {
        withProgress(message = "Zooming in...", {
          incProgress(0.2)

          heat_return <- deg_heat_sub(
            ht_brush = input$ht_brush,
            ht = shiny_env$ht,
            ht_pos_main = shiny_env$ht_pos_main,
            heatmap_data = data(),
            bar = bar(),
            all_gene_names = all_gene_names()
          )
          incProgress(0.6)
          shiny_env$ht_select <- heat_return$ht_select
          shiny_env$submap_data <- heat_return$submap_data
          shiny_env$group_colors <- heat_return$group_colors
          shiny_env$column_groups <- heat_return$column_groups
          shiny_env$bar <- heat_return$bar

          shiny_env$ht_sub <- ComplexHeatmap::draw(
            shiny_env$ht_select,
            annotation_legend_side = "top",
            heatmap_legend_side = "top"
          )

          shiny_env$ht_pos_sub <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht_sub)

          return(shiny_env$ht_sub)
        })
      }
    })

    sub_heatmap_object <- reactive({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("Select a region on the heatmap to zoom in.", 0.5, 0.5)
      } else {
        withProgress(message = "Creating sub-heatmap", {
          incProgress(0.2)
          heat_return <- deg_heat_sub(
            ht_brush = input$ht_brush,
            ht = shiny_env$ht,
            ht_pos_main = shiny_env$ht_pos_main,
            heatmap_data = data(),
            bar = bar(),
            all_gene_names = all_gene_names()
          )
          incProgress(0.6)
          shiny_env$ht_select <- heat_return$ht_select
          shiny_env$submap_data <- heat_return$submap_data
          shiny_env$group_colors <- heat_return$group_colors
          shiny_env$column_groups <- heat_return$column_groups
          shiny_env$bar <- heat_return$bar

          return(
            ComplexHeatmap::draw(
              shiny_env$ht_select,
              annotation_legend_side = "top",
              heatmap_legend_side = "top"
            )
          )
        })
      }
    })

    dl_heatmap_sub <- ottoPlots::mod_download_figure_server(
      id = "dl_heatmap_sub",
      filename = "heatmap_zoom",
      figure = reactive({
        sub_heatmap_object()
      }),
      width = 5,
      height = 7,
      label = "Right"
    )
    # Sub Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      # zoomed in, but not clicked
      if (is.null(input$ht_click) &&
        !is.null(shiny_env$ht_sub) &&
        !is.null(input$ht_brush)
      ) {
        p <- '<br><p style="color:red;text-align:right;">Click on the sub-heatmap &#10230;</p>'
        html <- GetoptLong::qq(p)
        return(HTML(html))
      }

      if (is.null(input$ht_click) ||
        is.null(shiny_env$ht_sub) ||
        is.null(input$ht_brush)
      ) {
        return(NULL)
      }

      deg_click_info(
        click = input$ht_click,
        ht_sub = shiny_env$ht_sub,
        ht_sub_obj = shiny_env$ht_select,
        ht_pos_sub = shiny_env$ht_pos_sub,
        sub_groups = shiny_env$column_groups,
        group_colors = shiny_env$group_colors,
        bar = shiny_env$bar,
        data = shiny_env$submap_data
      )
    })
  })
}

## To be copied in the UI
# mod_12_heatmap_ui("12_heatmap_1")

## To be copied in the server
# mod_12_heatmap_server("12_heatmap_1")
