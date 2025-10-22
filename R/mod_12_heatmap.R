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
        width = 4,
        plotOutput(
          outputId = ns("main_heatmap"),
          height = "450px",
          width = "100%",
          click = ns("ht_main_click"),
          brush = ns("ht_brush")
        ),
        tippy::tippy_this(
          ns("main_heatmap"),
          "Tip: Drag over any region of the heatmap to zoom into that selection. Click on the main heatmap for sample details.",
          theme = "light"
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
            ),
            align = "right"
          )
        ),
        uiOutput(
          outputId = ns("ht_click_content")
        )
      ),
      column(
        width = 8,
        plotOutput(
          outputId = ns("sub_heatmap"),
          height = "650px",
          width = "100%"
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
                                  heatmap_color,
                                  select_gene_id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Interactive heatmap environment
    shiny_env <- new.env()

    output$main_heatmap <- renderPlot(
      {
        req(!is.null(data()))
        withProgress(message = "Creating Heatmap", {
          incProgress(0.2)
          # Assign heatmap to be used in multiple components
          shiny_env$ht <- main_heatmap_object()

          # Use heatmap position in multiple components
          if (!is.null(shiny_env$ht)) {
            tryCatch({
              shiny_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht)
            }, error = function(e) {
              shiny_env$ht_pos_main <- NULL
            })
          } else {
            shiny_env$ht_pos_main <- NULL
          }
        })
        return(shiny_env$ht)
      },
      width = 250
    )

    main_heatmap_object <- reactive({
      req(!is.null(data()))
      withProgress(message = "Creating Heatmap for export", {
        incProgress(0.2)
        deg_heatmap(
          df = data(),
          bar = bar(),
          heatmap_color_select = unlist(strsplit(heatmap_color(), "-")),
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
      label = NULL
    )

    main_heatmap_group_info <- reactive({
      req(!is.null(data()))

      heatmap_data <- data()
      column_groups <- detect_groups(colnames(heatmap_data))
      group_count <- length(unique(column_groups))
      group_palette <- gg_color_hue(2 + group_count)

      column_colors <- setNames(
        group_palette[seq_len(group_count)],
        unique(column_groups)
      )

      bar_labels <- NULL
      change_colors <- NULL
      if (!is.null(bar())) {
        bar_vec <- bar()
        if (!is.null(bar_vec)) {
          bar_labels <- bar_vec
          bar_labels[bar_labels == -1] <- "Negative"
          bar_labels[bar_labels == 1] <- "Positive"
          bar_labels <- as.character(bar_labels)
          names(bar_labels) <- rownames(heatmap_data)

          change_levels <- unique(bar_labels)
          if (length(change_levels) > 0) {
            change_colors <- group_palette[
              (group_count + 1):(group_count + length(change_levels))
            ]
            change_colors <- setNames(change_colors, change_levels)
          }
        }
      }

      list(
        column_groups = column_groups,
        group_colors = if (is.null(change_colors)) {
          column_colors
        } else {
          c(column_colors, change_colors)
        },
        bar_labels = bar_labels
      )
    })

    # depending on the number of genes selected
    # change the height of the sub heatmap
    height_sub_heatmap <- reactive({
      if (is.null(input$ht_brush)) {
        return(400)
      }

      if (is.null(shiny_env$ht) || is.null(shiny_env$ht_pos_main)) {
        return(400)
      }

      # Get the row ids of selected genes
      tryCatch({
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
        # convert to height, pixels
        height1 <- max(
          400, # minimum
          min(
            30000, # maximum
            12 * length(row_index)
          )
        )
        return(height1) # max width is 1000
      }, error = function(e) {
        return(400)
      })
    })


    output$sub_heatmap <- renderPlot(
      {
        if (is.null(input$ht_brush)) {
          grid::grid.newpage()
          return(invisible(NULL))
        }

        shinybusy::show_modal_spinner(
          spin = "orbit",
          text = "Creating Heatmap",
          color = "#000000"
        )
        on.exit(shinybusy::remove_modal_spinner(), add = TRUE)

        heat_return <- sub_heatmap_calc()
        if (is.null(heat_return)) {
          grid::grid.newpage()
          return(invisible(NULL))
        }

        shiny_env$submap_data <- heat_return$submap_data

        ComplexHeatmap::draw(
          heat_return$ht_select,
          annotation_legend_side = "top",
          heatmap_legend_side = "top"
        )
      },
      height = reactive(height_sub_heatmap())
    )

    sub_heatmap_calc <- reactive({
      req(!is.null(data()))
      if (is.null(shiny_env$ht) || is.null(shiny_env$ht_pos_main)) {
        return(NULL)
      }
      
      submap_return <- tryCatch({ # tolerates error; otherwise stuck with spinner
        deg_heat_sub(
          ht_brush = input$ht_brush,
          ht = shiny_env$ht,
          ht_pos_main = shiny_env$ht_pos_main,
          heatmap_data = data(),
          bar = bar(),
          all_gene_names = all_gene_names(),
          select_gene_id = select_gene_id()
        )},
        error = function(e) {e$message}
      )
      
      if ("character" %in% class(submap_return)){
        submap_return <- NULL
      }

      if (!is.null(dim(submap_return$ht_select))){
        if (nrow(submap_return$ht_select) == 0 || 
            ncol(submap_return$ht_select) == 0) {
          submap_return <- NULL
        }
      }

      if (!is.null(submap_return)) {
        shiny_env$submap_data <- submap_return$submap_data
      } else {
        shiny_env$submap_data <- NULL
      }

      return(submap_return)
    })

    sub_heatmap_object <- reactive({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("Select a region on the heatmap to zoom in.", 0.5, 0.5)
      } else {
        heat_return <- sub_heatmap_calc()
        if (is.null(heat_return)) {
          grid::grid.newpage()
          grid::grid.text("Select a region on the heatmap to zoom in.", 0.5, 0.5)
        } else {
          shiny_env$submap_data <- heat_return$submap_data
          ComplexHeatmap::draw(
            heat_return$ht_select,
            annotation_legend_side = "top",
            heatmap_legend_side = "top"
          )
        }
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
      label = "\u2192"
    )
    # Sub Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      main_ready <- !is.null(shiny_env$ht) &&
        !is.null(shiny_env$ht_pos_main) &&
        length(shiny_env$ht@ht_list) >= 1

      if (!is.null(input$ht_main_click) && main_ready) {
        info <- tryCatch(main_heatmap_group_info(), error = function(e) NULL)
        if (!is.null(info)) {
          matrix_obj <- shiny_env$ht@ht_list[[1]]
          if (!is.null(matrix_obj@matrix)) {
            click_data_main <- matrix_obj@matrix
            column_groups <- info$column_groups
            if (is.factor(column_groups)) {
              column_groups <- as.character(column_groups)
            }
            bar_labels <- info$bar_labels
            if (!is.null(bar_labels) && is.factor(bar_labels)) {
              bar_labels <- as.character(bar_labels)
            }
            if (!is.null(bar_labels)) {
              bar_labels <- as.character(bar_labels)
            }
            return(deg_click_info(
              click = input$ht_main_click,
              ht_sub = shiny_env$ht,
              ht_sub_obj = matrix_obj,
              ht_pos_sub = shiny_env$ht_pos_main,
              sub_groups = column_groups,
              group_colors = info$group_colors,
              bar = bar_labels,
              data = click_data_main
            ))
          }
        }
      }

      if (!is.null(input$ht_brush) && is.null(input$ht_main_click)) {
        note <- '<br><p style="color:red;text-align:right;">Click on the main heatmap &#10230;</p>'
        html <- GetoptLong::qq(note)
        return(HTML(html))
      }

      return(NULL)
    })
  })
}

## To be copied in the UI
# mod_12_heatmap_ui("12_heatmap_1")

## To be copied in the server
# mod_12_heatmap_server("12_heatmap_1")
