#' Gene Expression Plot Module UI
#'
#' Renders a reusable plot with controls identical to the DEG tab gene popup.
#'
#' @param id Module ID.
#' @param plot_height Plot height passed to `plotOutput`.
#' @param show_download If `TRUE`, include download controls.
#'
#' @noRd
mod_gene_expression_plot_ui <- function(id,
                                        plot_height = "500px",
                                        show_download = TRUE) {
  ns <- NS(id)

  download_col <- NULL
  data_col_width <- 6
  plot_col_width <- 6

  if (isTRUE(show_download)) {
    download_col <- column(
      width = 2,
      ottoPlots::mod_download_figure_ui(
        id = ns("download_plot")
      )
    )
    data_col_width <- 5
    plot_col_width <- 5
  }

  tagList(
    plotOutput(
      outputId = ns("gene_expression_plot"),
      width = "100%",
      height = plot_height
    ),
    fluidRow(
      download_col,
      column(
        width = data_col_width,
        selectInput(
          inputId = ns("data_type"),
          label = NULL,
          choices = c(
            "Raw data" = "raw",
            "Transformed expression" = "normalized"
          ),
          selected = "raw",
          selectize = FALSE
        )
      ),
      column(
        width = plot_col_width,
        selectInput(
          inputId = ns("plot_type"),
          label = NULL,
          choices = c(
            "Sample bar plot" = "bar",
            "Group boxplot" = "box",
            "Group violin" = "violin"
          ),
          selected = "bar",
          selectize = FALSE
        )
      )
    )
  )
}

#' Gene Expression Plot Module Server
#'
#' @param id Module ID.
#' @param plot_data Reactive that returns a list with `data` (data.frame) and `display_name`.
#' @param palette_name Reactive returning a palette identifier.
#' @param plot_grid_lines Reactive logical passed to `refine_ggplot2`.
#' @param ggplot2_theme Reactive ggplot2 theme name passed to `refine_ggplot2`.
#' @param counts_are_counts Reactive logical indicating if raw values represent counts.
#' @param download_filename Filename for the download handler.
#' @param default_plot_type Default plot type ("bar", "box", or "violin").
#' @param default_data_type Default data type ("raw" or "normalized").
#'
#' @noRd
mod_gene_expression_plot_server <- function(id,
                                            plot_data,
                                            palette_name,
                                            plot_grid_lines,
                                            ggplot2_theme,
                                            counts_are_counts = reactive(FALSE),
                                            download_filename = "gene_expression_plot",
                                            default_plot_type = "bar",
                                            default_data_type = "raw") {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Update available choices whenever data changes so the UI
    # mirrors whether raw counts exist.
    observeEvent(plot_data(), {
      payload <- plot_data()
      if (is.null(payload)) {
        return()
      }
      expr_df <- payload$data
      has_raw <- "raw_data" %in% colnames(expr_df) &&
        any(!is.na(expr_df$raw_data))

      data_choices <- if (has_raw) {
        c("Raw data" = "raw", "Transformed expression" = "normalized")
      } else {
        c("Transformed expression" = "normalized")
      }

      selected_data_type <- input$data_type
      if (is.null(selected_data_type) || !(selected_data_type %in% data_choices)) {
        selected_data_type <- if (has_raw && default_data_type == "raw") "raw" else "normalized"
      }

      updateSelectInput(
        session = session,
        inputId = "data_type",
        choices = data_choices,
        selected = selected_data_type
      )

      updateSelectInput(
        session = session,
        inputId = "plot_type",
        selected = default_plot_type
      )
    }, ignoreNULL = FALSE)

    plot_object <- reactive({
      payload <- plot_data()
      req(payload)

      expr_df <- payload$data
      req(nrow(expr_df) > 0)

      has_raw_col <- "raw_data" %in% colnames(expr_df)
      if (has_raw_col) {
        valid_rows <- !(is.na(expr_df$expression) & is.na(expr_df$raw_data))
      } else {
        valid_rows <- !is.na(expr_df$expression)
      }
      expr_df <- expr_df[valid_rows, , drop = FALSE]
      req(nrow(expr_df) > 0)

      expr_df$group <- droplevels(expr_df$group)
      palette_choice <- palette_name()
      if (is.null(palette_choice) || length(palette_choice) == 0) {
        palette_choice <- "Set1"
      }
      has_raw <- "raw_data" %in% colnames(expr_df) &&
        any(!is.na(expr_df$raw_data))

      chosen_data <- input$data_type
      if (is.null(chosen_data) || (!has_raw && identical(chosen_data, "raw"))) {
        chosen_data <- "normalized"
      }
      use_raw <- identical(chosen_data, "raw") && has_raw

      plot_type <- input$plot_type
      if (is.null(plot_type) || !(plot_type %in% c("bar", "box", "violin"))) {
        plot_type <- default_plot_type
      }

      display_name <- payload$display_name
      if (is.null(display_name)) {
        display_name <- ""
      }

      counts_flag <- isTRUE(counts_are_counts())
      data_type_to_use <- if (use_raw) "raw" else "normalized"
      p <- build_gene_expression_plot(
        expr_df = expr_df,
        display_name = display_name,
        palette_name = palette_choice,
        plot_type = plot_type,
        data_type = data_type_to_use,
        counts_are_counts = counts_flag
      )
      refine_ggplot2(
        p = p,
        gridline = plot_grid_lines(),
        ggplot2_theme = ggplot2_theme()
      )
    })

    output$gene_expression_plot <- renderPlot({
      req(plot_object())
      print(plot_object())
    })

    ottoPlots::mod_download_figure_server(
      id = "download_plot",
      filename = download_filename,
      figure = reactive({
        plot_object()
      }),
      label = ""
    )

    return(plot_object)
  })
}
