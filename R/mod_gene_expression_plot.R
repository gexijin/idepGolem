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
      group_levels <- levels(expr_df$group)
      if (is.null(group_levels)) {
        group_levels <- unique(expr_df$group)
        expr_df$group <- factor(expr_df$group, levels = group_levels)
      }

      color_palette <- generate_colors(
        n = max(1, length(group_levels)),
        palette_name = palette_choice
      )

      has_raw <- "raw_data" %in% colnames(expr_df) &&
        any(!is.na(expr_df$raw_data))

      chosen_data <- input$data_type
      if (is.null(chosen_data) || (!has_raw && identical(chosen_data, "raw"))) {
        chosen_data <- "normalized"
      }
      use_raw <- identical(chosen_data, "raw") && has_raw

      expr_df$plot_value <- if (use_raw) {
        as.numeric(expr_df$raw_data)
      } else {
        as.numeric(expr_df$expression)
      }
      req(any(!is.na(expr_df$plot_value)))

      expression_label <- if (use_raw) {
        "Raw data"
      } else {
        "Normalized Expression"
      }

      plot_type <- input$plot_type
      if (is.null(plot_type) || !(plot_type %in% c("bar", "box", "violin"))) {
        plot_type <- default_plot_type
      }
      use_bar <- identical(plot_type, "bar")
      use_violin <- identical(plot_type, "violin")

      base_theme <- ggplot2::theme_light() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title.y = ggplot2::element_text(size = 14, color = "black"),
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(size = 12),
          axis.text.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 12),
          legend.title = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "grey60", fill = NA)
        )

      display_name <- payload$display_name
      if (is.null(display_name)) {
        display_name <- ""
      }

      if (use_bar) {
        expr_df <- expr_df[order(expr_df$group, expr_df$sample), , drop = FALSE]

        p <- ggplot2::ggplot(
          expr_df,
          ggplot2::aes(x = sample, y = plot_value, fill = group)
        ) +
          ggplot2::geom_col(color = "black", width = 0.7, na.rm = TRUE) +
          ggplot2::labs(
            title = display_name,
            y = expression_label
          ) +
          ggplot2::scale_fill_manual(values = color_palette) +
          base_theme +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
            legend.position = "right"
          )

        sample_count <- length(unique(expr_df$sample))
        if (sample_count < 50) {
          is_counts_upload <- isTRUE(counts_are_counts())
          expr_df$plot_label <- if (use_raw) {
            ifelse(
              is.na(expr_df$raw_data),
              NA_character_,
              if (is_counts_upload) {
                prettyNum(expr_df$raw_data, big.mark = ",", preserve.width = "none")
              } else {
                formatC(expr_df$raw_data, format = "f", digits = 3)
              }
            )
          } else {
            ifelse(
              is.na(expr_df$plot_value),
              NA_character_,
              formatC(expr_df$plot_value, format = "f", digits = 3)
            )
          }
          label_df <- expr_df[!is.na(expr_df$plot_label), , drop = FALSE]
          if (nrow(label_df) > 0) {
            max_y <- max(expr_df$plot_value, na.rm = TRUE)
            if (is.finite(max_y) && max_y > 0) {
              p <- p + ggplot2::expand_limits(y = max_y * 1.08)
            }
            p <- p +
              ggplot2::geom_text(
                data = label_df,
                mapping = ggplot2::aes(label = plot_label),
                vjust = -0.25,
                size = 3.2,
                color = "black"
              ) +
              ggplot2::coord_cartesian(clip = "off")
          }
        }
      } else {
        jitter_position <- ggplot2::position_jitter(width = 0.15, height = 0)

        p <- ggplot2::ggplot(
          expr_df,
          ggplot2::aes(x = group, y = plot_value, fill = group)
        ) +
          ggplot2::labs(
            title = display_name,
            y = expression_label
          ) +
          ggplot2::scale_fill_manual(values = color_palette) +
          base_theme +
          ggplot2::theme(legend.position = "none")

        if (use_violin) {
          p <- p +
            ggplot2::geom_violin(
              alpha = 0.7,
              color = "black",
              trim = FALSE,
              na.rm = TRUE
            )
        } else {
          p <- p +
            ggplot2::geom_boxplot(
              width = 0.65,
              alpha = 0.9,
              color = "black",
              outlier.shape = NA,
              na.rm = TRUE
            )
        }

        p <- p +
          ggplot2::geom_point(
            ggplot2::aes(color = group),
            position = jitter_position,
            size = 2.5,
            alpha = 0.85,
            na.rm = TRUE
          ) +
          ggplot2::scale_color_manual(values = color_palette) +
          ggplot2::guides(color = "none")
      }

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
