#' Build a gene expression plot shared across the app and reports.
#'
#' @param expr_df A data.frame with at least the columns `sample`, `group`,
#'   `expression`, and optionally `raw_data`.
#' @param display_name Character title shown on the plot.
#' @param palette_name Name of the palette to pass to `generate_colors()`.
#' @param plot_type One of `"bar"`, `"box"`, or `"violin"`.
#' @param data_type Either `"raw"` or `"normalized"`. Falls back to `"normalized"`
#'   if raw data are unavailable.
#' @param counts_are_counts Logical used to format integer labels when raw
#'   counts are displayed.
#' @param label_threshold Maximum number of samples to label on bar plots.
#'
#' @return A `ggplot` object.
#' @export
build_gene_expression_plot <- function(expr_df,
                                       display_name = "",
                                       palette_name = "Set1",
                                       plot_type = c("bar", "box", "violin"),
                                       data_type = c("raw", "normalized"),
                                       counts_are_counts = FALSE,
                                       label_threshold = 50) {
  stopifnot(is.data.frame(expr_df))
  if (!all(c("sample", "group", "expression") %in% colnames(expr_df))) {
    stop("`expr_df` must contain sample, group, and expression columns.", call. = FALSE)
  }

  plot_type <- match.arg(plot_type)
  data_type <- match.arg(data_type)

  df <- expr_df
  df$group <- droplevels(df$group)

  palette_choice <- palette_name
  if (is.null(palette_choice) || length(palette_choice) == 0 || !nzchar(palette_choice)) {
    palette_choice <- "Set1"
  }
  group_levels <- levels(df$group)
  if (is.null(group_levels)) {
    group_levels <- unique(df$group)
    df$group <- factor(df$group, levels = group_levels)
  }
  color_palette <- generate_colors(
    n = max(1, length(group_levels)),
    palette_name = palette_choice
  )

  has_raw <- "raw_data" %in% colnames(df) &&
    any(!is.na(df$raw_data))
  use_raw <- identical(data_type, "raw") && has_raw
  df$plot_value <- if (use_raw) {
    as.numeric(df$raw_data)
  } else {
    as.numeric(df$expression)
  }

  if (!any(!is.na(df$plot_value))) {
    stop("No numeric values available for plotting.", call. = FALSE)
  }

  expression_label <- if (use_raw) {
    "Raw data"
  } else {
    "Normalized Expression"
  }

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

  if (identical(plot_type, "bar")) {
    df <- df[order(df$group, df$sample), , drop = FALSE]
    p <- ggplot2::ggplot(
      df,
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

    sample_count <- length(unique(df$sample))
    if (sample_count < label_threshold) {
      is_counts_upload <- isTRUE(counts_are_counts)
      df$plot_label <- if (use_raw) {
        ifelse(
          is.na(df$raw_data),
          NA_character_,
          if (is_counts_upload) {
            prettyNum(df$raw_data, big.mark = ",", preserve.width = "none")
          } else {
            formatC(df$raw_data, format = "f", digits = 3)
          }
        )
      } else {
        ifelse(
          is.na(df$plot_value),
          NA_character_,
          formatC(df$plot_value, format = "f", digits = 3)
        )
      }

      label_df <- df[!is.na(df$plot_label), , drop = FALSE]
      if (nrow(label_df) > 0) {
        max_y <- max(df$plot_value, na.rm = TRUE)
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
      df,
      ggplot2::aes(x = group, y = plot_value, fill = group)
    ) +
      ggplot2::labs(
        title = display_name,
        y = expression_label
      ) +
      ggplot2::scale_fill_manual(values = color_palette) +
      base_theme +
      ggplot2::theme(legend.position = "none")

    if (identical(plot_type, "violin")) {
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

  p
}

default_marker_gene_definitions <- function() {
  list(
    list(symbol = "GAPDH", description = "Housekeeping gene"),
    list(symbol = "ACTB", description = "Housekeeping gene"),
    list(symbol = "H2AC6", description = "Histone mRNAs lack poly(A) tails"),
    list(symbol = "MT-CO1", description = "Mitochondrial mRNA"),
    list(symbol = "MT-RNR2", description = "16s mitochondrial rRNA"),
    list(symbol = "UTY", description = "Male-specific"),
    list(symbol = "XIST", description = "Female-specific")
  )
}

#' Prepare marker gene plot payloads for Shiny modules and reports.
#'
#' @param expr_matrix Transformed expression matrix or SummarizedExperiment.
#' @param sample_info Data frame describing the samples.
#' @param raw_counts Raw counts matrix, if available.
#' @param converted_data Converted counts used as a fallback for raw data.
#' @param all_gene_info Data frame with at least `symbol`, `ensembl_gene_id`,
#'   and optionally `description`.
#' @param marker_definitions Optional custom marker definitions.
#'
#' @return A list of payloads; each payload contains `symbol`, `description`,
#'   `data`, and `display_name`.
#' @export
marker_gene_plot_payloads <- function(expr_matrix,
                                      sample_info = NULL,
                                      raw_counts = NULL,
                                      converted_data = NULL,
                                      all_gene_info = NULL,
                                      marker_definitions = default_marker_gene_definitions()) {
  expr_mat <- normalize_expression_input(expr_matrix)
  if (is.null(expr_mat) || nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    return(list())
  }

  payloads <- lapply(marker_definitions, function(def) {
    payload <- marker_gene_payload_for_symbol(
      target_symbol = def$symbol,
      expr_matrix = expr_mat,
      sample_info = sample_info,
      raw_counts = raw_counts,
      converted_data = converted_data,
      all_gene_info = all_gene_info
    )
    if (is.null(payload)) {
      return(NULL)
    }
    payload$symbol <- def$symbol
    payload$description <- def$description
    payload
  })
  purrr::compact(payloads)
}

marker_gene_payload_for_symbol <- function(target_symbol,
                                           expr_matrix,
                                           sample_info,
                                           raw_counts,
                                           converted_data,
                                           all_gene_info) {
  target_symbol_clean <- trimws(target_symbol)
  if (!nzchar(target_symbol_clean)) {
    return(NULL)
  }
  target_symbol_upper <- toupper(target_symbol_clean)

  sample_names <- colnames(expr_matrix)
  if (length(sample_names) == 0) {
    return(NULL)
  }

  gene_info <- all_gene_info
  symbol_matches <- NULL
  if (!is.null(gene_info) &&
    "symbol" %in% colnames(gene_info)) {
    gene_symbols <- trimws(as.character(gene_info$symbol))
    symbol_matches <- gene_info[
      !is.na(gene_symbols) &
        toupper(gene_symbols) == target_symbol_upper,
      ,
      drop = FALSE
    ]
  }

  candidate_gene_ids <- character(0)
  symbol_label <- target_symbol_clean
  gene_description <- NULL

  if (!is.null(symbol_matches) && nrow(symbol_matches) > 0) {
    if ("ensembl_gene_id" %in% colnames(symbol_matches)) {
      candidate_gene_ids <- unique(na.omit(symbol_matches$ensembl_gene_id))
      candidate_gene_ids <- candidate_gene_ids[
        !is.na(candidate_gene_ids) &
          nzchar(candidate_gene_ids)
      ]
    }
    first_symbol <- trimws(as.character(symbol_matches$symbol[1]))
    if (!is.na(first_symbol) && nzchar(first_symbol)) {
      symbol_label <- first_symbol
    }
    if ("description" %in% colnames(symbol_matches)) {
      gene_description <- symbol_matches$description[1]
    }
  }

  rownames_expr <- rownames(expr_matrix)
  available_rows <- intersect(candidate_gene_ids, rownames_expr)
  if (length(available_rows) == 0) {
    cleaned_row_names <- trimws(toupper(rownames_expr))
    available_rows <- rownames_expr[cleaned_row_names == target_symbol_upper]
  }
  if (length(available_rows) == 0) {
    return(NULL)
  }
  gene_row <- available_rows[1]

  expr_values <- as.numeric(expr_matrix[gene_row, sample_names, drop = TRUE])
  if (all(is.na(expr_values))) {
    return(NULL)
  }

  sample_groups <- detect_groups(sample_names, sample_info)
  empty_group <- is.na(sample_groups) | sample_groups == ""
  sample_groups[empty_group] <- sample_names[empty_group]

  raw_matrix <- normalize_expression_input(raw_counts)
  if (is.null(raw_matrix)) {
    raw_matrix <- normalize_expression_input(converted_data)
  }

  raw_values <- rep(NA_real_, length(sample_names))
  if (!is.null(raw_matrix) && is.matrix(raw_matrix)) {
    available_samples <- intersect(sample_names, colnames(raw_matrix))
    if (length(available_samples) > 0) {
      candidate_ids <- unique(c(
        gene_row,
        candidate_gene_ids,
        target_symbol_upper,
        target_symbol_clean
      ))
      candidate_ids <- candidate_ids[
        !is.na(candidate_ids) &
          nzchar(candidate_ids)
      ]
      for (candidate in candidate_ids) {
        if (candidate %in% rownames(raw_matrix)) {
          raw_vec <- as.numeric(raw_matrix[candidate, available_samples, drop = TRUE])
          idx <- match(available_samples, sample_names)
          raw_values[idx] <- raw_vec
          break
        }
      }
    }
  }

  df <- data.frame(
    sample = sample_names,
    expression = expr_values,
    group = factor(sample_groups, levels = unique(sample_groups)),
    raw_data = raw_values,
    stringsAsFactors = FALSE
  )

  has_values <- any(!is.na(df$expression)) || any(!is.na(df$raw_data))
  if (!has_values) {
    return(NULL)
  }

  symbol_label <- toupper(symbol_label)
  display_name <- symbol_label
  if (!is.null(gene_description) &&
    !is.na(gene_description) &&
    nzchar(gene_description)) {
    display_name <- paste0(symbol_label, ": ", gene_description)
  }

  list(
    data = df,
    display_name = display_name
  )
}

normalize_expression_input <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      return(NULL)
    }
    x <- SummarizedExperiment::assay(x)
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    return(NULL)
  }
  storage.mode(x) <- "numeric"
  x
}
