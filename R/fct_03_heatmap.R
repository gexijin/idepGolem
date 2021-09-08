#' fct_03_heatmap.R This file holds all of the main data analysis functions
#' associated with third tab of the iDEP website.
#'
#'
#' @section fct_03_heatmap.R functions:
#' \code{add_legend}
#'
#'
#' @name fct_03_heatmap.R
NULL


# adding sample legends to heatmap; this is for the main heatmap
# https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param ... DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
add_legend <- function(...) {
  opar <- par(
    fig = c(0, 1, 0, 1),
    oma = c(0, 0, 0, 0),
    mar = c(0, 0, 0, 6),
    new = TRUE
  )
  on.exit(par(opar))
  plot(x = 0, y = 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend(...)
}

#' Get color choice list for the heatmap
#'
#' get_heatmap_colors
#'
#' Calling the function will provide a list of heatmap colors
#' for the select input bar to use and gives the colors that the
#' heatmap will use for the selected color set.
#'
#' @return Returns a list containing color choices and a data frame
#' of hex color values.
get_heatmap_colors <- function() {
  hmcols <- colorRampPalette(rev(c(
    "#D73027", "#FC8D59", "#FEE090", "#FFFFBF",
    "#E0F3F8", "#91BFDB", "#4575B4"
  )))(75)
  heat_colors <- rbind(
    gplots::greenred(75), gplots::bluered(75),
    gplots::colorpanel(75, "green", "black", "magenta"),
    gplots::colorpanel(75, "blue", "yellow", "red"), hmcols
  )
  rownames(heat_colors) <- c(
    "Green-Black-Red", "Blue-White-Red", "Green-Black-Magenta",
    "Blue-Yellow-Red", "Blue-white-brown"
  )
  color_choices <- setNames(1:dim(heat_colors)[1], rownames(heat_colors))

  return(
    list(
      color_choices = color_choices,
      heat_colors = heat_colors
    )
  )
}

#' SD DISTRIBUTION PLOT
sd_density <- function(
  data,
  top
) {

  sds <- apply(data[, 1:dim(data)[2]], 1, sd)
  max_sd <- mean(sds) + 4 * sd(sds)
  sds[sds > max_sd] <- max_sd

  if (top > length(sds)) {
    top <- length(sds)
  }

  sds <- as.data.frame(sds)

  plot <- ggplot2::ggplot(sds, ggplot2::aes(x = sds)) +
    ggplot2::geom_density(color = "darkblue", fill = "lightblue") +
    ggplot2::labs(
      title = "Standard Deviations of All Genes",
      y = "Density",
      x = "Standard Deviation"
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      ),
      axis.text.x = ggplot2::element_text(
        size = 14
      ),
      axis.text.y = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      legend.text = ggplot2::element_text(size = 12)
    )

  if (top == 0) {
    return(plot)
  }

  cutoff <- sort(sds$sds, decreasing = TRUE)[top]
  plot <- plot +
  ggplot2::geom_vline(
      ggplot2::aes(xintercept = cutoff),
      color = "red",
      linetype = "dashed",
      size = 1
  ) +
  ggplot2::annotate(
    "text",
    x = cutoff + 0.4 * sd(sds[, 1]),
    y = 1,
    colour = "red",
    label = paste0("Top ", top)
  )

  return(plot)
}


#' DATA FOR HEATMAP PANEL
process_heatmap_data <- function(
  data,
  n_genes,
  gene_centering,
  gene_normalize,
  sample_centering,
  sample_normalize
) {
  if (n_genes < 10) {
    n_genes <- 10
  } else if (n_genes > nrow(data)) {
    n_genes <- nrow(data)
  }

  # Center by genes, cutoff really large values ---------
  if (gene_centering) {
    data <-
      data[1:n_genes, ] - apply(data[1:n_genes, ], 1, mean)
  }

  # Standardize by gene ---------
  if (gene_normalize) {
    data <- data / apply(data, 1, sd)
  }

  # Row centering and normalize ----------
  data <- scale(
    data,
    center = sample_centering,
    scale = sample_normalize
  )

  data <- data[1:n_genes, ]

  return(data)
}

#' HEATMAP FUNCTION
heatmap_main <- function(
  data,
  n_genes,
  heatmap_cutoff,
  sample_info,
  select_factors_heatmap,
  dist_funs,
  dist_function,
  hclust_function,
  no_sample_clustering,
  heatmap_color_select
) {

  if (ncol(data) < 20) {
    cex_factor <- 12
  } else if (ncol(x) < 31) {
    cex_factor <- 10
  } else {
    cex_factor <- 8
  }

  cutoff <- median(unlist(data)) + heatmap_cutoff * sd(unlist(data))
  data[data > cutoff] <- cutoff
  cutoff <- median(unlist(data)) - heatmap_cutoff * sd(unlist(data))
  data[data < cutoff] <- cutoff

  groups <- detect_groups(colnames(data))

  if (!is.null(sample_info) && !is.null(select_factors_heatmap)) {
    if (select_factors_heatmap == "Sample_Name") {
      groups <- detect_groups(colnames(data))
    } else {
      ix <- match(select_factors_heatmap, colnames(sample_info))
      groups <- sample_info[, ix]
    }
  }

  groups_colors <- gg_color_hue(length(unique(groups)))

  if (n_genes > 1000) {
    row_dend <- FALSE
    show_gene_names <- FALSE
  } else if (n_genes > 60) {
    row_dend <- TRUE
    show_gene_names <- FALSE
  } else {
    row_dend <- TRUE
    show_gene_names <- TRUE
  }

  if (min(data) < 0) {
    col_fun <- circlize::colorRamp2(
      c(min(data), 0, max(data)),
      heatmap_color_select
    )
  } else {
    col_fun <- circlize::colorRamp2(
      c(min(data), median(data), max(data)),
      heatmap_color_select
    )
  }

  if (length(groups) < 30) {
    show_group_leg <- TRUE
  } else {
    show_group_leg <- FALSE
  }


  heat_ann <- ComplexHeatmap::HeatmapAnnotation(
    Group = groups,
    col = list(Group = setNames(groups_colors, unique(groups))),
    annotation_legend_param = list(
      Group = list(nrow = 1, title = NULL)
    ),
    show_annotation_name = list(Group = FALSE),
    show_legend = show_group_leg
  )

  if (no_sample_clustering) {
    heat <- ComplexHeatmap::Heatmap(
      data,
      name = "Expression",
      col = col_fun,
      row_order = rownames(data),
      column_order = colnames(data),
      show_column_dend = TRUE,
      show_row_dend = row_dend,
      row_dend_side = "left",
      row_dend_width = grid::unit(2, "cm"),
      top_annotation = heat_ann,
      show_row_names = show_gene_names,
      heatmap_legend_param = list(
        direction = "horizontal",
        legend_width = grid::unit(6, "cm"),
        title = "Color Key",
        title_position = "topcenter"
      ),
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = cex_factor)
    )
  } else {
    heat <- ComplexHeatmap::Heatmap(
      data,
      name = "Expression",
      col = col_fun,
      clustering_method_rows = hclust_function,
      clustering_method_columns = hclust_function,
      clustering_distance_rows = function(x) {
        dist_funs[[as.numeric(dist_function)]](x)
      },
      clustering_distance_columns = function(x) {
        dist_funs[[as.numeric(dist_function)]](x)
      },
      show_column_dend = TRUE,
      show_row_dend = row_dend,
      row_dend_side = "left",
      row_dend_width = grid::unit(2, "cm"),
      top_annotation = heat_ann,
      show_row_names = show_gene_names,
      heatmap_legend_param = list(
        direction = "horizontal",
        legend_width = grid::unit(6, "cm"),
        title = "Color Key",
        title_position = "topcenter"
      ),
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = cex_factor)
    )
  }

  ComplexHeatmap::draw(
    heat,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "top"
  )
}
