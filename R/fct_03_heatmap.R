#' fct_03_heatmap.R This file holds all of the main data analysis functions
#' associated with third tab of the iDEP website.
#'
#'
#' @section fct_03_heatmap.R functions:
#'
#'
#'
#' @name fct_03_heatmap.R
NULL


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
  hmcols <- grDevices::colorRampPalette(rev(c(
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
sd_density <- function(data,
                       n_genes_max,
                       n_genes_min) {
  sds <- apply(data[, 1:dim(data)[2]], 1, sd)
  max_sd <- mean(sds) + 4 * sd(sds)
  sds[sds > max_sd] <- max_sd

  if (n_genes_max - n_genes_min == 0) {
    n_genes_max <- n_genes_min + 10
  }

  if (n_genes_max > length(sds)) {
    n_genes_max <- length(sds)
  }

  if (n_genes_min > length(sds)) {
    n_genes_min <- length(sds) - 10
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

  if (n_genes_max == 10 && n_genes_min == 0) {
    return(plot)
  } else if (n_genes_min == 0) {
    cutoff <- sort(sds$sds, decreasing = TRUE)[n_genes_max]

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
        label = paste0("Upper: ", n_genes_max)
      )

    return(plot)
  } else {
    cutoff_max <- sort(sds$sds, decreasing = TRUE)[n_genes_max]
    cutoff_min <- sort(sds$sds, decreasing = TRUE)[n_genes_min]

    plot <- plot +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = cutoff_max),
        color = "red",
        linetype = "dashed",
        size = 1
      ) +
      ggplot2::annotate(
        "text",
        x = cutoff_max - 0.4 * sd(sds[, 1]),
        y = 1,
        colour = "red",
        label = paste0("Upper: ", n_genes_max)
      ) +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = cutoff_min),
        color = "red",
        linetype = "dashed",
        size = 1
      ) +
      ggplot2::annotate(
        "text",
        x = cutoff_min + 0.4 * sd(sds[, 1]),
        y = 1,
        colour = "red",
        label = paste0("Lower: ", n_genes_min)
      )

    return(plot)
  }
}


#' Heatmap data process
#'
#' This function prepares the data from pre-processing
#' to be displayed in a heatmap. It takes in limits for
#' what genes to subset, what centering and standardizing
#' to perform, and what gene ID label to use.
#'
#' @param data Processed data matrix
#' @param n_genes_max Row number upper limit to display in heatmap
#' @param n_genes_min Row number lower limit to display in heatmap
#' @param gene_centering TRUE/FALSE subtract mean from gene rows
#' @param gene_normalize TRUE/FALSE divide by SD in gene rows
#' @param sample_centering TRUE/FALSE subtract mean from sample columns
#' @param sample_normalize TRUE/FALSE divide by SD in sample columns
#' @param all_gene_names Data frame of gene names
#' @param select_gene_id Desired ID type for heatmap labels
#'   (User_ID, ensembl_ID, symbol)
#'
#' @return Subsetted data matrix ([n_genes_min:n_genes_max, ]) with
#'   gene IDs as the select_gene_id
process_heatmap_data <- function(data,
                                 n_genes_max,
                                 n_genes_min,
                                 gene_centering,
                                 gene_normalize,
                                 sample_centering,
                                 sample_normalize,
                                 all_gene_names,
                                 select_gene_id) {
  data <- rowname_id_swap(
    data_matrix = data,
    all_gene_names = all_gene_names,
    select_gene_id = select_gene_id
  )

  data <- data[order(-apply(
    data[, 1:dim(data)[2]],
    1,
    sd
  )), ]

  if (n_genes_max - n_genes_min < 10) {
    n_genes_max <- n_genes_min + 10
  } else if (n_genes_max > nrow(data)) {
    n_genes_max <- nrow(data)
  } else if (n_genes_min > nrow(data)) {
    n_genes_min <- nrow(data) - 10
  }

  # Center by genes, cutoff really large values ---------
  if (gene_centering) {
    data <-
      data[(n_genes_min + 1):n_genes_max, ] -
      apply(data[(n_genes_min + 1):n_genes_max, ], 1, mean)
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

  if (gene_centering) {
    return(data)
  } else {
    data <- data[(n_genes_min + 1):n_genes_max, ]
  }

  return(data)
}

#' Draw a heatmap of processed data
#'
#' Uses the package ComplexHeatmaps to draw a heatmap of the
#' processed data that has been prepared for the heatmap. The
#' returned heatmap visualizes the processed expression of the
#' gene range provided in process_heatmap_data.
#'
#' @param data Data returned from process_heatmap_data
#' @param n_genes Amount of genes included in the heatmap
#' @param heatmap_cutoff Z score max to filter data
#' @param sample_info Experiment design information from load data
#' @param select_factors_heatmap Factor group for annotation legend
#' @param dist_funs List of distance functions to use in heatmap
#' @param dist_function The selected distance function to use
#' @param hclust_function Type of clustering to perform
#' @param no_sample_clustering TRUE/FALSE Specify whehter to cluster columns
#' @param heatmap_color_select Vector of colors for heatmap scale
#' @param row_dend TRUE/FALSE Hide row dendogram
#'
#' @return Heatmap of the processed data.
heatmap_main <- function(data,
                         n_genes,
                         heatmap_cutoff,
                         sample_info,
                         select_factors_heatmap,
                         dist_funs,
                         dist_function,
                         hclust_function,
                         no_sample_clustering,
                         heatmap_color_select,
                         row_dend) {
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

  if (n_genes > 60) {
    show_gene_names <- FALSE
  } else {
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
    cluster_rows = TRUE,
    cluster_columns = !(no_sample_clustering),
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
    column_names_gp = grid::gpar(
      fontsize = cex_factor
    ),
    column_names_rot = 90
  )

  ComplexHeatmap::draw(
    heat,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "top"
  )
}
