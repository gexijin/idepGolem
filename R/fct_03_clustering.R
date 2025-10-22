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


#' Density plot of data standard deviation
#'
#' Draw a density plot of the standard deviation in the
#' data. The function also adds a  vertical red lines for a specified range of
#' genes.
#'
#' @param data Data matrix that has been through pre-processing
#' @param n_genes_max Integer for the upper limit of gene range
#'
#' @export
#' @return Formatted density plot of the standard deviation
#' distribution.
#'
#' @family plots
#' @family clustering functions
sd_density <- function(data,
                       n_genes_max) {
  n_genes_min <- 5
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
        x = cutoff_max + 0.4 * sd(sds[, 1]),
        y = 1,
        colour = "red",
        label = paste0("Upper: ", n_genes_max)
      )

    return(plot)
  }
}


#' Heatmap data processing
#'
#' This function prepares the data from pre-processing
#' to be displayed in a heatmap. It takes in limits for
#' what genes to subset, what centering and standardizing
#' to perform, and what gene ID label to use.
#'
#' @param data Processed data matrix
#' @param n_genes_max Integer for for number upper limit to display in heatmap
#' @param gene_centering TRUE/FALSE subtract mean from gene rows
#' @param gene_normalize TRUE/FALSE divide by SD in gene rows
#' @param sample_centering TRUE/FALSE subtract mean from sample columns
#' @param sample_normalize TRUE/FALSE divide by SD in sample columns
#' @param all_gene_names Data frame of gene names
#' @param select_gene_id Character string designating desired ID type for
#'   heatmap labels (User_ID, ensembl_ID, symbol)
#'
#' @export
#' @return Subsetted data matrix ([n_genes_min:n_genes_max, ]) with
#'   gene IDs as the select_gene_id
#'
#' @family clustering functions
#' @family heatmaps
process_heatmap_data <- function(data,
                                 n_genes_max,
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

  if (n_genes_max > nrow(data)) {
    n_genes_max <- nrow(data)
  } else if (n_genes_max < 10) {
    n_genes_max <- 10
  }

  # Center by genes, cutoff really large values ---------
  if (gene_centering) {
    data <-
      data[1:n_genes_max, ] -
      apply(data[1:n_genes_max, ], 1, mean)
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
    data <- data[1:n_genes_max, ]
  }

  return(data)
}

#' Draw a heatmap of processed data
#'
#' Uses the package ComplexHeatmaps to draw a heatmap of the
#' processed data that has been prepared for the heatmap. The
#' returned heatmap visualizes the processed expression of the
#' gene range provided in \code{\link{process_heatmap_data}()}.
#'
#' @param data Data matrix returned from \code{\link{process_heatmap_data}()}
#' @param cluster_meth Integer indicating which clustering method to use 1 for
#'   hierarchical and 2 for kmeans.
#' @param heatmap_cutoff Numeric value for Z score max to filter data
#' @param sample_info Matrix of experiment design information from load data
#' @param select_factors_heatmap Factor group for annotation legend
#' @param dist_funs List of distance functions to use in heatmap
#' @param dist_function The selected distance function to use
#' @param hclust_function Type of clustering to perform
#' @param sample_clustering TRUE/FALSE Specify whether to cluster columns
#' @param heatmap_color_select Vector of colors for heatmap scale
#' @param row_dend TRUE/FALSE Hide row dendrogram,
#' @param k_clusters Number of clusters to use for k-means
#' @param re_run Integer for a seed to Re-run k-means with
#' @param selected_genes Character list of genes to label on the heatmap
#' @param group_pal Named list of colors and their corresponding categories
#' @param sample_color Selected colorspace color palette
#'
#' @export
#' @return Heatmap of the processed data.
#'
#' @family clustering functions
#' @family heatmaps
#' @seealso
#' * \code{\link{dist_functions}()} for available distance functions,
#' *  \code{\link{hcluster_functions}()} for available functions for hierarchical
#' @md
heatmap_main <- function(data,
                         cluster_meth,
                         heatmap_cutoff,
                         sample_info,
                         select_factors_heatmap,
                         dist_funs,
                         dist_function,
                         hclust_function,
                         sample_clustering,
                         heatmap_color_select,
                         row_dend,
                         k_clusters,
                         re_run,
                         selected_genes,
                         group_pal = NULL,
                         sample_color = NULL) {
  # Filter with max z-score
  cutoff <- median(unlist(data)) + heatmap_cutoff * sd(unlist(data))
  data[data > cutoff] <- cutoff
  cutoff <- median(unlist(data)) - heatmap_cutoff * sd(unlist(data))
  data[data < cutoff] <- cutoff
  
  # Color scale
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
  groups <- detect_groups(colnames(data))
  heat_ann <- NULL
  if (!is.null(select_factors_heatmap)) {
    if (select_factors_heatmap != "All factors") { # one factor-------
      show_legend <- TRUE
      # Annotation for groups
      if (!is.null(sample_info) && !is.null(select_factors_heatmap)) {
        if (select_factors_heatmap == "Names") {
          groups <- detect_groups(colnames(data))
          show_legend <- FALSE
        } else {
          ix <- match(select_factors_heatmap, colnames(sample_info))
          groups <- sample_info[, ix]
        }
      }
      groups_colors <- colorspace::qualitative_hcl(length(unique(groups)), 
                                                   palette = sample_color,
                                                   c = 70)
      heat_ann <- ComplexHeatmap::HeatmapAnnotation(
        Group = groups,
        col = list(Group = setNames(groups_colors, unique(groups))),
        show_annotation_name = list(Group = FALSE),
        show_legend = show_legend
      )
    } else { # more factors------------------------
      
      heat_ann <- ComplexHeatmap::HeatmapAnnotation(
        df = sample_info,
        col = group_pal,
        show_legend = TRUE
      )
    }
  }
  
  # Different heatmaps for hierarchical and k-means
  if (cluster_meth == 1) {
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
      cluster_columns = sample_clustering,
      show_column_dend = TRUE,
      show_row_dend = row_dend,
      row_dend_side = "left",
      row_dend_width = grid::unit(1, "cm"),
      top_annotation = heat_ann,
      show_row_names = FALSE,
      show_column_names = FALSE,
      heatmap_legend_param = list(
        direction = "horizontal",
        legend_width = grid::unit(6, "cm"),
        title = "Color Key",
        title_position = "topcenter"
      )
    )
  } else if (cluster_meth == 2) {
    set.seed(re_run)
    if (k_clusters > 10) {
      row_title <- 8
    } else {
      row_title <- 10
    }
    heat <- ComplexHeatmap::Heatmap(
      data,
      name = "Expression",
      col = col_fun,
      row_km = k_clusters,
      cluster_row_slices = FALSE,
      cluster_rows = TRUE,
      cluster_columns = sample_clustering,
      show_column_dend = TRUE,
      show_row_dend = row_dend,
      row_dend_side = "left",
      row_dend_width = grid::unit(1, "cm"),
      top_annotation = heat_ann,
      show_row_names = FALSE,
      show_column_names = FALSE,
      heatmap_legend_param = list(
        direction = "horizontal",
        legend_width = grid::unit(6, "cm"),
        title = "Color Key",
        title_position = "topcenter"
      ),
      row_title_gp = grid::gpar(fontsize = row_title)
    )
  }

  # mark selected genes on heatmap
  if (!is.null(selected_genes)) {
    
    if ("Top 5" %in% selected_genes){
      selected_genes <- c(rownames(data)[1:5], 
                          selected_genes[which(selected_genes != "Top 5")])
    } 
    if ("Top 10" %in% selected_genes){
      selected_genes <- c(rownames(data)[1:10], 
                          selected_genes[which(selected_genes != "Top 10")])
    } 
    if ("Top 15" %in% selected_genes){
      selected_genes <- c(rownames(data)[1:15], 
                          selected_genes[which(selected_genes != "Top 15")])
    }
    
    selected_genes <- unique(selected_genes)
    ids <- row.names(heat@matrix)[heat@row_order]
    ix <- which(ids %in% selected_genes)
    req(length(ix) > 0)
    heat <- heat + ComplexHeatmap::rowAnnotation(
      mark = ComplexHeatmap::anno_mark(
        at = ix,
        labels = ids[ix],
        link_width = ggplot2::unit(8, "points"),
        labels_gp = grid::gpar(fontsize = 9)
      )
    )
  }

  return(
    ComplexHeatmap::draw(
      heat,
      heatmap_legend_side = "bottom"
    )
  )
}

#' Draw a dendogram of data samples
#'
#' Create a clustered tree of the samples in the dataset based on hierarchical
#'   clustering
#'
#' @param tree_data Data that has been through pre-processing
#' @param gene_centering TRUE/FALSE subtract mean from gene rows
#' @param gene_normalize TRUE/FALSE divide by SD in gene rows
#' @param sample_centering TRUE/FALSE subtract mean from sample columns
#' @param sample_normalize TRUE/FALSE divide by SD in sample columns
#' @param hclust_funs Clustering functions defined in idep, \code{\link{hcluster_functions}()}
#' @param hclust_function String of chosen clustering method
#' @param dist_funs Distance functions defined in idep \code{\link{dist_functions}()}
#' @param dist_function Selected distance function
#'
#' @export
#' @return Dendogram plot of dataset samples
#'
#' @family plots
#' @family clustering functions
draw_sample_tree <- function(tree_data,
                             gene_centering,
                             gene_normalize,
                             sample_centering,
                             sample_normalize,
                             hclust_funs,
                             hclust_function,
                             dist_funs,
                             dist_function) {
  max_gene <- apply(tree_data, 1, max)
  # Remove bottom 25% lowly expressed genes, which inflate the PPC
  tree_data <- tree_data[which(max_gene > quantile(max_gene)[1]), ]
  # Center by gene
  if (gene_centering) {
    tree_data <- tree_data - apply(tree_data, 1, mean)
  }
  # Normalize by gene
  if (gene_normalize) {
    tree_data <- tree_data / apply(tree_data, 1, sd)
  }
  # Center and normalize by sample
  tree_data <- scale(
    tree_data,
    center = sample_centering,
    scale = sample_normalize
  )

  par(mar = c(5.1, 4.1, 4.1, 20))

  plot(
    stats::as.dendrogram(
      hclust_funs[[hclust_function]]
      (dist_funs[[as.numeric(dist_function)]](t(tree_data)))
    ),
    xlab = "",
    ylab = paste(
      names(dist_funs)[as.numeric(dist_function)], "(",
      hclust_function, "linkage", ")"
    ),
    type = "rectangle",
    # leaflab = "textlike"
    horiz = TRUE
  )
}

#' Draw an elbow plot for k-cluster selection
#'
#' This function takes in the processed heatmap data and
#' creates an elbow plot to guide the selection of the
#' number of clusters to create
#'
#' @param heatmap_data Matrix of processed heatmap data from
#'   \code{\link{process_heatmap_data}()}
#'
#' @export
#' @return \code{ggplot2} object of formatted elbow plot
#'
#' @family plots
#' @family clustering functions
k_means_elbow <- function(heatmap_data) {
  k.max <- 20

  validate(
    need(nrow(heatmap_data) > k.max,
      message = paste(paste("To create the elbow plot, please select at least", k.max + 1), "genes.")
    )
  )


  factoextra::fviz_nbclust(
    heatmap_data,
    kmeans,
    method = "wss",
    k.max = k.max
  ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        size = 12
      ),
      axis.text.y = ggplot2::element_text(
        size = 12
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      )
    ) +
    ggplot2::labs(
      title = "k-Means Elbow Plot",
      y = "Within Sum of Square",
      x = "Clusters"
    )
}


#' Create annotation for shiny subheatmap
#'
#' Use the heatmap data to make an annotation for the
#' submap that will also show the legend
#'
#' @param data Matrix of heatmap data
#' @param sample_info Matrix of experiment design information from load data
#' @param select_factors_heatmap Factor to group by in the samples.
#'   "All factors" will use all of the sample information.
#' @param group_pal Named list of colors and their corresponding categories
#' @param sample_color Selected colorspace color palette
#'
#' @export
#' @return A list containing a ComplexHeatmap annotation object,
#'  a ComplexHeatmap legend, list of groups, and list of group colors.
#'
#' @family clustering functions
#' @family heatmaps
sub_heat_ann <- function(data,
                         sample_info,
                         select_factors_heatmap,
                         group_pal = NULL,
                         sample_color = NULL) {
  groups <- detect_groups(colnames(data))

  if (select_factors_heatmap == "All factors") {

    # Use all factors instead of just the first one
    heat_sub_ann <- ComplexHeatmap::HeatmapAnnotation(
      df = sample_info,
      col = group_pal,
      show_legend = TRUE
    )
    # For groups, use the first factor (for backward compatibility with click info)
    groups <- sample_info[, 1]
    group_colors <- group_pal[[1]]

  } else {
    
    if (!is.null(sample_info) && !is.null(select_factors_heatmap)) {
      if (select_factors_heatmap == "Names") {
        groups <- detect_groups(colnames(data))
      } else {
        ix <- match(select_factors_heatmap, colnames(sample_info))
        groups <- sample_info[, ix]
      }
    }
    group_colors <- colorspace::qualitative_hcl(length(unique(groups)), 
                                                 palette = sample_color,
                                                 c = 70)
    group_colors <- setNames(group_colors, unique(groups))
    
    heat_sub_ann <- ComplexHeatmap::HeatmapAnnotation(
      Group = groups,
      col = list(Group = setNames(group_colors, unique(groups))),
      show_annotation_name = list(Group = FALSE),
      show_legend = FALSE
    )
  }
  if (select_factors_heatmap != "All factors" && length(unique(groups)) < 10) {
    lgd <- ComplexHeatmap::Legend(
      at = unique(groups),
      legend_gp = grid::gpar(fill = group_colors),
      nrow = 1
    )
    
  } else {
    lgd <- NULL
  }
  return(list(
    heat_sub_ann = heat_sub_ann,
    lgd = lgd,
    groups = groups,
    group_colors = group_colors
  ))
}

#' Interactive click text for subheatmap
#'
#' Create a text output to tell the user the cell
#' information for their click.
#'
#' @param click Click input from subheatmap
#' @param ht_sub Heatmap object of drawn subheatmap object
#' @param ht_sub_obj Heatmap object with mapping info
#' @param ht_pos_sub DataFrame object of position information from submap
#' @param sub_groups Vector of group labels from submap
#' @param group_colors Vector of colors for the group annotation
#' @param cluster_meth Integer indicating which clustering method to use 1 for
#'   hierarchical and 2 for kmeans.
#' @param click_data Data matrix to get the data value from
#'
#' @export
#' @return HTML code to produce a table with information
#'  about the selected cell.
#'
#' @family clustering functions
#' @family heatmaps
cluster_heat_click_info <- function(click,
                                    ht_sub,
                                    ht_sub_obj,
                                    ht_pos_sub,
                                    sub_groups,
                                    group_colors,
                                    cluster_meth,
                                    click_data) {
  pos1 <- InteractiveComplexHeatmap::getPositionFromClick(click)

  pos <- InteractiveComplexHeatmap::selectPosition(
    ht_sub,
    mark = FALSE,
    pos = pos1,
    verbose = FALSE,
    ht_pos = ht_pos_sub
  )

  row_index <- pos[1, "row_index"]
  column_index <- pos[1, "column_index"]

  if (is.null(row_index)) {
    return("Select a cell in the heatmap.")
  }

  if (cluster_meth == 1) {
    value <- click_data[row_index, column_index]
    col <- ComplexHeatmap::map_to_colors(ht_sub_obj@matrix_color_mapping, value)
    sample <- colnames(click_data)[column_index]
    gene <- rownames(click_data)[row_index]
  } else if (cluster_meth == 2) {
    clust <- pos[1, "heatmap"]
    sub_click_data <- click_data[[clust]]
    value <- sub_click_data[row_index, column_index]
    col <- ComplexHeatmap::map_to_colors(
      ht_sub_obj@ht_list$heat_1@matrix_color_mapping,
      value
    )
    sample <- colnames(sub_click_data)[column_index]
    gene <- rownames(sub_click_data)[row_index]
  }
  group_name <- sub_groups[column_index]
  group_col <- group_colors[[group_name]]

  # HTML for info table
  # Pulled from https://github.com/jokergoo/InteractiveComplexHeatmap/blob/master/R/shiny-server.R
  # Lines 1669:1678
  html <- GetoptLong::qq("
<div>
<pre>
@{gene}  Expression: @{round(value, 2)} <span style='background-color:@{col};width=50px;'>    </span>
Sample: @{sample},  Group: @{group_name} <span style='background-color:@{group_col};width=50px;'>    </span>
</pre></div>")
  HTML(html)
}

#' Draw sub heatmap from brush input
#'
#' Use the brush input from the main heatmap to
#' create a larger subheatmap.
#'
#' @param ht_brush Brush input from the main heatmap
#' @param ht Main heatmap object
#' @param ht_pos_main Position of brush on main heatmap
#' @param heatmap_data Matrix of data for the heatmap from
#'   \code{\link{process_heatmap_data}()}
#' @param sample_info Matrix of experiment design file information
#' @param select_factors_heatmap Group design to label by
#' @param cluster_meth Integer indicating which clustering method to use 1 for
#'   hierarchical and 2 for kmeans.
#' @param group_pal Named list of colors and their corresponding categories
#' @param sample_color Selected colorspace color palette
#'
#' @export
#' @return A list containing a Heatmap from the brush selection
#'  of the main heatmap, the submap data matrix, the groups for
#'  the submap, the submap legend, and data for the click info.
#'
#' @family heatmaps
#' @family clustering functions
heat_sub <- function(ht_brush,
                     ht,
                     ht_pos_main,
                     heatmap_data,
                     sample_info,
                     select_factors_heatmap,
                     cluster_meth,
                     group_pal = NULL,
                     sample_color = NULL) {
  max_gene_ids <- 2000
  lt <- InteractiveComplexHeatmap::getPositionFromBrush(ht_brush)
  pos1 <- lt[[1]]
  pos2 <- lt[[2]]

  pos <- InteractiveComplexHeatmap::selectArea(
    ht,
    mark = FALSE,
    pos1 = pos1,
    pos2 = pos2,
    verbose = FALSE,
    ht_pos = ht_pos_main
  )
  column_index <- unlist(pos[1, "column_index"])

  # Annotation, groups, and legend
  sub_heat <- sub_heat_ann(
    data = heatmap_data,
    sample_info = sample_info,
    select_factors_heatmap = select_factors_heatmap,
    group_pal = group_pal,
    sample_color = sample_color
  )
  sub_ann <- sub_heat$heat_sub_ann[column_index]
  sub_groups <- sub_heat$groups[column_index]
  lgd <- sub_heat$lgd
  group_colors <- sub_heat$group_colors

  if (cluster_meth == 1) {
    row_index <- unlist(pos[1, "row_index"])
    m <- ht@ht_list[[1]]@matrix
    if (length(row_index) > max_gene_ids) {
      show_rows <- FALSE
    } else {
      show_rows <- TRUE
    }
    submap_data <- m[row_index, column_index, drop = FALSE]
    click_data <- submap_data

    ht_select <- ComplexHeatmap::Heatmap(
      m[row_index, column_index, drop = FALSE],
      col = ht@ht_list[[1]]@matrix_color_mapping@col_fun,
      show_heatmap_legend = FALSE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = show_rows,
      top_annotation = sub_ann,
      name = "heat_1"
    )
  } else if (cluster_meth == 2) {
    sub_heats <- c()
    all_rows <- c()
    click_data <- c()

    for (i in 1:nrow(pos)) {
      all_rows <- c(all_rows, unlist(pos[i, "row_index"]))
    }

    if (length(all_rows) > max_gene_ids) {
      show_rows <- FALSE
    } else {
      show_rows <- TRUE
    }
    m <- ht@ht_list[[1]]@matrix
    submap_data <- m[all_rows, column_index, drop = FALSE]

    for (i in 1:nrow(pos)) {
      row_index <- unlist(pos[i, "row_index"])

      click_data[[paste0("heat_", i)]] <- m[row_index, column_index, drop = FALSE]

      sub_heats[[i]] <- ComplexHeatmap::Heatmap(
        m[row_index, column_index, drop = FALSE],
        col = ht@ht_list[[1]]@matrix_color_mapping@col_fun,
        show_heatmap_legend = FALSE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = show_rows,
        name = paste0("heat_", i)
      )
      if (i == 1) {
        sub_heats[[i]] <- ComplexHeatmap::add_heatmap(
          sub_ann,
          sub_heats[[i]],
          direction = "vertical"
        )
      } else if (i >= 2) {
        sub_heats[[i]] <- ComplexHeatmap::add_heatmap(
          sub_heats[[i - 1]],
          sub_heats[[i]],
          direction = "vertical"
        )
      }
      ht_select <- sub_heats[[i]]
    }
  }

  return(list(
    ht_select = ht_select,
    submap_data = submap_data,
    sub_groups = sub_groups,
    lgd = lgd,
    group_colors = group_colors,
    click_data = click_data
  ))
}

#' Create a correlation plot from heatmap data
#'
#' Creates a correlation matrix heatmap from the
#' heatmap data to demonstrate the correlation
#' between samples.
#'
#' @param data Heatmap data from \code{\link{process_heatmap_data}()}
#' @param label_pcc TRUE/FALSE  t0 Label with correlation coefficient
#' @param heat_cols Vector of colors to use in the correlation matrix
#' @param text_col Color to make the text labels in the plot
#'
#' @export
#' @return \code{ggplot2} object as heatmap of correlation matrix
#'
#' @family plots
cor_plot <- function(data,
                     label_pcc,
                     heat_cols,
                     text_col) {
  # remove bottom 25% lowly expressed genes, which inflate the PPC
  max_gene <- apply(data, 1, max)
  data <- data[which(max_gene > quantile(max_gene)[1]), ]
  low_col <- heat_cols[[1]]
  mid_col <- heat_cols[[2]]
  high_col <- heat_cols[[3]]

  melted_cormat <- reshape2::melt(round(cor(data), 2), na.rm = TRUE)

  ggheatmap <- ggplot2::ggplot(
    melted_cormat,
    ggplot2::aes(Var2, Var1, fill = value)
  ) +
    ggplot2::geom_tile(color = text_col) +
    ggplot2::scale_fill_gradient2(
      low = low_col,
      high = high_col,
      mid = mid_col,
      space = "Lab",
      limit = c(
        min(melted_cormat[, 3]),
        max(melted_cormat[, 3])
      ),
      midpoint = median(melted_cormat[, 3]),
      name = "Pearson's \nCorrelation"
    ) +
    ggplot2::theme_minimal() + # minimal theme
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        vjust = 1,
        size = 14,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = 14)
    ) +
    ggplot2::coord_fixed()

  if (label_pcc && ncol(data) < 20) {
    ggheatmap <- ggheatmap +
      ggplot2::geom_text(
        ggplot2::aes(Var2, Var1, label = value),
        color = text_col,
        size = 4
      )
  }

  ggheatmap + ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(
      color = "black",
      size = 14
    ),
    legend.text = ggplot2::element_text(
      color = "black",
      size = 9,
      angle = 0,
      hjust = .5,
      vjust = .5
    ),
    legend.title.align = 0.5,
    legend.position = "right"
  )
}


#' GET RID OF LISTS IN A DATA FRAME
data_frame_with_list <- function(data_object) {
  set_lists_to_chars <- function(x) {
    if (class(x) == "list") {
      y <- paste(unlist(x[1]), sep = "", collapse = ", ")
    } else {
      y <- x
    }
    return(y)
  }
  new_frame <- data.frame(
    lapply(data_object, set_lists_to_chars),
    stringsAsFactors = F
  )
  return(new_frame)
}


#' Prep heatmap data for download
#'
#' Prep heatmap data for download or additional analysis by merging gene ids
#' with the clusters from kmean clustering or hierarchical clustering.
#'
#' @param heatmap Heatmap object from the \code{\link{heatmap_main}()} function.
#' @param heatmap_data Matrix of heatmap data from the
#'   \code{\link{process_heatmap_data}()} function
#' @param cluster_meth Integer designating the clustering method used. 1 for
#'   hierarchical and 2 for kmeans
#'
#' @return A dataframe with the heatmap data and associated clusters.
#' @export
prep_download <- function(heatmap,
                          heatmap_data,
                          cluster_meth = c(1, 2)) {
  # kmeans clustering
  if (cluster_meth == 2) {
    row_ord <- ComplexHeatmap::row_order(heatmap)

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

    rownames(clusts) <- rownames(heatmap_data[clusts$row_order, ])
    clusts <- clusts |>
      dplyr::select(-c(row_order))

    data <- merge(heatmap_data, clusts, by = "row.names", all = TRUE)
    rownames(data) <- data$Row.names
    data <- data |>
      dplyr::select(-c(Row.names))

    return(data)

    # hierarchical clustering - NOT CURRENTLY USED
    # this code adds clusters to hierarchical clustering, but will not be
    # implemented at this time for simplicity, the function just returns the
    # heatmap data for heirarchical
    # must add num_clust back into parameters if used in future
  } else if (FALSE) {
    if (num_clust > dim(heatmap_data)[1]) {
      stop(
        paste0(
          "The number of clusters must be between 1 and ",
          dim(heatmap_data)[1]
        )
      )
    }

    hclust <- as.hclust(ComplexHeatmap::row_dend(heatmap))
    clusters <- as.data.frame(cutree(hclust, num_clust))
    colnames(clusters) <- "cluster"

    data <- merge(heatmap_data, clusters, by = "row.names", all = TRUE)
    rownames(data) <- data$Row.names
    data <- data |>
      dplyr::select(-c(Row.names))

    return(data)
  } else {
    return(heatmap_data)
  }
}

#' Prepare Word Cloud Data
#' 
#' Prepares words in pathway and corresponding frequencies for 
#' constructing word clouds
#'
#' @param gene_lists List of gene data within each cluster
#' @param cluster Selected cluster from k-means clustering
#' @param select_org Selected organism
#' @param gmt_file Optional custom GMT file
#' @param idep_data iDEP data 
#' @param gene_info Gene info from pre-processing step
#' @param cloud_go GO selected for word cloud. KEGG, GOBP, etc. 
#' @param converted Converted data from pre-processing
#'
#' @returns Returns data frame of words from pathways in the selected cluster
#' and their frequencies.
#' @export
#'
prep_cloud_data <- function(gene_lists,
                            cluster,
                            cloud_go,
                            select_org,
                            converted,
                            gmt_file,
                            idep_data,
                            gene_info){
  
  # Retrieve pathway information
  paths1 <- read_pathway_sets(
    all_gene_names_query = gene_lists[[cluster]],
    converted = converted,
    go = cloud_go,
    select_org = select_org,
    gmt_file = gmt_file,
    idep_data = idep_data,
    gene_info = gene_info
  )
  
  # If null, either error or no pathways found
  if (is.null(paths1)) {
    return("Pathways Not Found")
  }
  
  # Select description, remove path ID
  paths1 <- data.frame(
    Descr = remove_pathway_id(names(paths1$pathway_table$gene_sets), cloud_go), 
    n = t(as.data.frame(lapply(paths1$pathway_table$gene_sets, length))), 
    row.names = NULL
  )
  
  # Remove common words/punctuation
  words1 <- paths1 |>
    dplyr::mutate(Descr = gsub("[-[:punct:]]", " ", Descr),
                  Descr = gsub("\\s+", " ", Descr),
                  Descr = trimws(Descr))|>
    tidytext::unnest_tokens(word, Descr) |> # Tokenize the descriptions
    dplyr::filter(!word %in% c("pathway", "pathways"),
                  nchar(word) > 2) |>
    dplyr::anti_join(tidytext::stop_words, by = "word") |>
    dplyr::group_by(word) |>
    dplyr::summarise(n = sum(n, na.rm = TRUE)) |>
    dplyr::arrange(-n)

  colnames(words1)[2] <- paste0("Cluster", cluster)
  return(words1)
}
