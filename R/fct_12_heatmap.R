#' Heatmap for significant contrast genes
#'
#' Create a ComplexHeatmap of the processed expression data for
#' the genes that were significantly expressed in the selected
#' comparison. The data for this heatmap comes from the
#' deg_heat_data function.
#'
#' @param data Submatrix of the processed data matrix from the
#'  deg_heta_data function
#' @param bar Vector to signify a positive (1) expression fold
#'  change or a negative (-1) change
#' @param heatmap_color_select Color vector to use for the
#'  heatmap expression scale
#' @param cluster_row Boolean to indicate whether or not to cluster rows
#'  (TRUE/FALSE)
#'
#' @export
#' @return A drawn heatmap from the filtered data.
deg_heatmap <- function(data,
                        bar,
                        heatmap_color_select,
                        cluster_rows) {
  # Number of genes to show
  n_genes <- as.character(table(bar))

  data <- as.matrix(data) - apply(data, 1, mean)
  cutoff <- median(unlist(data)) + 3 * sd(unlist(data))
  data[data > cutoff] <- cutoff
  cutoff <- median(unlist(data)) - 3 * sd(unlist(data))
  data[data < cutoff] <- cutoff

  # sometimes one row is all zeroes or the same value, 
  # this causes error for the complexHeatmap
  ix <- which(abs(apply(data, 1, sd)) < 1e-20)
  if(length(ix) > 0) {
    data <- data[-1 * ix, ]
    if(!is.null(bar)) {
      bar <- bar[-1 * ix]
    }
  }

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
  group_count <- length(unique(groups))
  groups_colors <- gg_color_hue(2 + group_count)

  top_ann <- ComplexHeatmap::HeatmapAnnotation(
    Group = groups,
    col = list(
      Group = setNames(groups_colors[1:group_count], unique(groups))
    ),
    annotation_legend_param = list(
      Group = list(nrow = 1, title = NULL)
    ),
    show_annotation_name = list(Group = FALSE),
    show_legend = FALSE
  )

  row_ann <- NULL
  if (!is.null(bar)) {
    bar[bar == -1] <- "Negative"
    bar[bar == 1] <- "Positive"
    groups <- bar

    row_ann <- ComplexHeatmap::rowAnnotation(
      Change = groups,
      col = list(
        Change = setNames(
          groups_colors[(group_count + 1):length(groups_colors)], unique(groups)
        )
      ),
      annotation_legend_param = list(
        Change = list(nrow = 1, title = NULL)
      ),
      show_annotation_name = list(Change = FALSE),
      show_legend = FALSE
    )
  }

  heat <- ComplexHeatmap::Heatmap(
    data,
    name = "Expression",
    col = col_fun,
    cluster_rows = cluster_rows,
    clustering_method_rows = "average",
    clustering_distance_rows = function(x) {
      as.dist(
        1 - cor(t(x), method = "pearson")
      )
    },
    cluster_columns = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    left_annotation = row_ann,
    top_annotation = top_ann,
    show_row_names = FALSE,
    show_column_names = FALSE,
    heatmap_legend_param = list(
      direction = "horizontal",
      legend_width = grid::unit(6, "cm"),
      title = "Color Key",
      title_position = "topcenter"
    )
  )
  return(
    heatmap = ComplexHeatmap::draw(
      heat,
      heatmap_legend_side = "bottom"
    )
  )
}

#' Plot brush selection from main heatmap
#'
#' Create a ComplexHeatmap object from the User brush selection
#' that is a sub plot of the main plot.
#'
#' @param ht_brush Input from the user creating a brush selection
#'  on the main heatmap
#' @param ht Main heatmap from the deg_heatmap function
#' @param ht_pos_main Main heatmap position information to use
#'  for the sub heatmap
#' @param heatmap_data Original data matrix that was plotted in
#'  the main heatmap
#' @param bar, groups of genes for colar bar on the left side
#' @param all_gene_names Data matrix of all the mapped gene names
#'
#' @export
#' @return A ComplexHeatmap object of the brushed selection from
#'  the main heatmap.
deg_heat_sub <- function(ht_brush,
                         ht,
                         ht_pos_main,
                         heatmap_data,
                         bar,
                         all_gene_names,
                         select_gene_id) {
  max_genes <- 2000
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

  # Annotations ----------
  column_groups <- detect_groups(colnames(heatmap_data))
  group_count <- length(unique(column_groups))
  groups_colors <- gg_color_hue(2 + group_count)

  top_ann <- ComplexHeatmap::HeatmapAnnotation(
    Group = column_groups,
    col = list(
      Group = setNames(
        groups_colors[1:group_count],
        unique(column_groups)
      )
    ),
    annotation_legend_param = list(
      Group = list(nrow = 1, title = NULL)
    ),
    show_annotation_name = list(Group = FALSE),
    show_legend = TRUE
  )

  row_ann <- NULL
  if (!is.null(bar)) {
    bar[bar == -1] <- "Down"
    bar[bar == 1] <- "Up"
    row_groups <- bar

    row_ann <- ComplexHeatmap::rowAnnotation(
      Change = row_groups,
      col = list(
        Change = setNames(
          groups_colors[(group_count + 1):length(groups_colors)],
          unique(row_groups)
        )
      ),
      annotation_legend_param = list(
        Change = list(nrow = 1, title = NULL)
      ),
      show_annotation_name = list(Change = FALSE),
      show_legend = TRUE
    )
    group_col_return <- setNames(
      groups_colors,
      c(unique(column_groups), unique(row_groups))
    )
  } else {
    group_col_return <- setNames(
      groups_colors,
      c(unique(column_groups))
    )
  }



  # End annotation ---------

  column_index <- unlist(pos[1, "column_index"])
  row_index <- unlist(pos[1, "row_index"])
  top_ann <- top_ann[column_index]
  if (!is.null(bar)) {
    row_ann <- row_ann[row_index]
  }
  column_groups <- column_groups[column_index]
  m <- ht@ht_list[[1]]@matrix

  bar_return <- bar[row_index]

  if (length(row_index) > max_genes) {
    show_rows <- FALSE
  } else {
    show_rows <- TRUE
  }


  if (ncol(all_gene_names) == 3) {
    genes <- rowname_id_swap(
      data_matrix = m[row_index, column_index, drop = FALSE],
      all_gene_names = all_gene_names,
      select_gene_id = select_gene_id
    )
  } else {
    genes <- m[row_index, column_index, drop = FALSE]
  }
  submap_data <- genes

  ht_select <- ComplexHeatmap::Heatmap(
    genes,
    col = ht@ht_list[[1]]@matrix_color_mapping@col_fun,
    show_heatmap_legend = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = show_rows,
    top_annotation = top_ann,
    left_annotation = row_ann,
    name = "heat_1"
  )

  # Show a subset of gene names when more than 50 is selected.
  # causes problem with returned click info.
  #  if (length(row_index) > max_genes) {
  #    loci <- seq(
  #      from = 1,
  #      to = nrow(genes),
  #      by = round(nrow(m) / 30, 0)
  #    )
  #    anno <- ComplexHeatmap::anno_mark(
  #      at = loci,
  #      labels = row.names(m)[loci],
  #      which = "row",
  #      labels_gp = grid::gpar(fontsize = 10),
  #      padding = ggtree::unit(.5, "mm")
  #    )
  #    ht_select <- ht_select + ComplexHeatmap::rowAnnotation(mark = anno)
  #  }

  return(list(
    ht_select = ht_select,
    submap_data = submap_data,
    group_colors = group_col_return,
    column_groups = column_groups,
    bar = bar_return
  ))
}

#' HTML code for sub-heatmap selected cell
#'
#' Create HTML code for a cell of information on the cell of the
#' sub-heatmap that the User clicks on. The cell contains the
#' expression value, the sample, the gene, the group and the
#' direction of the fold change.
#'
#' @param click Information fro what cell is clicked in the
#'  sub-heatmap
#' @param ht_sub The drawn sub-heatmap
#' @param ht_sub_obj The sub-heatmap ComplexHeatmap object
#' @param ht_pos_sub Position information for the sub-heatmap
#' @param sub_groups Vector of the groups that the samples
#'  belong to
#' @param group_colors The color of the top annotation that
#'  is used for each group and the side annotation that denotes
#'  the direction of the expression regulation
#' @param bar Vector to signify a positive (1) expression fold
#'  change or a negative (-1) change
#' @param data Sub data matrix that is plotted in the sub-heatmap
#'
#' @export
#' @return HTML code that will be used in the shiny UI to tell
#'  the user the information of the cell they selected.
deg_click_info <- function(click,
                           ht_sub,
                           ht_sub_obj,
                           ht_pos_sub,
                           sub_groups,
                           group_colors,
                           bar,
                           data) {
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

  value <- data[row_index, column_index]
  col <- ComplexHeatmap::map_to_colors(ht_sub_obj@matrix_color_mapping, value)
  sample <- colnames(data)[column_index]
  gene <- rownames(data)[row_index]

  group_name <- sub_groups[column_index]
  group_col <- group_colors[[group_name]]

  # HTML for info table
  # Pulled from https://github.com/jokergoo/InteractiveComplexHeatmap/blob/master/R/shiny-server.R
  # Lines 1669:1678
  p <- "
<div>
<pre>
@{gene}
Value: @{round(value, 2)} <span style='background-color:@{col};width=50px;'>    </span>
Sample: @{sample}
Group: @{group_name} <span style='background-color:@{group_col};width=50px;'>    </span>
"

  if (!is.null(bar)) {
    up_down <- bar[row_index]
    up_down_col <- group_colors[[up_down]]
    p <- paste0(
      p,
      "Regulation: @{up_down} <span style='background-color:@{up_down_col};width=50px;'>    </span>"
    )
  }
  p <- paste0(p, "</pre></div>")
  html <- GetoptLong::qq(p)

  return(HTML(html))
}
