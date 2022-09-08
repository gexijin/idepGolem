#' Dendogram of enriched pathways
#' 
#' Create a dendogram plot of the enriched pathways to illustrate
#' which paths contain similar genes.
#' 
#' @param go_table Enrichment table from the pathway analysis, the last column Genes
#' contains lists
#' @param group  Selected group int the Direction column 
#' @param right_margin Control the size of the dendogram labels
#' 
#' @export
#' @return A dendogram plot that shows the users what pathways are
#'  that are enriched share genes.
enrichment_tree_plot <- function(
  go_table,
  group,
  right_margin = 10
) {
  # a program for ploting enrichment results by highlighting the similarities among terms
  # must have columns: Direction, adj.Pval   Pathways Genes
  #  Direction  adj.Pval  nGenes  Pathways    Genes
  #Down regulated  3.58E-59  131  Ribonucleoprotein complex biogenesis  36  Nsun5 Nhp2 Rrp15 
  #Down regulated  2.55E-57  135  NcRNA metabolic process  23  Nsun5 Nhp2 Rrp15 Emg1 Ddx56 Rsl1d1
  # Up or down regulation is color-coded
  # gene set size if represented by the size of marker

  req(!is.null(go_table))
  req(!is.null(group))
  data <- go_table
  if(class(data) != "data.frame") {
    return(NULL)
  }
  # only one term or less
  if(nrow(data) <=1 || is.null(data)) {
    return(NULL)
  }

  # only use selected group
  if(group != "All Groups") {
    data <- data[data$Direction == group, ]
  }
  # this is unneccessary, but works
  gene_lists <- lapply(
    data$Genes,
    function(x) unlist(strsplit(as.character(x), ", "))
  )
  names(gene_lists) <- data$Pathways

  # Compute overlaps percentage--------------------
  n <- length(gene_lists)
  w <- matrix(NA, nrow = n, ncol = n)
  
  # Compute overlaps among all gene lists
  for(i in 1:n) {
    for (j in i:n) {
      u <- unlist(gene_lists[i])
      v <- unlist(gene_lists[j])
      w[i, j] <- length(intersect(u, v)) / length(unique(c(u, v)))
    }
  }
  # The lower half of the matrix filled in based on symmetry
  for(i in 1:n) {
    for(j in 1:(i-1)) {
      w[i, j] <- w[j, i]
    } 
  }

  Terms <- paste(
    sprintf("%-1.0e",
    as.numeric(data$adj_p_val)), 
    names(gene_lists)
  )
  rownames(w) <- Terms
  colnames(w) <- Terms

  # A large margin for showing 
  par(mar = c(0, 0, 1, right_margin))

  dend <- stats::as.dist(1 - w) |>
    stats::hclust(method = "average")
  # Permutated order of leaves 
  ix <- dend$order 

  leaf_type <- as.factor(data$Direction[ix])
  leaf_colors <- gg_color_hue(length(unique(data$Direction)))
  # Leaf size represent P values
  leaf_size <- -log10(as.numeric(data$adj_p_val[ix]))
  leaf_size <- 1.5 * leaf_size / max(leaf_size) + .2
  
  dend |> 
    stats::as.dendrogram(hang = -1) |>
    # Type of marker
    dendextend::set("leaves_pch", 19) |>
    # Size
    dendextend::set("leaves_cex", leaf_size) |>
    # up or down genes
    dendextend::set("leaves_col", leaf_colors[leaf_type]) |>
    plot(horiz = TRUE)
  
  # Add legend using a second layer
  par(lend = 1)
  add_legend(
    "top",
    pch = 19,
    col = leaf_colors,
    legend = levels(leaf_type),
    bty = "n",
    horiz = T 
  )

  return(recordPlot())
}

#' Generate barplot for enrichment results
#' 
#' Used to translate enrichment analysis results into bar plot like those in ShinyGO
#' 
#' @param enrichment_dataframe A data frame of pathways, P values ect.
#' @param pathway_order Sort pathway list
#' @param order_x x-axis order
#' @param plot_size size mapping vairalbe
#' @param plot_color color mapping vairable
#' @param plot_font_size font size
#' @param plot_marker_size marker size
#' @param plot_high_color High color
#' @param plot_low_color  low color
#' @param threshold_wald_test whether to use threshold-based Wald test 
#' @param chart_type barplot, lollipop, or dotplot
#' @param aspect_ratio aspect ratio of plot
#' @param select_cluster  which cluster is selected
#' 
#' @export
#' @return A ggplot2 object
enrich_barplot <- function(
  enrichment_dataframe,
  pathway_order,
  order_x,
  plot_size,
  plot_color,
  plot_font_size,
  plot_marker_size,
  plot_high_color,
  plot_low_color,
  chart_type,
  aspect_ratio,
  select_cluster
) {

  if(is.null(enrichment_dataframe)) {
    return(NULL)
  }
  req(pathway_order)
  req(order_x)
  req(plot_size)
  req(plot_color)
  req(plot_font_size)
  req(plot_marker_size)
  req(plot_high_color)
  req(plot_low_color)
  req(chart_type)
  req(aspect_ratio)
  req(select_cluster)

  fake <- data.frame(a = 1:3, b = 1:3)
  blank <- ggplot2::ggplot(fake, ggplot2::aes(x = a, y = b)) +
    ggplot2::geom_blank() +
    ggplot2::annotate(
    "text",
    x = 2,
    y = 2,
    label = "Select a group of genes from above",
    size = 13
    ) +
    ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )

  df <- enrichment_dataframe

  # filter by group
  if(select_cluster != "All Groups") {
    df <- subset(df, group == select_cluster)
  }

  # if "All Groups"
  if(length(unique(df$group)) > 1) {
    return(blank)
  }

  # Remove spaces in col names
  colnames(df) <- gsub(" ", "", colnames(df))

  df <- subset(df, select = -group)
  colnames(df)[1:5] <- c(
    "EnrichmentFDR", "nGenes",
    "PathwayGenes", "FoldEnrichment", "Pathway"
  )

  # why some pathways appear twice?
  df <- df[!duplicated(df$Pathway), ]


  df$EnrichmentFDR <- as.numeric(df$EnrichmentFDR)
  df$nGenes <- as.numeric(df$nGenes)
  df$PathwayGenes <- as.numeric(df$PathwayGenes)
  df$FoldEnrichment <- as.numeric(df$FoldEnrichment)

  x  <- order_x
  size  <- plot_size
  color_by <- plot_color
  font_size <- plot_font_size
  marker_size <- plot_marker_size
  # validate values; users can input any numeric value outside the range
  if(font_size < 1 || font_size >= 20) {
    font_size <- 12
  }
  if(marker_size < 0 || marker_size > 20) {
    marker_size <- 4
  }

  # convert to vector so that we can look up the readable names of columns 
  columns <- unlist(column_selection)

  df$EnrichmentFDR <- -log10(df$EnrichmentFDR)
  ix <- which(colnames(df) == pathway_order)

  # sort the pathways
  if(ix >0 && ix < dim(df)[2]) {
    df <- df[order(df[, ix], decreasing = TRUE), ]
  }

  # convert to factor so that the levels are not reordered by ggplot2
  df$Pathway <- factor(df$Pathway, levels = rev(df$Pathway))

  p <- ggplot2::ggplot(df,
    ggplot2::aes_string(
      x = x,
      y = "Pathway",
      size = size,
      color = color_by)
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_continuous(
      low = plot_low_color, 
      high = plot_high_color,
      name = names(columns)[columns == color_by],
      guide = ggplot2::guide_colorbar(reverse = TRUE)
    ) +
    ggplot2::scale_size(range = c(1, marker_size)) +
    ggplot2::xlab(names(columns)[columns == x]) +
    ggplot2::ylab(NULL) +
    ggplot2::guides(
      size  = ggplot2::guide_legend(
        order = 2, 
        title = names(columns)[columns == size]
      ),
      color = ggplot2::guide_colorbar(order = 1)
    ) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = font_size),
      axis.title = ggplot2::element_text(size = 12)
    ) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 12), # decrease legend font
      legend.text = ggplot2::element_text(size = 12)
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(size = 5))
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 5))
    )

  if(chart_type == "dotplot") {
    p <- p
  } else if(chart_type == "lollipop") {
    p <- p +
    ggplot2::geom_segment(
      ggplot2::aes_string(
        x = 0,
        xend = x,
        y = "Pathway",
        yend = "Pathway"
      ),
      size=1
    )
  } else if(chart_type == "barplot") {
    p <- ggplot2::ggplot(df, 
    ggplot2::aes_string(x = x, y = "Pathway", fill = color_by)) +
    ggplot2::geom_col(width = 0.8,
      position = ggplot2::position_dodge(0.7)
    ) +
    ggplot2::scale_fill_continuous(
      low = plot_low_color,
      high=plot_high_color,
      name = names(columns)[columns == color_by],
      guide = ggplot2::guide_colorbar(reverse = TRUE)
    ) +
    ggplot2::xlab(names(columns)[columns == x]) +
    ggplot2::ylab(NULL) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = font_size))
  }
  return(p)
}

#' VisNetwork data
#' 
#' Create VisNetwork data that can be inputted in the vis_network_plot
#' function to create an interactive network of enriched pathways.
#' 
#' @param network GO table from the pathway analysis
#' @param up_down_reg_deg Plot just up/down or both
#' @param wrap_text_network_deg Wrap the text from the pathway description
#' @param layout_vis_deg BUtton to reset the layout of the network
#' @param edge_cutoff_deg P-value to cutoff enriched pathways
#' 
#' @export
#' @return Data that can be inputted in the vis_network_plot function
#'  to create an interactive network.
network_data <- function(
  network,
  up_down_reg_deg,
  wrap_text_network_deg,
  layout_vis_deg,
  edge_cutoff_deg
) {
  if(up_down_reg_deg != "All Groups") {
    network <- network[network$Direction == up_down_reg_deg, ]
  }
  if(dim(network)[1] == 0) {
    return(NULL)
  }
  
  if(wrap_text_network_deg) {
    # Wrap long pathway names using default width of 30
    network$Pathways <- wrap_strings(network$Pathways)
  }

  g <- enrichment_network(
    network,
    layout_button = layout_vis_deg,
    edge_cutoff = edge_cutoff_deg
  )

  if(is.null(g)) {
    return(NULL)
  }

  vis_net <- visNetwork::toVisNetworkData(g)
    
  # Color codes: https://www.rapidtables.com/web/color/RGB_Color.html
  vis_net$nodes$shape <- "dot"
    
  vis_net$nodes$size <- 5 + vis_net$nodes$size^2 
    
  return(vis_net)
}



#' ENRICHMENT NETWORK FUNCTION
enrichment_network <- function(go_table,
                               layout_button = 0,
                               edge_cutoff = 5) {
  req(!is.null(go_table))
  gene_lists <- lapply(go_table$Genes, function(x) unlist(strsplit(as.character(x), " ")))
  names(gene_lists) <- go_table$Pathways
#  go_table$Direction <- gsub(" .*", "", go_table$Direction)

  g <- enrich_net(
    data = go_table,
    gene_set = gene_lists,
    node_id = "Pathways",
    num_char = 100,
    p_value = "adj.Pval",
    p_value_cutoff = 1,
    degree_cutoff = 0,
    n = 200,
    group = go_table$Direction,
    group_color = gg_color_hue(2 + length(unique(go_table$Direction))),
    vertex.label.cex = 1,
    vertex.label.color = "black",
    show_legend = FALSE,
    layout_button = layout_button,
    edge_cutoff = edge_cutoff
  )
}

#' numChar=100 maximum number of characters
#' n=200  maximum number of nodes
#' degree.cutoff = 0    Remove node if less connected
#' from PPInfer
enrich_net <- function(data,
                       gene_set,
                       node_id,
                       node_name = node_id,
                       p_value,
                       n = 50,
                       num_char = NULL,
                       p_value_cutoff = 0.05,
                       edge_cutoff = 0.05,
                       degree_cutoff = 0,
                       edge_width = function(x) {
                         5 * x^2
                       },
                       node_size = function(x) {
                         2.5 * log10(x)
                       },
                       group = FALSE,
                       group_color = c("green", "red"),
                       group_shape = c("circle", "square"),
                       legend_parameter = list("topright"),
                       show_legend = TRUE,
                       plotting = TRUE,
                       layout_button = 0,
                       ...) {
  set.seed(layout_button)
  data <- data.frame(data, group)
  colnames(data)[length(colnames(data))] <- "Group"
  data <- data[as.numeric(data[, "adj_p_val"]) < p_value_cutoff, ]
  data <- data[order(data[, "adj_p_val"]), ]
  n <- min(nrow(data), n)
  if (n == 0) {
    stop("no enriched term found...")
  }
  data <- data[1:n, ]
  index <- match(data[, node_id], names(gene_set))
  gene_sets_list <- list()
  for (i in 1:n) {
    gene_sets_list[[i]] <- gene_set[[index[i]]]
  }
  names(gene_sets_list) <- data[, node_name]

  if (is.null(num_char)) {
    num_char <- max(nchar(as.character(data[, node_name])))
  } else {
    if (length(unique(substr(data[, node_name], 1, num_char))) < nrow(data)) {
      num_char <- max(nchar(as.character(data[, node_name])))
      message("Note : numChar is too small.", "\n")
    }
  }
  data[, node_name] <- paste(
    substr(data[, node_name], 1, num_char),
    ifelse(nchar(as.character(data[, node_name])) > num_char, "...", ""),
    sep = ""
  )
  w <- matrix(NA, nrow = n, ncol = n)

  no_overlap <- TRUE
  for (i in 1:n) {
    for (j in i:n) {
      u <- unlist(gene_sets_list[i])
      v <- unlist(gene_sets_list[j])
      w[i, j] <- length(intersect(u, v)) / length(unique(c(u, v)))
      if(i != j && w[i,j] > edge_cutoff){
        no_overlap <- FALSE
      }
    }
  }

  if(no_overlap) {
    return(NULL)
  }
  list_edges <- stack(data.frame(w))
  list_edges <- cbind(
    list_edges[, 1],
    rep(data[, node_name], n),
    rep(data[, node_name], each = n)
  )
  list_edges <- list_edges[list_edges[, 2] != list_edges[, 3], ]
  list_edges <- list_edges[!is.na(list_edges[, 1]), ]
  g <- igraph::graph.data.frame(list_edges[, -1], directed = FALSE)
  igraph::E(g)$width <- edge_width(as.numeric(list_edges[, 1]))
  igraph::V(g)$size <- node_size(lengths(gene_sets_list))
  g <- igraph::delete.edges(g, igraph::E(g)[as.numeric(list_edges[, 1]) < edge_cutoff])
  index_deg <- igraph::degree(g) >= degree_cutoff
  g <- igraph::delete.vertices(g, igraph::V(g)[!index_deg])
  data <- data[index_deg, ]
  index <- index[index_deg]
  if (length(igraph::V(g)) == 0) {
    stop("no categories greater than degree_cutoff...")
  }
  n <- min(nrow(data), n)
  data <- data[1:n, ]
  group_level <- sort(unique(group))
  p_values <- log10(as.numeric(data[, "adj_p_val"]))

  for (i in 1:length(group_level)) {
    index <- data[, "Group"] == group_level[i]
    igraph::V(g)$shape[index] <- group_shape[1] #group_shape[i]
    group_p_values <- p_values[index]

    if (length(group_p_values) > 0) {
      if (max(group_p_values) == min(group_p_values)) {
        igraph::V(g)$color[index] <- grDevices::adjustcolor(
          group_color[i],
          alpha.f = 0.5
        )
      } else {
        igraph::V(g)$color[index] <- sapply(
          1 - .9 * (group_p_values - min(group_p_values)) /
            (max(group_p_values) - min(group_p_values)),
          function(x) {
            grDevices::adjustcolor(group_color[i], alpha.f = .1 + x)
          }
        )
      }
    }
  }
  if (plotting) {
    plot(g, , vertex.label.dist = 1.2, ...)
    if (show_legend) {
      legend.parameter$legend <- group.level
      legend.parameter$text.col <- group.color
      legend.parameter$bty <- "n"
      do.call(legend, legend.parameter)
    }
  }

  return(g)
}

#' VIS NETWORK FUNCTION
vis_network_plot <- function(network_data) {
  if(is.null(network_data)) {
    return(NULL)
  }
  visNetwork::visNetwork(
    nodes = network_data$nodes,
    edges = network_data$edges,
    height = "700px",
    width = "700px"
  ) |>
    visNetwork::visIgraphLayout(layout = "layout_with_fr") |>
    visNetwork::visNodes(
      color = list(
        border = "#000000",
        highlight = "#FF8000"
      ),
      font = list(
        color = "#000000",
        size = 20
      ),
      borderWidth = 1,
      shadow = list(
        enabled = TRUE,
        size = 10
      )
    ) |>
    visNetwork::visEdges(
      shadow = FALSE,
      color = list(
        color = "#A9A9A9",
        highlight = "#FFD700"
      )
    ) |>
    visNetwork::visExport(
      type = "jpeg",
      name = "export-network",
      float = "left",
      label = "Export as an image (only what's visible on the screen!)",
      background = "white",
      style = ""
    )
}
