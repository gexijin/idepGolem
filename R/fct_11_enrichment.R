
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
    label = "Please select a gene list first!",
    size = 15
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