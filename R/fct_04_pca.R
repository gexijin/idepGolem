#' fct_05_pca.R This file holds all of the main data analysis functions
#' associated with fifth tab of the iDEP website.
#'
#'
#' @section fct_05_pca.R functions:
#' \code{my_pgsea}
#'
#'
#' @name fct_05_pca.R
NULL

#' Prepare PC data
#' @param data Data that has been through pre-processing
#' @param sample_info Matrix array with experiment info
#'
#' @export
#' @return pca data ready for plotting
get_pc <- function(data,
                   sample_info) {

  pca.object <- prcomp(t(data))

  # 5 pc's or number of columns if <5
  npc <- min(5, ncol(data))
  pcaData <- as.data.frame(pca.object$x[, 1:npc])

  groups <- detect_groups(sample_names = colnames(data), sample_info = sample_info)
  # Missing design clause
  if (is.null(sample_info)) {
    pcaData <- cbind(pcaData, detect_groups(colnames(data), sample_info))
  } else {
    pcaData <- cbind(pcaData, detect_groups(colnames(data), sample_info), sample_info)
  }
  # dim(pcaData)[2]
  colnames(pcaData)[npc + 1] <- "Names"

  # return data frame
  return(pcaData)
}

#' Get variance percentages
#' @param data Data
#' @export
#' @return importance of each pc
get_pc_variance <- function(data) {
  # subset data if more than 100 columns
  if (ncol(data) > 100) {
    part <- 1:100
    data <- data[, part]
  }
  # pca
  pca.object <- prcomp(t(data))

  # var proportions vector
  prop_var <- summary(pca.object)$importance[2, ] * 100

  return(prop_var |> round(1))
}



#' Principal component analysis plot
#'
#' Draw a PCA plot with designated PCA components on axis
#'
#' @param data Matrix of gene data that has been through
#'  \code{\link{pre_process}()}
#' @param sample_info Matrix of sample information from experiment design file
#' @param PCAx Integer designating the PC to be plotted on the x axis
#' @param PCAy Integer designating the PC to be plotted on on the y axis
#' @param selected_color String designating factor to color points by. Should be
#'  one of the design factors from the design file or "Names" as default which
#'  automatically detects groups from gene data file
#' @param selected_shape String designating factor to shape points by.
#'  Should be one of the design factors from the design file or "Names" as
#'  default which automatically detects groups from gene data file
#' @param plots_color_select Vector of colors for plots
#'
#' @export
#' @return A \code{ggplot} object as a PCA plot
#'
#' @seealso \code{\link{PCA_plot_3d}()} for three-dimensional version
#'
#' @family PCA functions
#' @family plots
PCA_plot <- function(data,
                     sample_info,
                     PCAx = 1,
                     PCAy = 2,
                     selected_color = "Names",
                     selected_shape = "Names",
                     plots_color_select) {
  # no design file
  if (is.null(selected_color)) {
    selected_color <- "Names"
  }
  if (is.null(selected_shape)) {
    selected_shape <- "Names"
  }
  memo <- ""

  if (ncol(data) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }

  # get groups
  groups <- detect_groups(sample_names = colnames(data), sample_info = sample_info)

  # get data
  pcaData <- get_pc(data, sample_info)

  levels <- length(unique(groups))
  # plot color scheme
  color_palette <- generate_colors(n = levels, palette_name = plots_color_select)
  
  nshapes <- (length(unique(groups)) / 8) + 1

  # hide legend for large or no groups levels
  if (levels <= 1 | levels > 20) {
    group_fill <- NULL
    legend <- "none"
  } else {
    group_fill <- groups
    legend <- "right"
  }

  # adjust axis label size
  if (ncol(data) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }


  # Set point & text size based on number of sample
  point_size <- 6
  if (ncol(data) >= 40) {
    point_size <- 3
  }
  
  pcaData$tooltip_text <- paste0(
    "Name: ", rownames(pcaData),
    "<br>PC",PCAx,": ", round(pcaData[[paste0("PC", PCAx)]], 2),
    "<br>PC",PCAy,": ", round(pcaData[[paste0("PC", PCAy)]], 2)
  )

  plot_PCA <- ggplot2::ggplot(
    data = pcaData,
    mapping = ggplot2::aes_string(
      x = paste0("PC", PCAx),
      y = paste0("PC", PCAy),
      color = selected_color,
      shape = selected_shape,
      group = selected_shape,
      text = "tooltip_text" # Add this line to specify the tooltip content
    )
  ) +
    # Preferred shapes
    ggplot2::scale_shape_manual(
      values = c(
        rep(c(15, 16, 18, 0, 1, 3:5), nshapes)
      )
    ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right", # TODO no legend for large data
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
      ),
      axis.text.y = ggplot2::element_text(
        size = 16
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      )
    ) +
    ggplot2::labs(
      title = memo,
      y = "Dimension 2",
      x = "Dimension 1"
    ) +
    ggplot2::scale_color_manual(values = color_palette)
    #+
  # removed - causes plot legend to be missing shapes
  # ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = 15)))9


  # selected principal components
  PCAxy <- c(as.integer(PCAx), as.integer(PCAy))
  percentVar <- get_pc_variance(data)[PCAxy] # round(100 * summary(pca.object)$importance[2, PCAxy], 0)
  plot_PCA <- plot_PCA + ggplot2::xlab(paste0("PC", PCAx, ": ", percentVar[1], "% Variance"))
  plot_PCA <- plot_PCA + ggplot2::ylab(paste0("PC", PCAy, ": ", percentVar[2], "% Variance"))
  return(plot_PCA)
}





#' Principal component analysis 3D plot
#'
#' Draw a 3D PCA plot with designated PCA components on axis
#'
#' @param data Matrix of gene data that has been through
#'  \code{\link{pre_process}()}
#' @param sample_info Matrix of sample information from experiment design file
#' @param PCAx Integer designating the PC to be plotted on the x axis
#' @param PCAy Integer designating the PC to be plotted on the y axis
#' @param PCAz Integer designating the PC to be plotted on the z axis
#' @param selected_color String designating factor to color points by. Should be
#'  one of the design factors from the design file or "Names" as default which
#'  automatically detects groups from gene data file
#' @param selected_shape String designating factor to shape points by.
#'  Should be one of the design factors from the design file or "Names" as
#'  default which automatically detects groups from gene data file
#' @param plots_color_select Vector of colors for plots
#'
#' @export
#' @return Formatted PCA plot
#'
#' @family PCA functions
#' @seealso \code{\link{PCA_plot}()}
#'
PCA_plot_3d <- function(data,
                        sample_info,
                        PCAx = 1,
                        PCAy = 2,
                        PCAz = 3,
                        selected_color = "Names",
                        selected_shape = "Names",
                        plots_color_select) {
  # no design file
  if (is.null(selected_color)) {
    selected_color <- "Names"
  }
  if (is.null(selected_shape)) {
    selected_shape <- "Names"
  }
  counts <- data
  memo <- ""

  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }

  pca.object <- prcomp(t(data))

  # 5 pc's or number of columns if <5
  npc <- min(5, ncol(data))
  pcaData <- as.data.frame(pca.object$x[, 1:npc])

  groups <- detect_groups(sample_names = colnames(data), sample_info = sample_info)
  
  levels <- length(unique(groups))
  
  # Missing design clause
  if (is.null(sample_info)) {
    pcaData <- cbind(pcaData, detect_groups(colnames(data), sample_info))
  } else {
    pcaData <- cbind(pcaData, detect_groups(colnames(data), sample_info), sample_info)
  }
  colnames(pcaData)[npc + 1] <- "Names"
  
  nshapes <- (length(unique(groups)) / 8) + 1

  # plot color scheme
  color_palette <- generate_colors(n = levels, palette_name = plots_color_select)
  
  # selected principal components
  PCAxyz <- c(as.integer(PCAx), as.integer(PCAy), as.integer(PCAz))
  percentVar <- get_pc_variance(data)[PCAxyz]
  plot_PCA <- plotly::plot_ly(pcaData,
    x = pcaData[, as.integer(PCAx)],
    y = pcaData[, as.integer(PCAy)],
    z = pcaData[, as.integer(PCAz)],
    color = pcaData[, selected_color],
    colors = color_palette[sort(as.factor(pcaData[, selected_color]))],
    symbol = pcaData[,selected_shape],
    symbols = rep(c("square", "circle", "diamond", "square-open", "circle-open",
                    "cross", "x", "diamond-open"), nshapes),
    text = rownames(pcaData),
    hovertemplate = ~paste(
      "<b>",rownames(pcaData),"</b><br><br>",
      "PC ", PCAy, ":%{y:.3f}<br>",
      "PC ", PCAx, ":%{x:.3f}<br>",
      "PC ", PCAz, ":%{z:.3f}<br>",
      "<extra></extra>"
    ),
    type = "scatter3d",
    mode = "markers"
    )
  plot_PCA <- plotly::layout(
    p = plot_PCA,
    legend = list(title = list(text = "Names")),
    plot_bgcolor = "#e5ecf6",
    scene = list(
      xaxis = list(
        title = paste0("PC", PCAx, ": ", percentVar[1], "% Variance"),
        automargin = TRUE
      ),
      yaxis = list(
        title = paste0("PC", PCAy, ": ", percentVar[2], "% Variance"),
        automargin = TRUE
      ),
      zaxis = list(
        title = paste0("PC", PCAz, ": ", percentVar[3], "% Variance"),
        automargin = TRUE
      )
    )
  )
  return(plot_PCA)
}

#' t-SNE plot
#'
#' Draw a t-distributed stochastic neighbor embedding (t-SNE) plot
#'
#' @param data Matrix of gene data that has been through
#'  \code{\link{pre_process}()}
#' @param sample_info Matrix of sample information from experiment design file
#' @param selected_color String designating factor to color points by. Should be
#'  one of the design factors from the design file or "Names" as default which
#'  automatically detects groups from gene data file
#' @param selected_shape String designating factor to shape points by.
#'  Should be one of the design factors from the design file or "Names" as
#'  default which automatically detects groups from gene data file
#' @param plots_color_select Vector of colors for plots
#'
#' @export
#' @return A \code{ggplot} object formatted t-SNE plot
#'
#' @family PCA functions
#' @family plots
#'
t_SNE_plot <- function(data,
                       sample_info,
                       selected_color = "Names",
                       selected_shape = "Names",
                       plots_color_select) {
  # no design file
  if (is.null(selected_color)) {
    selected_color <- "Names"
  }
  if (is.null(selected_shape)) {
    selected_shape <- "Names"
  }

  counts <- data
  memo <- ""

  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }


  x <- data
  y <- sample_info
  tsne <- Rtsne::Rtsne(t(x), dims = 2, perplexity = 1, verbose = FALSE, max_iter = 400)
  pcaData <- as.data.frame(tsne$Y)
  rownames(pcaData) <- rownames(t(x))
  
  # Missing design clause
  if (is.null(sample_info)) {
    pcaData <- cbind(pcaData, detect_groups(colnames(x), y))
  } else {
    pcaData <- cbind(pcaData, detect_groups(colnames(x), y), sample_info)
  }

  colnames(pcaData)[1:3] <- c("x1", "x2", "Names")
  
  pcaData$tooltip_text <- paste0(
    "Name: ", rownames(pcaData),
    "<br>Dimension 1: ", round(pcaData$x1, 2),
    "<br>Dimension 2: ", round(pcaData$x2, 2)
  )

  nshapes <- (length(unique(as.factor(pcaData$Names))) / 8) + 1
  
  # Set point size based on number of sample
  point_size <- 6
  if (ncol(x) >= 40) {
    point_size <- 3
  }

  # plot color scheme
  color_palette <- generate_colors(n = nlevels(as.factor(pcaData$Names)), palette_name = plots_color_select)

  # Generate plot
  plot_t_SNE <- ggplot2::ggplot(
    data = pcaData,
    ggplot2::aes_string(
      x = "x1",
      y = "x2",
      color = selected_color,
      shape = selected_shape,
      text = "tooltip_text"
    )
  ) +
    ggplot2::geom_point(size = point_size) +
    # Preferred shapes
    ggplot2::scale_shape_manual(values = rep(c(15, 16, 18, 0, 1, 3:5), nshapes)) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
      ),
      axis.text.y = ggplot2::element_text(
        size = 16
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      )
    ) +
    ggplot2::labs(
      title = memo,
      y = "Dimension 2",
      x = "Dimension 1"
    ) +
    ggplot2::scale_color_manual(values = color_palette)

  return(plot_t_SNE)
}

#' MDS plot
#'
#' Draw a multidimensional scaling (MDS) plot
#'
#'
#' @param data Matrix of gene data that has been through
#'  \code{\link{pre_process}()}
#' @param sample_info Matrix of sample information from experiment design file
#' @param selected_color String designating factor to color points by. Should be
#'  one of the design factors from the design file or "Names" as default which
#'  automatically detects groups from gene data file
#' @param selected_shape String designating factor to shape points by.
#'  Should be one of the design factors from the design file or "Names" as
#'  default which automatically detects groups from gene data file
#' @param plots_color_select Vector of colors for plots
#'
#' @export
#' @return A \code{ggplot} object formatted PCA plot
#'
#' @family PCA functions
#'
MDS_plot <- function(data,
                     sample_info,
                     selected_shape,
                     selected_color,
                     plots_color_select) {
  # no design file
  if (is.null(selected_color)) {
    selected_color <- "Names"
  }
  if (is.null(selected_shape)) {
    selected_shape <- "Names"
  }

  counts <- data
  memo <- ""

  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }

  x <- data
  y <- sample_info

  fit <- cmdscale(
    # dist_functions()$pearson_correlation(t(x)),
    dist_functions()$Euclidean(t(x)),
    eig = T,
    k = 2
  )
  pcaData <- as.data.frame(fit$points[, 1:2])

  # Missing design clause
  if (is.null(sample_info)) {
    pcaData <- cbind(pcaData, detect_groups(colnames(x), y))
  } else {
    pcaData <- cbind(pcaData, detect_groups(colnames(x), y), sample_info)
  }
  colnames(pcaData)[1:3] <- c("x1", "x2", "Names")
  
  pcaData$tooltip_text <- paste0(
    "Name: ", rownames(pcaData),
    "<br>Dimension 1: ", round(pcaData$x1, 2),
    "<br>Dimension 2: ", round(pcaData$x2, 2)
  )
  
  nshapes <- (length(unique(as.factor(pcaData$Names))) / 8) + 1

  # Set point & text size based on number of sample
  point_size <- 6

  if (ncol(x) >= 40) {
    point_size <- 3
    # text_size <- 16
  }

  # plot color scheme
  color_palette <- generate_colors(n = nlevels(as.factor(pcaData$Names)), palette_name = plots_color_select)

  p <- ggplot2::ggplot(
    data = pcaData,
    ggplot2::aes_string(
      x = "x1",
      y = "x2",
      color = selected_color,
      shape = selected_shape,
      text = "tooltip_text"
    )
  )
  p <- p + ggplot2::geom_point(size = point_size) +
    # Preferred shapes
    ggplot2::scale_shape_manual(values = rep(c(15, 16, 18, 0, 1, 3:5), nshapes)) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
      ),
      axis.text.y = ggplot2::element_text(
        size = 16
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      )
    ) +
    ggplot2::labs(
      title = memo,
      y = "Dimension 2",
      x = "Dimension 1"
    ) +
    ggplot2::scale_color_manual(values = color_palette)

  return(p)
}


#' Correlations Between Principle Components and Factors
#'
#' This function calculates the correlation between the principle components
#'  and factors and return the result in a formatted string.
#'
#' @param data Matrix of gene data that has been through
#'  \code{\link{pre_process}()}
#' @param sample_info  Matrix of sample information from experiment design file
#'
#' @export
#' @return A string
#'
#' @family PCA functions
#'
pc_factor_correlation <- function(data,
                                  sample_info) {
  x <- data
  y <- sample_info
  if (is.null(y)) {
    y <- as.matrix(detect_groups(colnames(data)))
  }

  if (dim(y)[2] == 1) {
    return(NULL)
  }
  pca.object <- prcomp(t(x))

  # 5 pc's or number of columns if <5
  npc <- min(5, ncol(data))

  pcaData <- as.data.frame(pca.object$x[, 1:npc])
  pvals <- matrix(1, nrow = npc, ncol = ncol(y))
  for (i in 1:npc) {
    for (j in 1:ncol(y)) {
      pvals[i, j] <- summary(
        aov(
          pcaData[, i] ~ as.factor(y[, j])
        )
      )[[1]][["Pr(>F)"]][1]
    }
  }
  pvals <- pvals * npc * ncol(y) # correcting for multiple testing
  pvals[pvals > 1] <- 1
  colnames(pvals) <- colnames(y)
  rownames(pvals) <- paste0("PC", 1:npc)
  a <- ""
  nchar0 <- nchar(a)
  for (i in 1:npc) {
    j <- which.min(pvals[i, ])
    if (pvals[i, j] < 0.05) {
      a <- paste0(
        a, rownames(pvals)[i],
        " is correlated with ", colnames(pvals)[j],
        " (p=", sprintf("%-3.2e", pvals[i, j]), ")."
      )
    }
  }
  return(a)
}



#' Principal Component Analysis with PCAtools package
#'
#' Draw a PCA plot using PCAtools package
#'
#' @param data Matrix of gene data that has been through
#'  \code{\link{pre_process}()}
#' @param sample_info Matrix of sample information from experiment design file
#' @param select_gene_id String indicating which gene id to use, default is
#'  "symbol"
#' @param all_gene_names Dataframe of gene names from
#'  \code{\link{get_all_gene_names}}
#' @param selected_x String indicating x axis selection, eg "PC1"
#' @param selected_y String indicating y axis selection, eg "PC2"
#' @param encircle TRUE/FALSE to draw shapes in plot, default is true
#' @param encircleFill TRUE/FALSE to fill shapes in plot, default is TRUE
#' @param showLoadings TRUE/FALSE to draw gene vectors onto plot, default is
#'  TRUE
#' @param pointlabs TRUE/FALSE to show column names on points, default is TRUE
#' @param point_size Positive number value to control point size
#' @param ui_color String designating factor to color points by.
#'  Should be one of the design factors from the design file
#' @param ui_shape String designating factor to shape points by.
#'  Should be one of the design factors from the design file
#'
#' @export
#' @return A \code{ggplot} object formatted as a PCA plot using PCAtools package
#'
#' @family PCA functions
#'
#' @seealso \code{\link[PCAtools]{pca}()} and \code{\link[PCAtools]{biplot}()}
#'  for the original functions from the PCAtools package
#'
PCA_biplot <- function(data,
                       sample_info,
                       select_gene_id = "symbol",
                       all_gene_names,
                       selected_x = "PC1",
                       selected_y = "PC2",
                       encircle = TRUE,
                       encircleFill = TRUE,
                       showLoadings = TRUE,
                       pointlabs = TRUE,
                       point_size = 4.0,
                       ui_color = NULL,
                       ui_shape = NULL) {
  # missing design
  if (is.null(sample_info)) {
    meta_data <- as.data.frame(colnames(data))
    rownames(meta_data) <- colnames(data)
  } else {
    meta_data <- sample_info
  }

  # Swap rownames
  data <- rowname_id_swap(
    data_matrix = data,
    all_gene_names = all_gene_names,
    select_gene_id = select_gene_id
  )

  pca_obj <- PCAtools::pca(data, metadata = meta_data, removeVar = 0.1)

  if (pointlabs == TRUE) {
    show_point_labels <- rownames(pca_obj$metadata)
  } else {
    show_point_labels <- NULL
  }

  PCAtools::biplot(
    pcaobj = pca_obj,
    x = selected_x,
    y = selected_y,
    colby = ui_color,
    shape = ui_shape,
    # colLegendTitle = 'Color?',
    encircle = encircle,
    encircleFill = encircleFill,
    showLoadings = showLoadings,
    lab = show_point_labels,
    legendPosition = "right",
    legendLabSize = 16,
    legendIconSize = 8.0,
    pointSize = point_size,
    title = "Principal Component Scores"
  )
}

#' Principal Component Analysis with PCAtools package
#'
#' Draw a Scree plot with Horn's and Elbow suggestion for cutoffs using
#' PCAtools package
#'
#' @param data Matrix of gene data that has been through
#'  \code{\link{pre_process}()}
#'
#' @export
#' @return A \code{ggplot} object formatted as a Scree plot
#'
#' @family PCA functions
#'
#' @seealso \code{\link[PCAtools]{screeplot}()} for PCAtools documentation
#'
PCA_Scree <- function(processed_data) {
  suppressWarnings(
    pca_obj <- PCAtools::pca(mat = processed_data, removeVar = 0.1)
  )
  suppressWarnings(
    horn <- PCAtools::parallelPCA(processed_data)
  )

  suppressWarnings(
    elbow <- PCAtools::findElbowPoint(pca_obj$variance)
  )

  p <- PCAtools::screeplot(
    pca_obj,
    vline = c(horn$n, elbow)
  )
  p <- p +
    ggplot2::geom_label(
      ggplot2::aes(
        x = horn$n + .1,
        y = 60,
        label = "Horn's",
        vjust = .5,
        hjust = .5,
        size = 8
      )
    )
  p <- p + ggplot2::geom_label(ggplot2::aes(
    x = elbow + .1,
    y = 70,
    label = "Elbow",
    vjust = .5,
    hjust = .5,
    size = 8
  ))
  return(p)
}

#' Principal Component Analysis with PCAtools package
#'
#' Generates a plot showing correlations between Principal Components and design factors
#'
#' @param data Matrix of gene data that has been through
#'  \code{\link{pre_process}())}
#' @param sample_info Matrix of sample information from experiment design file
#'
#' @return A \code{ggplot} object as a formatted plot generated with PCAtools
#'  package
#' @export
#'
#' @family PCA functions
#' @seealso \code{\link[PCAtools]{pca}()} and
#'   \code{\link[PCAtools]{eigencorplot}()} for specific documentation from the
#'   PCAtools package
#'
PCAtools_eigencorplot <- function(processed_data,
                                  sample_info) {
  # missing design
  if (is.null(sample_info)) {
    return(NULL)
    # meta_data <- as.data.frame(colnames(processed_data))
    # colnames(meta_data)[1] <- "Sample_Name"
  } else {
    meta_data <- sample_info
    
    # Design Factors must be converted to numeric
    meta_data <- as.data.frame(meta_data)
    meta_data <- sapply(meta_data, function(x) as.numeric(factor(x)))
    meta_data <- as.data.frame(meta_data)

    # maintain rownames
    rownames(meta_data) <- rownames(sample_info)

    # create PCA object
    pca_obj <- PCAtools::pca(processed_data, metadata = meta_data, removeVar = 0.1)

    # plot
    p <- PCAtools::eigencorplot(pca_obj, metavars = colnames(meta_data))
    return(p)
  }
}

#'  Principal Component Analysis with PCAtools package
#'  
#'  Convert PCA loadings to percent contribution values by squaring and plot
#'  them
#'
#' @param processed_data Matrix of gene data that has been through
#'  \code{\link{pre_process}())}
#' @param all_gene_names vector of gene names from 
#' \code{\link{pre_process}())}
#' @param PC # Selected principal component
#'
#' @return A \code{ggplot} object as a formatted plot generated with PCAtools
#'  package
#'  
#' @export
#' 
#' @family PCA Functions
#'
var_imp_plots <- function(data,
                          all_gene_names,
                          PC){
  
  # Swap rownames with gene symbols
  data <- rowname_id_swap(
    data_matrix = data,
    all_gene_names = all_gene_names,
    select_gene_id = "symbol"
  )
  
  #create pca object
  pca_obj <- PCAtools::pca(mat = data)
  
  #Convert loadings to importance values
  loadings <- as.data.frame(PCAtools::getLoadings(pca_obj))
  imp_vals <- data.frame(vals = ((loadings[[paste0(PC)]])^2) * 100, 
                         gene = rownames(loadings))
  
  #Sort and trim importance values greatest to least
  imp_vals <- imp_vals[order(imp_vals$vals, decreasing = TRUE), ]
  imp_vals <- imp_vals[c(1:10),]
  
  imp_plot <- ggplot2::ggplot(data = imp_vals)+
    ggplot2::aes(
      x = reorder(gene, -vals), 
      y = vals, 
      label = paste(round(vals, 2), "%")
    )+
    ggplot2::geom_bar(stat = "identity", color = "black", fill = "dodgerblue")+
    ggplot2::geom_text(vjust = 2, size = 5)+
    ggplot2::theme_minimal()+
    ggplot2::labs(
      x = "Genes", 
      y = "Contribution (%)", 
      title = paste0("Contribution of Genes to ", PC, " (Top 10)")
    )+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       size = 12,
                                                       vjust = 0.5),
                   axis.text.y = ggplot2::element_text(size = 12),
                   plot.title = ggplot2::element_text(size = 20),
                   axis.title.x = ggplot2::element_text(size = 12,
                                                        vjust = 0.8),
                   axis.title.y = ggplot2::element_text(size = 12)
    )
  
  return(imp_plot)
  
}


#' Gets plot width dimensions from session$clientdata
#'
#' @param client_data session info from shiny
#' @param plot_name Name of plot object
#' @param tab Tab name
#'
#' @return string with dimension info
#' @export
#'
get_plot_width <- function(client_data,
                           plot_name,
                           tab) {
  tryCatch(
    {
      width <- paste0("output_", tab, "-", plot_name, "_width")
      height <- paste0("output_", tab, "-", plot_name, "_height")
      w <- client_data[[width]]
      # h <- client_data[[height]]
    },
    Warning = function(w) {
      print("Warning in get_plot_width")
    },
    error = function(e) {
      print("Error in get_plot_width")
    }
  )

  # return width in inches -- default to 6.5
  return(6.5)
}


#' Gets plot height size from session$clientdata
#'
#' @param client_data session info from shiny
#' @param plot_name Name of plot object
#' @param tab Tab name
#'
#' @return height in inches of plot
#' @export
#'
get_plot_height <- function(client_data,
                            plot_name,
                            tab) {
  tryCatch(
    {
      h <- 5.0099
      width <- paste0("output_", tab, "-", plot_name, "_width")
      height <- paste0("output_", tab, "-", plot_name, "_height")
      w <- client_data[[width]]
      h <- client_data[[height]]
      aspect_ratio <- h / w
      h <- aspect_ratio * 6.5
    },
    warning = function(w) {
      print("Warning in get_plot_height")
    },
    error = function(e) {
      print("Error in get_plot_height")
    }
  )


  return(round(h, 3))
}
