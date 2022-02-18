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



#' Principal Component Analysis
#'
#' Draw a PCA plot where user selects which PCs on axes
#'
#' @param data Data that has been through pre-processing
#' @param sample_info Matrix array with experiment info
#' @param PCAx PC on x axis
#' @param PCAy PC on y axis
#'
#' @return Formatted PCA plot 
#' 
PCA_plot <- function(
  data,
  sample_info,
  PCAx = 1,
  PCAy = 2
) {
  counts <- data
  memo <- ""
  
  if (ncol(counts) > 100) {
    part <- 1:100
    counts <- counts[, part]
    memo <- paste("(only showing 100 samples)")
  }
  
  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }
  
  ## x<-processed_data$data
  ## y<-loaded_data$sample_info
  x <- data
  y <- sample_info
  pca.object <- prcomp(t(x))
  npc <- 5
  pcaData <- as.data.frame(pca.object$x[, 1:npc])
  # pvals <- matrix(1, nrow = npc, ncol = ncol(y))
  # for (i in 1:npc) {
  #   for (j in 1:ncol(y)) {
  #     pvals[i, j] <- summary(aov(pcaData[, i] ~ as.factor(y[, j])))[[1]][["Pr(>F)"]][1]
  #   }
  # }
  # pvals <- pvals * npc * ncol(y) # correcting for multiple testing
  # pvals[pvals > 1] <- 1
  # colnames(pvals) <- colnames(y)
  # rownames(pvals) <- paste0("PC", 1:npc)
  # a <- "<h4>Correlation between Principal Components (PCs) with factors </h4>"
  # nchar0 <- nchar(a)
  # for (i in 1:npc) {
  #   j <- which.min(pvals[i, ])
  #   if (pvals[i, j] < 0.05) {
  #     a <- paste0(
  #       a, rownames(pvals)[i],
  #       " is correlated with ", colnames(pvals)[j],
  #       " (p=", sprintf("%-3.2e", pvals[i, j]), ").<br>"
  #     )
  #   }
  # }
  # if (nchar(a) == nchar0) {
  #   return(NULL)
  # } else {
    groups <- detect_groups(sample_names = colnames(data), sample_info = sample_info)
  
  
  if (nlevels(groups) <= 1 | nlevels(groups) > 20) {
    group_fill <- NULL
    legend <- "none"
  } else {
    group_fill <- groups
    legend <- "right"
  }
  
  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }
  
  P1 <- paste("PC", PCAx, sep = "")
  P2 <- paste("PC", PCAy, sep = "")
  
  # Set point & text size based on number of sample
  point_size <- 6
  if (ncol(x) >= 40) {
    point_size <- 3
  }
  
  plot_PCA <- ggplot2::ggplot(
    data = pcaData, 
    ggplot2::aes(
      x = pcaData[, as.integer(PCAx)],
      y = pcaData[, as.integer(PCAy)],
      color = groups,
      shape = groups
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
      title = paste("Principal Component Analysis (PCA) ", memo),
      y = "Dimension 2",
      x = "Dimension 1"
    )
  # selected principal components
  PCAxy <- c(as.integer(PCAx), as.integer(PCAy))
  percentVar <- round(100 * summary(pca.object)$importance[2, PCAxy], 0)
  plot_PCA <- plot_PCA + ggplot2::xlab(paste0("PC", PCAx, ": ", percentVar[1], "% Variance"))
  plot_PCA <- plot_PCA + ggplot2::ylab(paste0("PC", PCAy, ": ", percentVar[2], "% Variance"))
  return(plot_PCA)
}

#' TSNE FUNCTION 
#'
#' Draw a t-sne plot where user selects which PCs on axes
#'
#'
#' @param data Data that has been through pre-processing
#' @param sample_info Matrix array with experiment info
#'
#' @return Formatted T-sne plot
#'
#'
t_SNE_plot <- function(
  data,
  sample_info
) {
  counts <- data
  memo <- ""
  
  if (ncol(counts) > 100) {
    part <- 1:100
    counts <- counts[, part]
    memo <- paste("(only showing 100 samples)")
  }

  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }

  
  x <- data
  y <- sample_info
  tsne <- Rtsne::Rtsne(t(x), dims = 2, perplexity = 1, verbose = FALSE, max_iter = 400)
  pcaData <- as.data.frame(tsne$Y)
  pcaData <- cbind(pcaData, detect_groups(colnames(x), y))

  colnames(pcaData) <- c("x1", "x2", "Sample_Name")
  
  # Set point size based on number of sample
  point_size <- 6
  if (ncol(x) >= 40) {
    point_size <- 3
  }
  
  #Generate plot
  plot_t_SNE <- ggplot2::ggplot(
    data = pcaData,
    ggplot2::aes(
      x = x1,
      y = x2,
      color = Sample_Name,
      shape = Sample_Name)
    ) +
  ggplot2::geom_point(size = point_size) +
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
      title = paste("T-SNE ", memo),
      y = "Dimension 2",
      x = "Dimension 1"
    )

  return(plot_t_SNE)
}

#' MDS FUNCTION
#'
#' Draw a MDS plot
#'
#'
#' @param data Data that has been through pre-processing
#' @param sample_info Matrix array with experiment info
#'
#' @return Formatted PCA plot
#'
MDS_plot <- function(
  data,
  sample_info
) {
  counts <- data
  memo <- ""
  
  if (ncol(counts) > 100) {
    part <- 1:100
    counts <- counts[, part]
    memo <- paste("(only showing 100 samples)")
  }
  
  if (ncol(counts) < 31) {
    x_axis_labels <- 16
  } else {
    x_axis_labels <- 12
  }
  
  x <- data
  y <- sample_info

  fit <- cmdscale(
    dist_functions()$pearson_correlation(t(x)),
    eig = T,
    k = 2
  )
  pcaData <- as.data.frame(fit$points[, 1:2])
  pcaData <- cbind(pcaData, detect_groups(colnames(x), y))
  colnames(pcaData) <- c("x1", "x2", "Sample_Name")
  
  # Set point & text size based on number of sample
  point_size <- 6

  if (ncol(x) >= 40) {
    point_size <- 3
    #text_size <- 16
  }
  p <- ggplot2::ggplot(
    data = pcaData,
    ggplot2::aes(
      x = x1,
      y = x2,
      color = Sample_Name,
      shape = Sample_Name
    )
  )
  p <- p + ggplot2::geom_point(size = point_size) +
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
      title = paste("Multi-Dimensional Scaling (MDS) ", memo),
      y = "Dimension 2",
      x = "Dimension 1"
    )
  return(p)
}


#' Correlations Between Principle Components and Factors
#'
#' Desc
#'
#'
#' @param data Data that has been through pre-processing
#' @param sample_info Matrix array with experiment info
#'
#' @return text with correlation
#'
pc_factor_correlation <- function(
  data,
  sample_info
) {
  x <- data
  y <- sample_info
  if (is.null(y)){
    y <- as.matrix(detect_groups(colnames(data)))
  }
  
  if(dim(y)[2] == 1)
  {
    return("No design file uploaded")
  }
  pca.object <- prcomp(t(x))
  npc <- 5
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
  a <- "Correlation between Principal Components (PCs) with factors: "
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