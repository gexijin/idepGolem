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
){
  
  ##x<-processed_data$data
  ##y<-loaded_data$sample_info
  
  x <- data
  y <- sample_info
  
  
  pca_object <- prcomp(t(x))
  npc <- 5
  pca_data <- as.data.frame(pca_object$x[, 1:npc]) 
  pvals <- matrix(1, nrow = npc, ncol = ncol(y))
  for (i in 1:npc) {
    for (j in 1:ncol(y)) {
      pvals[i, j] <- summary(aov(pca_data[, i] ~ as.factor(y[, j])))[[1]][["Pr(>F)"]][1]
    }
  }
  
  # Correcting for multiple testing
  pvals <- pvals * npc * ncol(y)
  pvals[pvals > 1] <- 1
  
  colnames(pvals) <- colnames(y)
  rownames(pvals) <- paste0("PC", 1:npc)
  
  
  
  a <- "<h4> Correlation between Principal Components (PCs) with factors </h4>"
  nchar_0 <- nchar(a)
  for (i in 1:npc) {
    j <- which.min(pvals[i,])
    if (pvals[i,j]< 0.05) {
      a <- paste0(
        a,
        rownames(pvals)[i], 
        " is correlated with ",
        colnames(pvals)[j],
        " (p=",
        sprintf("%-3.2e", pvals[i, j]), 
        ").<br>"
      )
    }
  }
  if(nchar(a) == nchar_0) {
    return(NULL)
  } else {
    groups <- detect_groups(sample_names = colnames(data), sample_info = sample_info)

    P1 <- paste("PC", PCAx, sep = "")
    P2 <- paste("PC", PCAy, sep = "")
  
    p <- ggplot2::ggplot(
      pca_data,
      ggplot2::aes(
        x = pca_data[, as.integer(PCAx)],
        y = pca_data[, as.integer(PCAy)],
        color = groups,
        shape = groups
      )
    ) + ggplot2::geom_point()

    # Selected principal components
    PCA_xy <- c(as.integer(PCAx), as.integer(PCAy))
    percent_var <- round(100 * summary(pca_object)$importance[2, PCA_xy], 0)
  
    p <- p + ggplot2::xlab(paste0("PC", PCAx, ": ", percent_var[1], "% variance")) +
      ggplot2::ylab(paste0("PC", PCAy, ": ", percent_var[2], "% variance")) + 
      ggplot2::ggtitle("Principal Component Analysis (PCA)") + 
      ggplot2::coord_fixed(ratio = 1.0) + 
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16,hjust = 0.5),
        aspect.ratio = 1,
        axis.text.x = ggplot2::element_text(size = 16),
        axis.text.y = ggplot2::element_text(size = 16),
        axis.title.x = ggplot2::element_text(size = 16),
        axis.title.y = ggplot2::element_text(size = 16),
        legend.text = ggplot2::element_text(size=16)
      )
  
    return(p)
  }
}

#' TSNE FUNCTION 
t_SNE_plot <- function(
  data,
  sample_info
) {	 
  tsne <- Rtsne::Rtsne(
    t(data),
    dims = 2,
    perplexity = 1,
    verbose = FALSE,
    max_iter = 400
  )
  
  pca_data <- as.data.frame(tsne$Y);
  pca_data <- cbind(pca_data, detect_groups(colnames(data), sample_info))
  
  colnames(pca_data) <- c("x1", "x2", "Sample_Name")
  
  p <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(x1, x2, color = Sample_Name, shape = Sample_Name)
  ) + 
    ggplot2::geom_point() +
    ggplot2::xlab("Dimension 1") + 
    ggplot2::ylab("Dimension 2") +
    ggplot2::ggtitle("t-SNE plot") +
    ggplot2::coord_fixed(ratio=1.) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_text(size = 16),
      axis.text.y = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_text(size = 16),
      axis.title.y = ggplot2::element_text(size = 16),
      legend.text = ggplot2::element_text(size = 16)
    )
  
  return(p)
}

# MDS FUNCTION
MDS_plot <- function(
  data,
  sample_info
){
  fit <- cmdscale(
    dist_functions()$pearson_correlation(t(data)),
    eig = T,
    k = 2
  )
  
  pca_data <- as.data.frame(fit$points[, 1:2])
  pca_data <- cbind(pca_data, detect_groups(colnames(data), sample_info))
  
  colnames(pca_data) <- c("x1", "x2", "Sample_Name")		 	
  
  p <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(
      x1,
      x2,
      color = Sample_Name,
      shape = Sample_Name
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::xlab("Dimension 1") +
    ggplot2::ylab("Dimension 2") +
    ggplot2::ggtitle("Multidimensional scaling (MDS)") +
    ggplot2::coord_fixed(ratio = 1.) + 
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_text(size = 16),
      axis.text.y = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_text(size = 16),
      axis.title.y = ggplot2::element_text(size = 16),
      legend.text = ggplot2::element_text(size = 16)
    )
 
  return(p)
}
