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
  
  
  pca.object <- prcomp(t(x))
  npc=5
  pcaData = as.data.frame(pca.object$x[,1:npc]); 
  pvals = matrix(1,nrow=npc,ncol=ncol(y))
  for (i in 1:npc ){
    for (j in 1:ncol(y) )
      pvals[i,j] =summary( aov(pcaData[,i] ~ as.factor(y[,j])))[[1]][["Pr(>F)"]][1]
  }
  
  pvals = pvals * npc* ncol(y)   # correcting for multiple testing
  pvals[pvals>1] = 1
  
  colnames(pvals) = colnames(y)
  rownames(pvals) = paste0("PC",1:npc)
  
  
  
  a="<h4>Correlation between Principal Components (PCs) with factors </h4>"
  nchar0 = nchar(a)
  for ( i in 1:npc){
    j = which.min(pvals[i,])
    if(pvals[i,j]< 0.05) a=paste0(a,rownames(pvals)[i], 
                                  " is correlated with ", colnames(pvals)[j],
                                  " (p=",   sprintf("%-3.2e",pvals[i,j]),").<br>")
  }
  if(nchar(a) == nchar0 ) return(NULL) else 
    
  groups <- detect_groups(sample_names=colnames(data), sample_info=sample_info)

  P1 <- paste("PC", PCAx, sep="")
  P2 <- paste("PC", PCAy, sep="")
  
  
  
 # p <- ggplot2::ggplot(pcaData, ggplot2::aes(PC1, PC2, color=groups,shape=groups)) 
  p <- ggplot2::ggplot(pcaData, ggplot2::aes(x=pcaData[, as.integer(PCAx)],y = pcaData[, as.integer(PCAy)], color=groups,shape=groups)) 
  
  p <- p + ggplot2::geom_point()

  PCAxy <- c(as.integer(PCAx ),as.integer(PCAy) ) # selected principal components
  percentVar=round(100*summary(pca.object)$importance[2, PCAxy],0)
  
  p=p+ggplot2::xlab(paste0("PC", PCAx ,": ", percentVar[1],"% variance")) 
  p=p+ggplot2::ylab(paste0("PC", PCAy ,": ",percentVar[2],"% variance")) 
  p=p+ggplot2::ggtitle("Principal Component Analysis (PCA)")+ggplot2::coord_fixed(ratio=1.0)+ 
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16,hjust = 0.5)) + ggplot2::theme(aspect.ratio=1) +
    ggplot2::theme(axis.text.x = ggplot2::element_text( size = 16),
                   axis.text.y = ggplot2::element_text( size = 16),
                   axis.title.x = ggplot2::element_text( size = 16),
                   axis.title.y = ggplot2::element_text( size = 16) ) +
    ggplot2::theme(legend.text=ggplot2::element_text(size=16))
  print(p)
  
  
  return(p)
}






#PCA_plot(processed_data$data, loaded_data$sample_info)


t_SNE_plot <- function(data, sample_info){
  
  x <- data
  y <- sample_info		 
  tsne <- Rtsne::Rtsne(t(x), dims = 2, perplexity=1, verbose=FALSE, max_iter = 400)
  
  pcaData = as.data.frame(tsne$Y);
  pcaData = cbind(pcaData,detect_groups(colnames(x), y) )
  
  colnames(pcaData) = c("x1", "x2", "Sample_Name")
  
  
  p <- ggplot2::ggplot(pcaData, ggplot2::aes(x1, x2, color=Sample_Name, shape = Sample_Name)) 
  
  p<-	p+ggplot2::geom_point()
  
  p=p+ggplot2::xlab("Dimension 1") 
  p=p+ggplot2::ylab("Dimension 2") 
  p=p+ggplot2::ggtitle("t-SNE plot")+ ggplot2::coord_fixed(ratio=1.)+ 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::theme(aspect.ratio=1) +
    ggplot2::theme(axis.text.x = ggplot2::element_text( size = 16),
                   axis.text.y = ggplot2::element_text( size = 16),
                   axis.title.x = ggplot2::element_text( size = 16),
                   axis.title.y = ggplot2::element_text( size = 16) ) +
    ggplot2::theme(legend.text=ggplot2::element_text(size=16))
  
  return(p)
}





# MDS
MDS_plot <- function(
  data,
  sample_info
){
  x <- data
  y <- sample_info
  
  
  #fit <- cmdscale( dist_pcc(t(x) ), eig=T, k=2)
  fit <- cmdscale( dist_functions()$pearson_correlation(t(x) ), eig=T, k=2)
  

  pcaData = as.data.frame(fit$points[,1:2])
  pcaData = cbind(pcaData,detect_groups(colnames(x), y) )
  
  colnames(pcaData) = c("x1", "x2", "Sample_Name")		 	
  
  
  p <- ggplot2::ggplot(pcaData, ggplot2::aes(x1, x2, color=Sample_Name, shape = Sample_Name))  
  
  p <- p+ggplot2::geom_point()
  
  #	p <- p+	 scale_shape_manual(values= shapes)	 
  
  p=p+ggplot2::xlab("Dimension 1") 
  p=p+ggplot2::ylab("Dimension 2") 
  p=p+ggplot2::ggtitle("Multidimensional scaling (MDS)")+ ggplot2::coord_fixed(ratio=1.)+ 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::theme(aspect.ratio=1) +
    ggplot2::theme(axis.text.x = ggplot2::element_text( size = 16),
                   axis.text.y = ggplot2::element_text( size = 16),
                   axis.title.x = ggplot2::element_text( size = 16),
                   axis.title.y = ggplot2::element_text( size = 16) ) +
    ggplot2::theme(legend.text=ggplot2::element_text(size=16))
 
  return(p)
  
}




#Pathway





# 
# ## edit later
# # Runs pathway analysis using PGSEA; this is copied and revised from PGSEA package
# my_pgsea <- function(exprs, cl, range = c(25, 500), ref = NULL, center = TRUE,
#                      p.value = 0.005, weighted = TRUE, nPermutation = 100, enforceRange = TRUE, ...) {
#   if (is(exprs, "ExpressionSet")) {
#     exprs <- exprs(exprs)
#   }
#   if (!is.list(cl)) {
#     stop("cl need to be a list")
#   }
#   if (!is.null(ref)) {
#     if (!is.numeric(ref)) {
#       stop("column index's required")
#     }
#   }
#   if (!is.null(ref)) {
#     if (options()$verbose) {
#       cat("Creating ratios...", "\n")
#     }
#     ref_mean <- apply(exprs[, ref], 1, mean, na.rm = TRUE)
#     exprs <- sweep(exprs, 1, ref_mean, "-")
#   }
#   if (center) {
#     exprs <- scale(exprs, scale = FALSE)
#   } # column centering is done
#   results <- matrix(NA, length(cl), ncol(exprs))
#   rownames(results) <- names(cl)
#   colnames(results) <- colnames(exprs)
#   mode(results) <- "numeric"
#   Setsize <- c(rep(0, length(cl))) # gene set size vector
#   mean2 <- c(rep(0, length(cl))) # mean of the range of means
#   meanSD <- c(rep(0, length(cl))) # SD of the range of means
#   if (is.logical(p.value)) {
#     p.results <- results
#     mean.results <- results
#   }
#   for (i in 1:length(cl)) { # for each gene list
#     # cat("\nProcessing gene set",i);
#     if (class(cl[[i]]) == "smc") {
#       clids <- cl[[i]]@ids
#     } else if (class(cl[[i]]) %in% c("GeneColorSet", "GeneSet")) {
#       clids <- cl[[i]]@geneIds
#     } else {
#       clids <- cl[[i]]
#     }
#     if (options()$verbose) {
#       cat("Testing region ", i, "\n")
#     }
#     ix <- match(clids, rownames(exprs))
#     ix <- unique(ix[!is.na(ix)])
#     present <- sum(!is.na(ix))
#     Setsize[i] <- present
#     if (present < range[1]) {
#       if (options()$verbose) {
#         cat(
#           "Skipping region ", i, " because too small-",
#           present, ",\n"
#         )
#       }
#       next
#     }
#     if (present > range[2]) {
#       if (options()$verbose) {
#         cat(
#           "Skipping region ", i, " because too large-",
#           present, "\n"
#         )
#       }
#       next
#     }
#     texprs <- exprs[ix, ] # expression matrix for genes in gene set
#     if (any(is.na(texprs))) {
#       cat("Warning - 'NA' values within expression data, enrichment scores are estimates only.\n")
#     }
#     if (!is.matrix(texprs)) {
#       texprs <- as.matrix(texprs)
#     }
# 
#     stat <- try(apply(texprs, 2, t.test, ...))
#     means <- try(apply(texprs, 2, mean, trim = 0.1)) # trim mean
#     ps <- unlist(lapply(stat, function(x) x$p.value))
#     stat <- unlist(lapply(stat, function(x) x$statistic))
#     p.results[i, ] <- ps
#     mean.results[i, ] <- means
#     results[i, ] <- as.numeric(stat)
# 
#     # permutation of gene sets of the same size
#     if (nPermutation > 2) { # no permutation if <=2
#       meansRanges <- c(0, rep(nPermutation))
#       for (k in 1:nPermutation) {
#         ix <- sample.int(dim(exprs)[1], length(ix))
#         texprs <- exprs[ix, ]
#         means <- try(apply(texprs, 2, mean, trim = 0.1)) # trim mean
#         meansRanges[k] <- dynamicRange(means)
#       }
#       mean2[i] <- mean(meansRanges)
#       meanSD[i] <- sd(meansRanges, na.rm = TRUE) # NA are removed before calculating standard deviation
#     }
#   }
#   return(list(results = results, p.results = p.results, means = mean.results, size = Setsize, mean2 = mean2, meanSD = meanSD))
# }

