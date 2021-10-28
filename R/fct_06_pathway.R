#' fct_08_pathway.R This file holds all of the main data analysis functions
#' associated with eighth tab of the iDEP website.
#'
#'
#' @section fct_08_pathway.R functions:
#'
#'
#' @name fct_08_pathway.R
NULL

#' GAGE PATHWAY DATA
gage_data <- function(
  select_go,
  select_contrast,
  min_set_size,
  max_set_size,
  limma,
  gene_p_val_cutoff,
  gene_sets,
  absolute_fold,
  pathway_p_val_cutoff,
  n_pathway_show
) {
  if(select_go == "ID not recognized!") {
    return(as.data.frame("Gene ID not recognized."))
  }
  if(is.null(select_contrast)) {
    return(NULL)
  }
  my_range <- c(min_set_size, max_set_size)
  no_sig <- as.data.frame("No significant pathway found.")
  if(length(limma$top_genes) == 0 ) {
    return(no_sig)
  }
  
  if(length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]  
  } else {
    top <- limma$top_genes
    ix <- match(select_contrast, names(top))
    if(is.na(ix)) {
      return(no_sig)
    }
    top_1 <- top[[ix]] 
  }
  if(dim(top_1)[1] == 0) {
    return(no_sig)
  }
  colnames(top_1) <- c("Fold", "FDR")

  top_1 <- top_1[which(top_1$FDR < gene_p_val_cutoff), ]
 
  gmt <- gene_sets 
  if(length(gene_sets) == 0) {
    return(as.data.frame("No gene set found!"))
  }

  fold <- top_1[, 1]
  names(fold) <- rownames(top_1)
  if(absolute_fold) {
    fold <- abs(fold)
  }
  paths <- gage::gage(fold, gsets = gmt, ref = NULL, samp = NULL)

  paths <-  rbind(paths$greater, paths$less)
 
  if(dim(paths)[1] < 1 | dim(paths)[2] < 6 ) {
    return(no_sig)
  }
  top_1 <- paths[, c("stat.mean", "set.size", "q.val")]
  colnames(top_1) <- c("statistic", "Genes", "adj.Pval")
  top_1 <- top_1[order(top_1[, 3]), ]  
  if(length(which(top_1[, 3] <= pathway_p_val_cutoff)) == 0) {
    return(no_sig)
  }
  top_1 <- top_1[which(top_1[, 3] <= pathway_p_val_cutoff), , drop = FALSE]
  if(dim(top_1)[1] > n_pathway_show) {
    top1 <- top1[1:input$nPathwayShow, ,drop=FALSE]
  }		 
  top_1 <- as.data.frame(top_1)
  top_1 <- cbind(rep(select_contrast, dim(top_1)[1]), row.names(top_1), top_1) 
  top_1$statistic <- as.character(round(as.numeric(top_1$statistic), 4)) 
  top_1$adj.Pval <- sprintf("%-2.1e", as.numeric(top_1$adj.Pval))
  top_1[, 2] <- as.character(top_1[, 2])
  top_1[, 1] <- as.character(top_1[, 1])
  colnames(top_1)[1] <- "Direction"
  if(pathway_method == 1) {
    p.m <- "GAGE"
  } else if(pathway_method == 2) {
    p.m <- "PGSEA"
  } else if(pathway_method == 3) {
    p.m <- "GSEA"
  } else if(pathway_method == 4) {
    p.m <- "PGSEA_All"
  } else if(pathway_method == 5) {
    p.m <- "ReactomePA"
  }
  colnames(top_1)[2] <- paste(p.m, "analysis:", gsub("-", " vs ", select_contrast))
  top_1[which(top_1[, 3] > 0), 1] <- "Up"
  top_1[which(top_1[, 3] < 0), 1] <- "Down"
  top_1 <- top_1[order(top_1[, 1], -abs(as.numeric(top_1[, 3]))), ]
  top_1[duplicated(top_1[, 1]), 1] <- ""

  return(top_1)
}