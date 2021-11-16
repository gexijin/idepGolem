#' fct_08_bicluster.R This file holds all of the main data analysis functions
#' associated with eighth tab of the iDEP website.
#'
#'
#' @section fct_07_bicluster.R functions:
#'
#'
#' @name fct_08_bicluster.R
NULL

#' WGCNA DATA
get_wgcna <- function(
  data,
  n_genes,
  soft_power,
  min_module_size
) {
  max_gene_wgcna = 3000
  # http://pklab.med.harvard.edu/scw2014/WGCNA.html
  if(n_genes > dim(data)[1]) {
    # Max	as data
    n_genes <- dim(data)[1]
  }
  if(n_genes < 50) {
    return(NULL)
  }
  if(dim(data)[2] < 4) {
    return(NULL)
  }
  if(n_genes > max_gene_wgcna) {
    n <- max_gene_wgcna
  } 			

  dat_expr <- t(data[1:n_genes, ])
  sub_gene_names <- colnames(dat_expr)

  # Choosing a soft-threshold to fit a scale-free topology to the network
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  sft <- WGCNA::pickSoftThreshold(
    dat_expr,
    dataIsExpr = TRUE,
    powerVector = powers,
    corFnc = cor,
    corOptions = list(use = "p"),
    networkType = "unsigned"
  )

  # Calclute the adjacency matrix
  adj <- WGCNA::adjacency(
    dat_expr,
    type = "unsigned",
    power = soft_power
  )

  # Turn adjacency matrix into topological overlap to
  # minimize the effects of noise and spurious associations
  tom <- WGCNA::TOMsimilarityFromExpr(
    dat_expr,
    networkType = "unsigned",
    TOMType = "unsigned",
    power = soft_power
  )
  colnames(tom) <- sub_gene_names
  rownames(tom) <- sub_gene_names

  # Module detection
  gene_tree <- flashClust::flashClust(
    as.dist(1 - tom),
    method = "average"
  )

  # Module identification using dynamic tree cut
  dynamic_mods <- dynamicTreeCut::cutreeDynamic(
    dendro = gene_tree,
    method = "tree",
    minClusterSize = min_module_size
  )

  dynamic_colors <- WGCNA::labels2colors(dynamic_mods)
  module_info <- cbind(sub_gene_names, dynamic_colors, dynamic_mods)
  # Remove genes not in any modules
  module_info <- module_info[which(module_info[, 2] != "grey"), ]
  module_info <- module_info[order(module_info[, 3]), ]
  n_modules <- length(unique(dynamic_colors)) - 1
  n_genes <- dim(module_info)[1]	
			
  return(list(
    data = t(dat_expr),
    powers = powers,
    sft = sft,
    tom = tom,
    dynamic_colors = dynamic_colors,
    module_info = module_info,
    n_modules = n_modules,
    n_genes = n_genes
  ))
}

#' MODULE NETWORK PLOT
get_module_plot <- function(
  wgcna
) {
  diss <- 1 - wgcna$tom
  dynamic_colors <- wgcna$dynamic_colors
		
  hier <- flashCut::flashClust(as.dist(diss), method = "average")

  # Set the diagonal of the dissimilarity to NA 
  diag(diss) <- NA
  WGCNA::plotDendroAndColors(
    hier,
    dynamic_colors,
    "Dynamic Tree Cut",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = "Gene dendrogram and module colors"
  )
}