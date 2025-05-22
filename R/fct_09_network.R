#' fct_08_bicluster.R This file holds all of the main data analysis functions
#' associated with eighth tab of the iDEP website.
#'
#'
#' @section fct_07_bicluster.R functions:
#'
#'
#' @name fct_08_bicluster.R
NULL

#' Utilize WGCNA function
#'
#' Run the WGCNA package on the processed data.
#'
#' @param data Data matrix that has been through the pre-processing
#' @param n_genes Number of most variable genes to use in the function
#' @param soft_power Value between 1-20
#'  (https://support.bioconductor.org/p/87024/)
#' @param min_module_size For modules detected by WGCNA, set a minimum size
#'
#' @export
#' @return A list of 8 objects. \code{data} is a submatrix of the input
#'  parameter \code{data} which on contains the genes that were selected with
#'  \code{n_genes}. \code{powers} is a numeric vector from 1-10 and then the
#'  even numbers from 10-20. \code{sft} is a return from the WGCNA package that
#'  is a list containing \code{powerEstimate} and \code{fitIndices}. For
#'  information on these objects visit
#'  https://www.rdocumentation.org/packages/WGCNA/versions/1.70-3/topics/pickSoftThreshold.
#'  \code{tom} is another return from the WGCNA package, for information on this
#'  matrix visit https://www.rdocumentation.org/packages/WGCNA/versions/1.70-3/topics/TOMsimilarityFromExpr.
#'  \code{dynamic_colors} is a vector of colors created by the WGCNA package to
#'  correspond to the modules that were identified by the pacakage.
#'  \code{module_info} is a data frame with a gene name, a color, and the
#'  identified module for the gene. \code{n_modules} gives the number of modules
#'  calculated. \code{n_genes} tells the number of genes included in the
#'  modules.
#'
get_wgcna <- function(data,
                      n_genes,
                      soft_power,
                      min_module_size) {
  max_gene_wgcna <- 3000
  # http://pklab.med.harvard.edu/scw2014/WGCNA.html
  if (n_genes > dim(data)[1]) {
    # Max	as data
    n_genes <- dim(data)[1]
  }
  if (n_genes < 50) {
    return(NULL)
  }
  if (dim(data)[2] < 4) {
    return(NULL)
  }
  if (n_genes > max_gene_wgcna) {
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

#' Dendogram plot of WGCNA modules
#'
#' Create a dendogram of the WGCNA return that color codes the modules and
#' includes a hierarchical dendogram.
#'
#' @param wgcna List returned from the \code{get_wgcna}
#'
#' @export
#' @return A dendogram plot of hierarchical clustering with a color bar to
#'  identify the modules.
get_module_plot <- function(wgcna) {
  diss <- 1 - wgcna$tom
  dynamic_colors <- wgcna$dynamic_colors

  hier <- flashClust::flashClust(as.dist(diss), method = "average")

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

#' Network of top genes
#'
#' Create a network plot of the top genes found with the WGCNA package.
#'
#' @param select_wgcna_module The module to create a plot of the top genes for,
#'  options can be found with the \code{get_wgcna_modules} function
#' @param wgcna List returned from the \code{get_wgcna}
#' @param top_genes_network Number of genes to include in the network plot
#' @param select_org Organism the expression data is for
#' @param all_gene_info Gene info that was found from querying the database
#' @param edge_threshold Value from 1 to 0.1 (0.4 recommended)
#'
#' @export
#' @return An adjacency matrix of the top genes in the selected module
get_network <- function(select_wgcna_module,
            wgcna,
            top_genes_network,
            select_org,
            all_gene_info,
            edge_threshold) {
  module <- unlist(strsplit(select_wgcna_module, " "))[2]
  module_colors <- wgcna$dynamic_colors
  in_module <- (module_colors == module)
  
  if (select_wgcna_module == "Entire network") {
    in_module <- rep(TRUE, length(in_module))
  }
  dat_expr <- t(wgcna$data)
  probes <- colnames(dat_expr)
  mod_probes <- probes[in_module]

  mod_tom <- wgcna$tom[in_module, in_module]
  dimnames(mod_tom) <- list(mod_probes, mod_probes)

  n_top <- min(top_genes_network, 1000)
  im_conn <- WGCNA::softConnectivity(dat_expr[, mod_probes])
  top <- (rank(-im_conn) <= n_top)

  # Adding symbols
  probe_to_gene <- NULL
  dim_all_gene_info <- dim(all_gene_info)
  if (select_org != "NEW" &&
    !is.null(dim_all_gene_info) &&
    dim_all_gene_info[1] > 1) {
  # If more than 50% genes have symbol
  if (sum(is.na(all_gene_info$symbol)) / dim_all_gene_info[1] < .5) {
    probe_to_gene <- all_gene_info[, c("ensembl_gene_id", "symbol")]
    probe_to_gene$symbol <- gsub(" ", "", probe_to_gene$symbol)

    ix <- which(
    is.na(probe_to_gene$symbol) |
      nchar(probe_to_gene$symbol) < 2 |
      toupper(probe_to_gene$symbol) == "NA" |
      toupper(probe_to_gene$symbol) == "0"
    )
    # Use gene ID
    probe_to_gene[ix, 2] <- probe_to_gene[ix, 1]
  }
  }

  net <- mod_tom[top, top]

  if (!is.null(probe_to_gene)) {
    ix <- match(colnames(net), probe_to_gene[, 1])
    colnames(net) <- dplyr::case_when(
      !is.na(probe_to_gene[ix, 2]) ~ probe_to_gene[ix, 2],
      TRUE ~ colnames(net)
    )
    ix <- match(rownames(net), probe_to_gene[, 1])
    rownames(net) <- dplyr::case_when(
      !is.na(probe_to_gene[ix, 2]) ~ probe_to_gene[ix, 2],
      TRUE ~ rownames(net)
    )
  }

  return(net)
}

#' Network of top genes, plot
#'
#' Create a network plot of the top genes found with the WGCNA package.
#'
#' @param adjacency_matrix igraph object
#' @param edge_threshold Value from 1 to 0.1 (0.4 recommended)
#' @export
#' @return A function that can be stored as an object and then called to produce
#'  the plot that the function created. If it is not stored and called, the
#'  function will only return another function.
get_network_plot <- function(adjacency_matrix, edge_threshold) {
  # plot using igraph, use the correlation as edge weight
  adjacency_matrix <- adjacency_matrix > edge_threshold
  for (i in 1:dim(adjacency_matrix)[1]) {
    # Remove self-connection
    adjacency_matrix[i, i] <- FALSE
  }
  graph <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
  # http://www.kateto.net/wp-content/uploads/2016/01/NetSciX_2016_Workshop.pdf
  net_plot <- function() {
    plot(
      graph,
      vertex.label.color = "black",
      vertex.label.dist = 3,
      vertex.size = 7
    )
  }
  return(net_plot)
}




#' List WGCNA modules
#'
#' Get the options for modules to select from running the \code{get_wgcna}
#' function.
#'
#' @param wgcna List returned from the \code{get_wgcna}
#'
#' @export
#' @return A character vector with all the strings that can be filled into the
#'  inut parameter \code{select_wgcna_module} in other WGCNA functions.
get_wgcna_modules <- function(wgcna) {
  if (dim(wgcna$module_info)[1] == 0) {
    # If no module
    return(NULL)
  } else {
    modules <- unique(wgcna$module_info[, c("dynamic_mods", "dynamic_colors")])
    module_list <- apply(modules, 1, paste, collapse = ". ")
    module_list <- paste0(
      module_list,
      " (",
      table(wgcna$module_info[, "dynamic_mods"]),
      " genes)"
    )
    module_list <- c(module_list, "Entire network")

    return(module_list)
  }
}

#' Gene vector query for enrichment
#'
#' Select a module to create a vector of gene IDs to use in an enrichment
#' analysis.
#'
#' @param select_wgcna_modules The module to create a plot of hte top genes for,
#'  options can be found with the \code{get_wgcna_modules} function
#' @param wgcna List returned from the \code{get_wgcna}
#'
#' @export
#' @return A vector of genes that are included in the selected module.
network_enrich_data <- function(select_wgcna_module,
                                wgcna) {
  module <- unlist(strsplit(select_wgcna_module, " "))[2]
  module_colors <- wgcna$dynamic_colors
  in_module <- (module_colors == module)

  if (select_wgcna_module == "Entire network") {
    in_module <- rep(TRUE, length(in_module))
  }

  probes <- rownames(wgcna$data)
  query <- probes[in_module]
  return(query)
}

#' Scale independence plot
#'
#' Using the WGCNA return, create a ggplot of scale independence.
#'
#' @param wgcna List returned from the \code{get_wgcna}
#'
#' @export
#' @return A formatted ggplot displaying the scale independence for the
#'  \code{get_wgcna} function return.
plot_scale_independence <- function(wgcna) {
  sft <- wgcna$sft
  powers <- wgcna$powers

  scale_plot <- ggplot2::ggplot(
    data = sft$fitIndices,
    ggplot2::aes(
      x = sft$fitIndices[, 1],
      y = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
    )
  ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = sft$fitIndices[, 1],
        y = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        label = powers,
        color = "red"
      )
    ) +
    ggplot2::labs(
      x = "Soft Threshold (power)",
      y = "Scale Free Topology Model Fit, signed R^2",
      title = "Scale independence"
    ) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = .80, color = "red")) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "none",
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      )
    )

  return(scale_plot)
}

#' Mean connectivity plot
#'
#' Create a ggplot from the wgcna object to display the mean connectivity.
#'
#' @param wgcna List returned from the \code{get_wgcna}
#'
#' @export
#' @return A formatted ggplot displaying the mean connectivityfor the
#'  \code{get_wgcna} function return.
plot_mean_connectivity <- function(wgcna) {
  sft <- wgcna$sft
  powers <- wgcna$powers


  connectivity_plot <- ggplot2::ggplot(
    data = sft$fitIndices,
    ggplot2::aes(
      x = sft$fitIndices[, 1],
      y = sft$fitIndices[, 5]
    )
  ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = powers,
        color = "red"
      )
    ) +
    ggplot2::labs(
      x = "Soft Threshold (power)",
      y = "Mean Connectivity",
      title = "Mean connectivity"
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "none",
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      )
    )

  return(connectivity_plot)
}

#' prepare_module_csv
#'
#' Create a dataframe from the wgcna object to export as a CSV file, this code is from old idep
#'
#' @param wgcna List returned from the \code{get_wgcna}
#' @param select_org Organism the expression data is for
#' @param all_gene_info Gene info that was found from querying the database
#' 
#' @export
#' @return A dataframe containing for csv file
prepare_module_csv <- function(wgcna,
select_org,
all_gene_info) {
  df <- merge(wgcna$module_info, wgcna$data, by.y = "row.names", by.x = "sub_gene_names", all.x = TRUE)
  dim_all_gene_info <- dim(all_gene_info)
  if (select_org != "NEW" &&
    !is.null(dim_all_gene_info) &&
    dim_all_gene_info[1] > 1) {
      df <- merge(df, all_gene_info, by.x = "sub_gene_names", by.y = "ensembl_gene_id", all.x = T)
      rownames(df) <- paste0(df$symbol, "__", df$sub_gene_names)
      # Convert row names to a column and name it
      df$symbol__gene_id <- rownames(df)
      # Reorder the columns to have "symbol__gene_id" as the first column
      df <- df[, c(ncol(df), 1:(ncol(df) - 1))]
  }
  df <- df[order(df$dynamic_mods), ]
  colnames(df)[2:4] <- c("gene_id", "module_color", "module")
  return(df)
}

#' prepare_module_csv_filter
#'
#' Create a dataframe from the wgcna object to export as a CSV file, allow the user to filter the module
#'
#' @param wgcna List returned from the \code{get_wgcna}
#' @param module The module to create a plot of hte top genes for,
#'
#' @export
#' @return A dataframe containing for csv file
prepare_module_csv_filter <- function(
    module_data,
    module) {
  if (module == "ALL") {
    return(module_data)
  }
  return(subset(module_data, module == module))
}
