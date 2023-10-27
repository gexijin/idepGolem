#' fct_08_bicluster.R This file holds all of the main data analysis functions
#' associated with eighth tab of the iDEP website.
#'
#'
#' @section fct_07_bicluster.R functions:
#'
#'
#' @name fct_08_bicluster.R
NULL

#' Get bicuster information
#'
#' Utilize the biclust package to run biclustering on the processed data for a
#' selected number of genes.
#'
#' @param data Matrix of gene data that has been through the
#'  \code{\link{pre_process}()}
#' @param n_genes Number of most variable genes to include in the analysis
#' @param biclust_method String indicating the method of biclustering to perform.
#'  Should be one of "biclust::BCCC()", "QUBIC::BCQU()", "runibic::BCUnibic()",
#'  "biclust::BCXmotifs()", "biclust::BCPlaid()", "biclust::BCSpectral()",
#'  "biclust::BCBimax()", "biclust::BCQuest()"
#'
#' @export
#' @return A list containing two object. \code{data} is the submatrix from the
#'  processed data with the biclust function \code{discetize} performed on it.
#'  For information on this function visit
#'  https://rdrr.io/cran/biclust/man/discretize.html. \code{res} is the return
#'  from using the \code{biclust} function with the selected method. For
#'  information on this function visit
#'  https://www.rdocumentation.org/packages/biclust/versions/2.0.3/topics/biclust.
#'
#'
#' @family bicluster functions
#'
#' @seealso
#' * \code{\link[biclust]{BCCC}()}, \code{\link[biclust]{BCXmotifs}()},
#'  \code{\link[biclust]{BCPlaid}()}, \code{\link[biclust]{BCXmotifs}()},
#'  \code{\link[biclust]{BCPlaid}()}, \code{\link[biclust]{BCSpectral}()},
#'  \code{\link[biclust]{BCBimax}()}, \code{\link[biclust]{BCquest}()} for
#'  specific function documentation from the biclust package
#' * \code{\link[QUBIC]{BCQU}()}, \code{\link[QUBIC]{BCQU}()} for specific
#'  function documentation from the QUBIC package
#' * \code{\link[runibic]{BCUnibic}()}, \code{\link[runibic]{BCUnibic}()}
#'  for specific function documentation from the runibic package
#'
#' @md
get_biclustering <- function(data,
                             n_genes,
                             biclust_method) {
  if (n_genes > dim(data)[1]) {
    n_genes <- dim(data)[1]
  }
  if (n_genes < 10) {
    n_genes <- 10
  }
  if (n_genes > 2000) {
    n_genes <- 2000
  }
  data <- as.matrix(data[1:n_genes, ]) - apply(data[1:n_genes, ], 1, mean)

  if (biclust_method == "BCXmotifs()") {
    data <- biclust::discretize(data)
  }

  run_biclust <- paste(
    "res <- biclust::biclust(as.matrix(data), method = ", biclust_method, ")"
  )
  eval(parse(text = run_biclust))

  return(list(
    data = data,
    res = res
  ))
}

#' Summary of biclust
#'
#' Give the summary of the identified clusters from the \code{get_biclustering}
#' function which gives the number of clusters found and the samples used.
#'
#' @param biclustering Return list from
#'  \code{\link{get_biclustering}()}
#' @param select_bicluster Number of selected cluster to return information on
#'
#' @export
#' @return A character string message with the requested summary information.
#'
#' @family bicluster functions
#'
bicluster_summary_message <- function(biclustering,
                                      select_bicluster) {
  res <- biclustering$res
  if (res@Number == 0) {
    return(
      "No cluster found! Perhaps sample size is too small."
    )
  } else {
    data <- biclust::bicluster(
      biclustering$data,
      res,
      as.numeric(select_bicluster)
    )[[1]]
    return(
      paste(
        res@Number,
        "clusters found. Cluster ",
        select_bicluster,
        " has",
        dim(data)[1],
        "genes correlated across",
        dim(data)[2],
        "samples."
      )
    )
  }
}

#' Table of biclustering data
#'
#' This function processes results from the biclustering functions to be
#' displayed in a table.
#'
#' @param res The \code{res} object returned in the list output from
#'  \code{\link{get_biclustering}()}
#' @param biclust_data Matrix of data returned from
#'  \code{\link[biclust]{bicluster}()}
#' @param select_org String indicating selected organism for the expression
#'  data.
#' @param all_gene_info Gene information matched species in iDEP database from
#'  \code{\link{gene_info}()}
#'
#' @export
#' @returns A dataframe of results
#'
#' @family bicluster functions
#'
get_biclust_table_data <- function(res,
                                   biclust_data,
                                   select_org,
                                   all_gene_info) {
  if (res@Number == 0) {
    return(as.data.frame("No clusters found!"))
  }

  dim_all_gene_info <- dim(all_gene_info)
  if (select_org == "NEW" |
    (!is.null(dim_all_gene_info) && dim_all_gene_info[1] == 1)) {
    biclust_genes <- as.data.frame(rownames(biclust_data))
    colnames(biclust_genes) <- "Genes"
  } else {
    clust_info <- merge(
      biclust_data,
      all_gene_info,
      by.x = "row.names",
      by.y = "ensembl_gene_id",
      all.x = T
    )

    clust_info <- clust_info[!duplicated(clust_info$Row.names), ]

    if (sum(is.na(clust_info$band)) == dim(clust_info)[1]) {
      clust_info$chr <- clust_info$chromosome_name
    } else {
      clust_info$chr <- paste(
        clust_info$chromosome_name,
        clust_info$band,
        sep = ""
      )
    }

    gene_mean <- rowMeans(biclust_data)

    clust_info <- cbind(clust_info, gene_mean)

    clust_info <- clust_info[order(
      clust_info$gene_mean,
      decreasing = T
    ), ]

    clust_info <- clust_info[, c(
      "Row.names",
      "gene_mean",
      "symbol",
      "chr",
      "gene_biotype"
    )]

    clust_info$gene_mean <- sprintf(
      "%-4.2f",
      as.numeric(clust_info$gene_mean)
    )

    colnames(clust_info) <- c(
      "Ensembl ID",
      "Mean",
      "Symbol",
      "Chr",
      "Type"
    )

    if (sum(is.na(clust_info$Symbol)) == dim(clust_info)[1]) {
      clust_info <- clust_info[, -3]
    }
  }

  return(clust_info)
}
