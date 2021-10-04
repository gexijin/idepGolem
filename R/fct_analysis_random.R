#' fct_analysis_random.R Miscellaneous data analysis functions here,
#'  we find a place for them later
#'
#'
#' @section fct_analysis_random.R functions:
#' \code{gene_group_heatmap} heatmap with color bar define gene groups
#'
#' \code{find_overlap_gmt} Given a gene set,
#' finds significant overlaps with a gene set database.
#'
#'
#' @name fct_analysis_random.R
NULL

### work in process
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param query DESCRIPTION.
#' @param gene_set DESCRIPTION.
#' @param min_fdr DESCRIPTION.
#' @param min_size DESCRIPTION.
#' @param max_size DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
find_overlap_gmt <- function(
  query, 
  gene_set,
  min_fdr = .2,
  min_size = 2,
  max_size = 10000
) {
  total_elements <- 30000 # why 3000?
  min_overlap <- 1 # nolint
  max_terms <- 10 # max number of enriched terms should be user input ????
  no_sig <- as.data.frame("No significant enrichment found!")
  query <- clean_gene_set(gene_set = query) # convert to upper case, unique()
  query_length <- length(query)

  if (query_length <= 2 || length(gene_set) < 1) {
    return(no_sig)
  }

  # gene sets smaller than 1 is ignored!!!
  gene_set <- gene_set[which(sapply(X = gene_set, FUN = length) > min_size)]
  # gene sets smaller than 1 is ignored!!!
  gene_set <- gene_set[which(sapply(X = gene_set, FUN = length) < max_size)]

  foo <- function(x, query) {
    length(intersect(query, x))
  }
  result <- unlist(lapply(X = gene_set, FUN = foo, query = query))
  result <- cbind(unlist(lapply(X = gene_set, FUN = length)), result)
  result <- result[which(result[, 2] > min_overlap), , drop = FALSE]
  if (nrow(result) == 0) {
    return(no_sig)
  }
  xx <- result[, 2]
  nn <- total_elements - query_length
  kk <- result[, 1]
  pval_enrich <- phyper(xx - 1, query_length, nn, kk, lower.tail = FALSE)
  fdr <- p.adjust(pval_enrich, method = "fdr", n = length(gene_set))
  result <- as.data.frame(cbind(fdr, result))
  result <- result[, c(1, 3, 2)]
  result$pathway <- rownames(result)
  result$Genes <- "" # place holder just
  colnames(result) <- c(
    "Corrected P value (FDR)",
    "Genes in list", "Total genes in category", "Functional Category", "Genes"
  )
  result <- result[which(result[, 1] < min_fdr), , drop = FALSE]
  if (nrow(result) == 0 || min(fdr) > min_fdr) {
    return(no_sig)
  }
  result <- result[order(result[, 1]), ]
  if (nrows(result) > max_terms) {
    result <- result[1:max_terms, ]
  }
  return(result)
}

#' Find overlap for pathway analysis
#' 
#' Use the pathway table from the read_pathway_sets function
#' to calculate adjusted p-values. Adjusted p-values determine
#' the enriched pathways from the selected qeury.
#' 
#' @param pathway_table Results from the read_pathway_sets
#'   query. If NULL or 0 rows there is no significant
#'   enrichment
#' @param query_set The vector of IDs that the enrichment
#'   analysis is being performed on
#' @param total_genes Length of the query set subtracted from
#'   the total number of genes in the database. Could change
#'   within the function if the background set changes to the
#'   filtered genes.
#' @param processed_data Data that has been filtered and
#'   transformed in the pre_process function
#' @param gene_info The gene info from the converted IDs and
#'   the function gene_info()
#' @param go Section of the database to query for pathway
#'   analysis
#' @param idep_data Data built in to idep
#' @param sub_pathway_files GMT files in the database that
#'   contain information for the matched species
#' @param use_filtered_background T/F Use the genes that
#'   passed the pre_process filter as the backgrounf
#' @param select_org Input for what organism the IDs are 
#'   pertaining to
#' @param reduced T/F Remove gene sets that are redudant
#'   from the final result
#' 
#' @return If there is significant enrichment, the data frame
#'   that is returned has a pathway for each row with the
#'   total genes in the database mapping to it as well as the
#'   number of genes in the query that map to it. It also
#'   contains a column for the p-value and a list of the
#'   specific IDs included in the pathway from the query.
find_overlap <- function(
  pathway_table,
  query_set,
  total_genes,
  processed_data,
  gene_info,
  go,
  idep_data,
  sub_pathway_files,
  use_filtered_background,
  select_org,
  reduced = FALSE
) {
  max_pval_filter <- 0.3
  max_genes_background <- 30000
  min_genes_background <- 2000
  max_terms <- 15
  min_fdr <- .05
  if(reduced) {
    reduced <- .9
  }
  error_msg <- data.frame("Enrichment" = "No significant enrichment found!")

  if(select_org == "NEW" && is.null(pathway_table)) {
    return(data.frame("Enrichment" = "No GMT file provided!"))
  } else if (select_org == "NEW") {
    find_overlap_gmt(
      query = query_set,
      gene_set = pathway_table,
      min_fdr = .2,
      min_size = 2,
      max_size = 10000
    )
  }

  if(dim(pathway_table)[1] == 0 || is.null(pathway_table)) {
    return(error_msg)
  }

  pathway_table <- pathway_table[which(
    pathway_table$overlap / length(query_set) /
    (as.numeric(pathway_table$n) / total_genes) > 1
  ), ]

  pathway_table$pval <- stats::phyper(
    pathway_table$overlap - 1,
		length(query_set),
		total_genes - length(query_set),   
		as.numeric(pathway_table$n), 
		lower.tail=FALSE
  )
  pathway_table <- subset(pathway_table, pathway_table$pval < max_pval_filter)

  # Background genes -----------
  if(!is.null(use_filtered_background)) {
    if(
      use_filtered_background && 
      length(row.names(processed_data)) > min_genes_background &&
      length(row.names(processed_data)) < max_genes_background + 1
    ) {
      pathway_table_bg <- background_pathway_sets(
        processed_data = processed_data,
        gene_info = gene_info,
        sub_query = query_set,
        go = go,
        pathway_table = pathway_table,
        idep_data = idep_data,
        sub_pathway_files = sub_pathway_files
      )

      pathway_table$pval <- phyper(
        pathway_table_bg$overlap - 1,
        length(query_set),
        length(rownames(processed_data)) - length(query_set),   
        as.numeric(pathway_table_bg$overlap_bg),
        lower.tail=FALSE
      ) 
    }
  }

  pathway_table$fdr <- stats::p.adjust(pathway_table$pval, method = "fdr")
  pathway_table <- pathway_table[order(pathway_table$fdr), ]

  if(dim(pathway_table)[1] > max_terms) {
    pathway_table <- pathway_table[1:max_terms, ]
  }
  
  if(min(pathway_table$fdr) > min_fdr) {
    pathway_table <- error_msg
  } else {
    pathway_table <- pathway_table[which(pathway_table$fdr < min_fdr), ]

    pathway_table <- subset(pathway_table, select = c(fdr, overlap, n, description, gene_sets) )
    pathway_table$n <- as.numeric(pathway_table$n)
    pathway_table$fdr <-  formatC(pathway_table$fdr, format = "e", digits = 2)
    colnames(pathway_table) <- c(
      "Corrected P value (FDR)", "Genes in query", "Total genes in category",
      "Functional Category", "Genes"
    )

    # Remove redudant gene sets
		if(reduced != FALSE ){
			n <- nrow(pathway_table)
			tem <- rep(TRUE, n)
			gene_lists = pathway_table$Genes
			for(i in 2:n) {
        for(j in 1:(i-1)) { 
				  if(tem[j]) {
					  common_genes = length(intersect(gene_lists[[i]], gene_lists[[j]]))
					  if(common_genes/ length(gene_lists[[j]]) > reduced) {
              tem[j] = FALSE
            }	
				  }			
				}				
      }
								
			pathway_table <- pathway_table[which(tem), ]		
		}
  }

  return(pathway_table)
}