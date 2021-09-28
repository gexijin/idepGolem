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


##### work in process, want to rewrite other code that this function relays on first
### Can we removed this?????
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param x DESCRIPTION.
#' @param bar DESCRIPTION.
#' @param n DESCRIPTION.
#' @param mycolor DESCRIPTION.
#' @param clusterNames DESCRIPTION.
#' @param sideColors DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
gene_group_heatmap <- function(x, bar = NULL, n = -1, mycolor = 1, clusterNames = NULL, sideColors = NULL) {
  # number of genes to show
  ngenes <- as.character(table(bar))
  if (length(bar) > n && n != -1) {
    ix <- sort(sample(1:length(bar), n))
    bar <- bar[ix]
    x <- x[ix, ]
  }
  if (!is.null(bar)) {
    if (is.null(sideColors)) {
      sideColors <- mycolors
    }
  }

  # this will cutoff very large values, which could skew the color
  x <- as.matrix(x) - apply(x, 1, mean)
  cutoff <- median(unlist(x)) + 3 * sd(unlist(x))
  x[x > cutoff] <- cutoff
  cutoff <- median(unlist(x)) - 3 * sd(unlist(x))
  x[x < cutoff] <- cutoff
  # colnames(x)= detectGroups(colnames(x))
  if (is.null(bar)) { # no side colors
    heatmap.2(x,
      Rowv = F, Colv = F, dendrogram = "none",
      col = heatColors[as.integer(mycolor), ], density.info = "none", trace = "none", scale = "none", keysize = .3,
      key = F, labRow = F
      # ,RowSideColors = mycolors[bar]
      , margins = c(8, 24),
      srtCol = 45
    )
  } else {
    heatmap.2(x,
      Rowv = F, Colv = F, dendrogram = "none",
      col = heatColors[as.integer(mycolor), ], density.info = "none", trace = "none", scale = "none", keysize = .3,
      key = F, labRow = F,
      RowSideColors = sideColors[bar],
      margins = c(8, 24),
      srtCol = 45
    )
  }

  if (!is.null(bar)) {
    legend.text <- paste("Cluster ", toupper(letters)[unique(bar)], " (N=", ngenes, ")", sep = "")
    if (!is.null(clusterNames) && length(clusterNames) >= length(unique(bar))) {
      legend.text <- paste(clusterNames[1:length(unique(bar))], " (N=", ngenes, ")", sep = "")
    }

    par(lend = 1) # square line ends for the color legend
    legend("topright", # location of the legend on the heatmap plot
      legend = legend.text, # category labels
      col = sideColors, # color key
      lty = 1, # line style
      lwd = 10
    ) # line width
  }
}

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

#' FIND OVERLAP FUNCTION
#' Main function. Find a query set of genes enriched with functional category
find_overlap <- function(
  converted,
  all_gene_names,
  g_info,
  GO,
  selectOrg,
  minFDR,
  reduced = FALSE
) {
	
        
	if(length(query_set) == 0) {
    return(id_not_recognized)
  }
	
	ix = grep(converted$species[1,1], gmt_files)
	total_genes <- converted$species[1, 7]
	
	if (length(ix) == 0) {
    return(id_not_recognized)
  }
	
	# If selected species is not the default "bestMatch", use that species directly
	if(selectOrg != speciesChoice[[1]]) {  
		ix = grep(findSpeciesById(selectOrg)[1,1], gmtFiles )
		if (length(ix) == 0 ) {return(idNotRecognized )}
		totalGenes <- orgInfo[which(orgInfo$id == as.numeric(selectOrg)),7]
	}
	pathway <- dbConnect(sqlite,gmtFiles[ix],flags=SQLITE_RO)
	
		
	sqlQuery = paste( " select distinct gene,pathwayID from pathway where gene IN ('", paste(querySet,collapse="', '"),"')" ,sep="")
	
	#cat(paste0("HH",GO,"HH") )
	
	if( GO != "All") sqlQuery = paste0(sqlQuery, " AND category ='",GO,"'")
	result <- dbGetQuery( pathway, sqlQuery  )
	if( dim(result)[1] ==0) {return(as.data.frame("No matching species or gene ID file!" )) }

	# given a pathway id, it finds the overlapped genes, symbol preferred
	
	
	x0 = table(result$pathwayID)					
	x0 = as.data.frame( x0[which(x0>=Min_overlap)] )# remove low overlaps
	if(dim(x0)[1] <= 5 ) return(idNotRecognized) # no data
	colnames(x0)=c("pathwayID","overlap")
	pathwayInfo <- dbGetQuery( pathway, paste( " select distinct id,n,Description from pathwayInfo where id IN ('", 
							paste(x0$pathwayID,collapse="', '"),   "') ",sep="") )
	
	x = merge(x0,pathwayInfo, by.x='pathwayID', by.y='id')
	#browser() 
	test <- totalGenes - length(querySet)
	if (test < 0) {
	  test <- 0
	}
	x$Pval=phyper(x$overlap-1,length(querySet),test,as.numeric(x$n), lower.tail=FALSE )
	x$FDR = p.adjust(x$Pval,method="fdr")
	x <- x[ order( x$FDR)  ,]  # sort according to FDR
	
	if(dim(x)[1] > maxTerms ) x = x[1:maxTerms,]	
	
	if(min(x$FDR) > minFDR) x=as.data.frame("No significant enrichment found!") else {
		x <- x[which(x$FDR < minFDR),] 

		x= cbind(x,sapply( x$pathwayID, sharedGenesPrefered ) )
		colnames(x)[7]= "Genes"
		x <- subset(x,select = c(FDR,overlap,n,description,Genes) )
		colnames(x) = c("Corrected P value (FDR)", "Genes in list", "Total genes in category","Functional Category","Genes"  )
		
		# remove redudant gene sets
		if(reduced != FALSE ){  # reduced=FALSE no filtering,  reduced = 0.9 filter sets overlap with 90%
			n=  nrow(x)
			tem=rep(TRUE,n )
			geneLists = lapply(x$Genes, function(y) unlist( strsplit(as.character(y)," " )   ) )
			for( i in 2:n)
				for( j in 1:(i-1) ) { 
				  if(tem[j]) { # skip if this one is already removed
					  commonGenes = length(intersect(geneLists[i] ,geneLists[j] ) )
					  if( commonGenes/ length(geneLists[j] ) > reduced )
						tem[i] = FALSE	
				  }			
				}								
			x <- x[which(tem),]		
		}
		

	}
			
	dbDisconnect(pathway)
	return(x )
} 

#' SHARED GENES FUNCTION
shared_genes_prefered <- function(
  pathway_id,
  result) {
		filter_ids <- result[which(result[, 2] == pathway_id), 1]
		ix = match(filter_ids, converted$conversionTable$ensembl_gene_id) # convert back to original
		tem2 <- unique(converted$conversionTable$User_input[ix] )
		if(length(unique(gInfo$symbol) )/dim(gInfo)[1] >.7  ) # if 70% genes has symbol in geneInfo
		{ ix = match(tem, gInfo$ensembl_gene_id); 
		tem2 <- unique( gInfo$symbol[ix] )      }
	return( paste( tem2 ,collapse=" ",sep="") )}