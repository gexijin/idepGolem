#' fct_06_pathway.R This file holds all of the main data analysis functions
#' associated with eighth tab of the iDEP website.
#'
#'
#' @section fct_06_pathway.R functions:
#'
#'
#' @name fct_06_pathway.R
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
    top_1 <- top_1[1:n_pathway_show, ,drop=FALSE]
  }		 
  top_1 <- as.data.frame(top_1)
  top_1 <- cbind(rep(select_contrast, dim(top_1)[1]), row.names(top_1), top_1) 
  top_1$statistic <- as.character(round(as.numeric(top_1$statistic), 4)) 
  top_1$adj.Pval <- sprintf("%-2.1e", as.numeric(top_1$adj.Pval))
  top_1[, 2] <- as.character(top_1[, 2])
  top_1[, 1] <- as.character(top_1[, 1])
  colnames(top_1)[1] <- "Direction"
  colnames(top_1)[2] <- paste("GAGE analysis:", gsub("-", " vs ", select_contrast))
  top_1[which(top_1[, 3] > 0), 1] <- "Up"
  top_1[which(top_1[, 3] < 0), 1] <- "Down"
  top_1 <- top_1[order(top_1[, 1], -abs(as.numeric(top_1[, 3]))), ]
  top_1[duplicated(top_1[, 1]), 1] <- ""
  rownames(top_1) <- seq(1, nrow(top_1), 1)

  return(top_1)
}

#' PGSEA FUNCTION
pgsea_data <- function(
  processed_data,
  gene_sets,
  my_range,
  pathway_p_val_cutoff,
  n_pathway_show
){
  subtype = detect_groups(colnames(processed_data))
	
  # Cut off to report in PGSEA. Otherwise NA
	p_value <- 0.01  
	if(length(gene_sets) == 0) {
    return(list(pg3 = NULL, best = 1))
  }
	
	pg <- PGSEA::PGSEA(
    processed_data - rowMeans(processed_data),
    cl = gene_sets,
    range = my_range,
    p.value = TRUE,
    weighted = FALSE
  )
	
	pg_results <- pg$results
  # Remove se/wrts with all missing(non-signficant)
	pg_results <- pg_results[rowSums(is.na(pg_results)) < ncol(pg_results), ]
	if(dim(pg_results)[1] < 2) {
    return()
  }
	best <- max(abs(pg_results))
	
	if(length(subtype) < 4 || length(unique(subtype)) <2 ||
     length(unique(subtype)) == dim(processed_data)[2]) { 
	  pg_results <- pg_results[order(-apply(pg_results, 1, sd)), ]
	  return(list(pg_data = pg_results[1:top, ], best <- best ))
	} 
    
	cat("\nComputing P values using ANOVA\n");
	path_p_value <- function (
    k,
    pg_results,
    subtype
  ){
	  return(summary(aov(pg_results[k, ] ~ subtype))[[1]][["Pr(>F)"]][1])
	}
	p_values <- sapply(1:dim(pg_results)[1], function(x) {
    path_p_value(
      k = x,
      pg_results = pg_results,
      subtype = subtype)
  })
	p_values <- stats::p.adjust(p_values, "fdr")
	

  if(sort(p_values)[2] > pathway_p_val_cutoff) {
    return(list(pg_data = NULL, best = best)) 
  } else {  
    n_sig_t <- rowSums(pg$p.results < p_value)
	
	  result <- cbind(as.matrix(p_values), n_sig_t, pg_results) 
	  result <- result[order(result[, 1]), ]
    result <- result[which(result[, 1] < pathway_p_val_cutoff), , drop = F]
	
	  pg_results = result[, -2]

	  # When there is only 1 left in the matrix pg_results becomes a vector
	  if(sum(p_values < pathway_p_val_cutoff) == 1) {
      pg_data <- t(as.matrix(pg_results))
      pg_data <- rbind(pg_data, pg_data)
    } else {
      if(dim(results)[1] > n_pathway_show) {
        pg_data <- pg_results[1:n_pathway_show, ]
      } else {
        pg_data <- pg_results
      }
    }

	  rownames(pg_data) <- sapply(rownames(pg_data) , extract_under)
	  a <- sprintf("%-3.2e", pg_data[, 1])
	  rownames(pg_data) <- paste(a, rownames(pg_data), sep = " ")
	  pg_data <- pg_data[, -1]
	
    # Sort by SD
	  pg_data <- pg_data[order(-apply(pg_data, 1, sd)), ]
  
    return(list(
      pg_data = pg_data,
      best = best
    ))
  }
 }

#' PLOT PGSEA
plot_pgsea <- function(
  my_range,
  processed_data,
  contrast_samples,
  gene_sets,
  pathway_p_val_cutoff,
  n_pathway_show
) {
	genes <- processed_data[, contrast_samples]	
	if(length(gene_sets)  == 0)  {
    return(
      NULL
    )
  } else {
	  subtype = detect_groups(colnames(genes))
	  result <- pgsea_data(
      processed_data = genes,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show
    )
					 
	  if(is.null(result$pg_data)) {
      return(
        NULL
      )
    } else {
      PGSEA::smcPlot(
        result$pg_data,
        factor(subtype),
        scale = c(-max(result$pg_data), max(result$pg_data)),
        show.grid = T,
        margins = c(3,1, 13, 38),
        col = .rwb,
        cex.lab = 0.5
      )
    }
  }
}

#' FGSEA DATA
fgsea_data <- function(
  select_contrast,
  my_range,
  limma,
  gene_p_val_cutoff,
  gene_sets,
  absolute_fold,
  pathway_p_val_cutoff,
  n_pathway_show
) {
	no_sig <- as.data.frame("No significant pathway found.")
	if(length(limma$top_genes) == 0) {
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
	colnames(top_1) <- c("Fold","FDR")
	  
	# Remove some genes
	top_1 <- top_1[which(top_1$FDR < gene_p_val_cutoff), ]
	  
  if(length(gene_sets) == 0) {
    return(as.data.frame("No gene set found!"))
  }


	fold <- top_1[, 1]
  names(fold) <- rownames(top_1)
	
  # Use absolute value of fold change, disregard direction
  if(absolute_fold) {
    fold <- abs(fold)
  }
	 
  paths <- fgsea::fgsea(
    pathways = gene_sets, 
    stats = fold,
    minSize = my_range[1],
    maxSize = my_range[2],
    nPerm = 100000                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ,
		nproc = 6
  )
	
  if(dim(paths)[1] < 1) {
    return(no_sig)
  }
	paths <- as.data.frame(paths)
  # Sort by NES
  paths <- paths[order(-abs(paths[, 5])), ]
	top_1 <- paths[, c(1, 5, 7, 3)]
	colnames(top_1) <- c("Pathway", "NES", "Genes", "adj.Pval")
	  
	if(length(which(top_1[, 4] <= pathway_p_val_cutoff)) == 0) {
    return(no_sig)
  }
	top_1 <- top_1[which(top_1[, 4] <= pathway_p_val_cutoff), , drop = FALSE]
	
  if(dim(top_1)[1] > n_pathway_show) {
    top_1 <- top_1[1:n_pathway_show, , drop = FALSE]
  }
	
	top_1 <- as.data.frame(top_1)
	top_1 <- cbind(rep(select_contrast, dim(top_1)[1]), top_1) 
	top_1[, 4] <- as.character(round(as.numeric(top_1[, 4]), 4)) 
	top_1$adj.Pval <- sprintf("%-2.1e", as.numeric(top_1$adj.Pval))
	top_1[, 1] <- as.character(top_1[, 1])
	colnames(top_1)[1] <- "Direction"
	colnames(top_1)[2] <- paste("GSEA analysis:", gsub("-"," vs ", select_contrast))
	top_1[which(as.numeric(top_1[, 3]) > 0), 1] <- "Up"
  top_1[which(as.numeric(top_1[, 3]) < 0), 1] <- "Down"
	top_1 <- top_1[order(top_1[, 1], -abs(as.numeric(top_1[, 3]))), ]
	top_1[duplicated(top_1[, 1]), 1] <- ""	 
	top_1[, 3] <- as.character(round(as.numeric(top_1[, 3]), 4))

	return(top_1)
}

#' REACTOME DATA
reactome_data <- function(
  select_contrast,
  my_range,
  limma,
  gene_p_val_cutoff,
  converted,
  idep_data,
  pathway_p_val_cutoff,
  n_pathway_show,
  absolute_fold
) {
  browser()
  ensembl_species <- c(
    "hsapiens_gene_ensembl","rnorvegicus_gene_ensembl", "mmusculus_gene_ensembl",
	  "celegans_gene_ensembl","scerevisiae_gene_ensembl", "drerio_gene_ensembl",
    "dmelanogaster_gene_ensembl"
  )
  reactome_pa_species <- c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly" )
	no_sig <- as.data.frame("No significant pathway found.")
	if(length(limma$top_genes) == 0) {
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

	# Remove some genes
	top_1 <- top_1[which(top_1$FDR < gene_p_val_cutoff), ]
	 
	fold <- top_1[, 1]
  names(fold) <- rownames(top_1)
	if(absolute_fold) {
    # Use absolute value of fold change, disregard direction
    fold <- abs(fold) 
  }
  
  species <- converted$species[1, 1]
  ix <- match(species, ensembl_species)	
  
  if(is.na(ix)) {
    return(as.data.frame("Species not coverted by ReactomePA package!"))
  }
	  
	fold <- convert_ensembl_to_entrez(
    query = fold,
    species = species,
    org_info = idep_data$org_info
  )  
	
  fold <- sort(fold, decreasing = T)
	paths <- ReactomePA::gsePathway(
    fold,
    nPerm = 5000,
    organism = reactome_pa_species[ix],
    minGSSize = my_range[1], 
		maxGSSize = my_range[2],
		pvalueCutoff = 0.5,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  
  paths <- as.data.frame(paths)
	  
	if(is.null(paths)) {
    return(no_sig)
  }
	if(dim(paths)[1] == 0) {
    return(no_sig)
  }

  if(dim(paths)[1] < 1) {
    return(no_sig)
  }
  paths <- as.data.frame(paths)
  paths <- paths[order(-abs(paths[, 5])), ]

	top_1 <- paths[, c(2, 5, 3, 7)]

	colnames(top_1) <- c("Pathway", "NES", "Genes", "adj.Pval")
	  
	if(length(which(top_1[, 4] <= pathway_p_val_cutoff)) == 0) {
    return(no_sig)
  }  
  top_1 <- top_1[which(top_1[, 4] <= pathway_p_val_cutoff), , drop = FALSE]
	
  if(dim(top_1)[1] > n_pathway_show) {
    top_1 <- top_1[1:n_pathway_show, , drop = FALSE]
  }
 	
	top_1 <- as.data.frame(top_1)
	top_1 <- cbind(rep(select_contrast, dim(top_1)[1]), top_1) 
	top_1[, 4] <- as.character(round(as.numeric(top_1[, 4]), 4)) 
	top_1$adj.Pval <- sprintf("%-2.1e", as.numeric(top_1$adj.Pval))
	top_1[, 1] <- as.character(top_1[, 1])
	colnames(top_1)[1] <- "Direction"
	colnames(top_1)[2] <- paste(
    "ReactomePA analysis:",
    gsub("-"," vs ", select_contrast)
  )
	top_1[which(as.numeric(top_1[, 3]) > 0), 1] <- "Up"
	top_1[which(as.numeric(top_1[, 3]) < 0), 1] <- "Down"
	top_1 <- top_1[order(top_1[, 1], -abs(as.numeric(top_1[, 3]))), ]
	top_1[duplicated(top_1[, 1]), 1] <- ""	 
  top_1[, 3] <- as.character(round(as.numeric(top_1[, 3]), 4))
	
  return(top_1)
}