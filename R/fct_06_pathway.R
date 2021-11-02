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
	
	  pg_results <- result[, -2]

	  # When there is only 1 left in the matrix pg_results becomes a vector
	  if(sum(p_values < pathway_p_val_cutoff) == 1) {
      pg_data <- t(as.matrix(pg_results))
      pg_data <- rbind(pg_data, pg_data)
    } else {
      if(dim(pg_results)[1] > n_pathway_show) {
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
	  subtype <- detect_groups(colnames(genes))
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
        margins = c(3, 1, 13, 38),
        col = PGSEA::.rwb,
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

pgsea_plot_all <- function(
  go,
  my_range,
  data,
  select_contrast,
  gene_sets,
  pathway_p_val_cutoff,
  n_pathway_show
) {
  if(length(gene_sets)  == 0) {
    plot.new()
    text(0, 1, "No gene sets!")
  } else {
    subtype <- detect_groups(colnames(data)) 
	  result <- pgsea_data(
      processed_data = data,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show
    )
    if(is.null(result$pg_data)) {
      plot.new()
      text(0.5, 1, "No significant pathway found!")
    } else {
      PGSEA::smcPlot(
        result$pg_data,
        factor(subtype),
        scale = c(-max(result$pg_data), max(result$pg_data)),
        show.grid = T,
        margins = c(3, 1, 13, 38),
        col = PGSEA::.rwb,
        cex.lab = 0.5
      )
    }
  }
}

get_pgsea_plot_data <- function(
  my_range,
  data,
  select_contrast,
  gene_sets,
  sample_info,
  select_factors_model,
  select_model_comprions,
  pathway_p_val_cutoff,
  n_pathway_show
) {	
	# Find sample related to the comparison
	iz <- match(detect_groups(colnames(data)), unlist(strsplit(select_contrast, "-")))
	iz <- which(!is.na(iz))
	
  if(!is.null(sample_info) & !is.null(select_factors_model) & length(select_model_comprions) > 0 ) {
		# Strings like: "groups: mutant vs. control"
    comparisons <- gsub(".*: ", "", select_model_comprions)
		comparisons <- gsub(" vs\\. ", "-", comparisons)
    # Corresponding factors
		factors_vector <- gsub(":.*", "", select_model_comprions)
    # Selected contrast lookes like: "mutant-control"
		ik <- match(select_contrast, comparisons)
		if(is.na(ik)) {
      iz <- 1:(dim(data)[2]) 
    } else {
      # Interaction term, use all samples	
      # Corresponding factors
			selected_factor <- factors_vector[ik]
			iz <- match(sample_info[, selected_factor], unlist(strsplit(select_contrast, "-")))
			iz <- which(!is.na(iz))
		}
	}

	if(grepl("I:", select_contrast)) {
    # If it is factor design use all samples
    iz <- 1:(dim(data)[2])
  } 
	if(is.na(iz)[1] | length(iz) <= 1) {
    iz <- 1:(dim(data)[2]) 
  }
	
	genes <- data
	genes <- genes[, iz]

	subtype <- detect_groups(colnames(genes)) 
  
  if(length(gene_sets)  == 0) {
    return(as.data.frame("No significant pathway!"))
  } else {
    result <- pgsea_data(
      processed_data = genes,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show
    )
					 
	  if(is.null(result$pg_data)) {
      return(as.data.frame("No significant pathway!"))
    } else {
      return( as.data.frame(result$pg_data) )
    } 
  }
}

get_pgsea_plot_all_samples_data <- function(
  data,
  select_contrast,
  gene_sets,
  my_range,
  pathway_p_val_cutoff,
  n_pathway_show
) {
  genes <- data
	subtype <- detect_groups(colnames(genes)) 
  
  if(length(gene_sets)  == 0) {
    plot.new()
    text(0,1, "No gene sets!")
  } else {
	  result <- pgsea_data(
      processed_data = genes,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show
    )
					 
	  if(is.null(result$pg_data)) {
      return(as.data.frame("No significant pathway!"))
    } else {
      return(as.data.frame(result$pg_data))
    }
  }	
}

pathway_select_data <- function(
  sig_pathways,
  gene_sets,
  contrast_samples,
  data,
  select_org,
  all_gene_names
) {
  if(sig_pathways == "All") {
    return (NULL) 
  }
  # Find the gene set
	ix <- which(names(gene_sets) == sig_pathways)
	if(length(ix) == 0) {
    return(NULL)
  }
  # Retrieve genes
  genes <- gene_sets[[ix]]

	# Find related samples	
	iz <-contrast_samples
	x <- data[which(rownames(data) %in% genes), iz]
	if(ncol(all_gene_names) == 3) {
    x <- rowname_id_swap(
      data_matrix = x,
      all_gene_names = all_gene_names,
      select_gene_id = "symbol"
    )
  }
	
	return(x)
}

# SELECTED HEATMAP
pathway_heatmap <- function(
  data,
  heatmap_color_select
) {
  # Number of genes to show
	n_genes <- nrow(data)

  data <- as.matrix(data) - apply(data, 1, mean)
  cutoff <- median(unlist(data)) + 3 * sd(unlist(data)) 
	data[data > cutoff] <- cutoff
	cutoff <- median(unlist(data)) - 3 * sd(unlist(data)) 
	data[data < cutoff] <- cutoff
	
	data <- data[which(apply(data, 1, sd) > 0), ]

  # Color scale
  if(min(data) < 0) {
    col_fun <- circlize::colorRamp2(
      c(min(data), 0, max(data)),
      heatmap_color_select
    )
  } else {
    col_fun <- circlize::colorRamp2(
      c(min(data), median(data), max(data)),
      heatmap_color_select
    )
  }
  
	groups <- detect_groups(colnames(data))
	group_count <- length(unique(groups))
  groups_colors <- gg_color_hue(group_count)
  
  top_ann <- ComplexHeatmap::HeatmapAnnotation(
    Group = groups,
    col = list(
      Group = setNames(groups_colors, unique(groups))
    ),
    annotation_legend_param = list(
      Group = list(nrow = 1, title = NULL)
    ),
    show_annotation_name = list(Group = FALSE),
    show_legend = FALSE
  )

  heat <- ComplexHeatmap::Heatmap(
    data,
    name = "Expression",
    col = col_fun,
    cluster_rows = TRUE,
    clustering_method_rows = "average",
    clustering_distance_rows = function(x) as.dist(
      1 - cor(t(x), method = "pearson")
    ),
    cluster_columns = TRUE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    top_annotation = top_ann,
    show_row_names = FALSE,
    show_column_names = FALSE,
    heatmap_legend_param = list(
      direction = "horizontal",
      legend_width = grid::unit(6, "cm"),
      title = "Color Key",
      title_position = "topcenter"
    )
  )
  
  return(
    heatmap = ComplexHeatmap::draw(
      heat,
      heatmap_legend_side = "bottom"
    )
  )
}


#' SUBHEATMAP
path_heat_sub <- function(
  ht_brush,
  ht,
  ht_pos_main,
  heatmap_data
) {
  lt <- InteractiveComplexHeatmap::getPositionFromBrush(ht_brush)
  pos1 <- lt[[1]]
  pos2 <- lt[[2]]

  pos <- InteractiveComplexHeatmap::selectArea(
    ht,
    mark = FALSE,
    pos1 = pos1,
    pos2 = pos2,
    verbose = FALSE,
    ht_pos = ht_pos_main
  )

  # Annotations ----------
    column_groups <- detect_groups(colnames(heatmap_data))
    groups_colors <- gg_color_hue(length(unique(column_groups)))
  
    top_ann <- ComplexHeatmap::HeatmapAnnotation(
      Group = column_groups,
      col = list(
        Group = setNames(
          groups_colors,
          unique(column_groups)
        )
      ),
      annotation_legend_param = list(
        Group = list(nrow = 1, title = NULL)
      ),
      show_annotation_name = list(Group = FALSE),
      show_legend = TRUE
    )
    
    group_col_return <- setNames(
      groups_colors,
      c(unique(column_groups))
    )
  # End annotation ---------

  column_index <- unlist(pos[1, "column_index"])
  row_index <- unlist(pos[1, "row_index"])
  top_ann <- top_ann[column_index]
  column_groups <- column_groups[column_index]
  m <- ht@ht_list[[1]]@matrix

  if (length(row_index) > 50) {
    show_rows <- FALSE
  } else {
    show_rows <- TRUE
  }

  submap_data <- m[row_index, column_index, drop = FALSE]

  ht_select <- ComplexHeatmap::Heatmap(
    submap_data,
    col = ht@ht_list[[1]]@matrix_color_mapping@col_fun,
    show_heatmap_legend = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = show_rows,
    top_annotation = top_ann,
    name = "heat_1"
  )

  return(list(
    ht_select = ht_select,
    submap_data = submap_data,
    group_colors = group_col_return,
    column_groups = column_groups
  ))
}

#' PATH SUB CLICK INFO
path_click_info <- function(
  click,
  ht_sub,
  ht_sub_obj,
  ht_pos_sub,
  sub_groups,
  group_colors,
  data
) {
  pos1 <- InteractiveComplexHeatmap::getPositionFromClick(click)
    
  pos <- InteractiveComplexHeatmap::selectPosition(
    ht_sub,
    mark = FALSE,
    pos = pos1,
    verbose = FALSE,
    ht_pos = ht_pos_sub
  )
  
  row_index <- pos[1, "row_index"]
  column_index <- pos[1, "column_index"]

  if (is.null(row_index)) {
    return("Select a cell in the heatmap.")
  }

  value <- data[row_index, column_index]
  col <- ComplexHeatmap::map_to_colors(ht_sub_obj@matrix_color_mapping, value)
  sample <- colnames(data)[column_index]
  gene <- rownames(data)[row_index]
  group_name <- sub_groups[column_index]
  group_col <- group_colors[[group_name]]

  # HTML for info table
  # Pulled from https://github.com/jokergoo/InteractiveComplexHeatmap/blob/master/R/shiny-server.R
  # Lines 1669:1678
  html <- GetoptLong::qq("
<div>
<pre>
Value: @{round(value, 2)} <span style='background-color:@{col};width=50px;'>    </span>
Sample: @{sample}
Gene: @{gene} 
Group: @{group_name} <span style='background-color:@{group_col};width=50px;'>    </span>   
</pre></div>"
)

 return(HTML(html))
}

get_pathway_list_data <- function(
  pathway_method,
  gage_pathway_data,
  fgsea_pathway_data,
  pgsea_plot_data,
  pgsea_plot_all_samples_data,
  go,
  select_org,
  gene_info
) {
	pathways <- NULL
	if(pathway_method == 1) {
    if(!is.null(gage_pathway_data)) {
      if(dim(gage_pathway_data)[2] > 1) { 
				pathways <- gage_pathway_data
				colnames(pathways)[2] <- "Pathways" 	
				colnames(pathways)[4] <- "nGenes"
			}
    }
  }
	if(pathway_method == 3) {
    if(!is.null(fgsea_pathway_data)) {
      if(dim(fgsea_pathway_data)[2] > 1) {
				pathways <- fgsea_pathway_data
				colnames(pathways)[2] <- "Pathways" 	
				colnames(pathways)[4] <- "nGenes" 
			}
    }
  }
	if(pathway_method == 2) {
    if(!is.null(pgsea_plot_data)) {
      if(dim(pgsea_plot_data)[2] > 1) {
				pathways <- as.data.frame(pgsea_plot_data)
				pathways$Pathways <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
				pathways$adj.Pval <- gsub(" .*", "", rownames(pathways))
				pathways$Direction <- "Diff"
			}
    }
  }
	if(pathway_method == 4) {
    if(!is.null(pgsea_plot_all_samples_data)) {
      if(dim(pgsea_plot_all_samples_data)[2] >1 ) {
				pathways <- as.data.frame(pgsea_plot_all_samples_data)
				pathways$Pathways <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
				pathways$adj.Pval <- gsub(" .*", "", rownames(pathways))
				pathways$Direction <- "Diff"	
			}
    }
  }	
	if(is.null(pathways)) {
    return(NULL)
  }	
	# if no gene set data, return pathway list
	if(is.null(gene_sets)) {
    return(pathways)
  } 
	
	pathways$adj.Pval <- as.numeric(pathways$adj.Pval)

  # Sometimes only one pathway is in the table
	if(nrow(pathways) > 1) {
    for(i in 2:nrow(pathways)) {
      if(nchar(pathways$Direction[i]) <= 1) {
			  pathways$Direction[i] = pathways$Direction[i-1]
      }
    }
	}	

	# Gene symbol matching symbols 
	probe_to_gene <- NULL
	if(go != "ID not recognized!" & select_org != "NEW") {
    # If more than 50% genes has symbol
    if(sum(is.na(gene_info$symbol)) / dim(gene_info)[1] < .5) { 
		  probe_to_gene <- gene_info[, c("ensembl_gene_id", "symbol")]
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
			

	
	pathways$Genes <- ""
	# looking up genes for each pathway
	for(i in 1:nrow(pathways)) {
    # Find the gene set
		ix <- which(names(gene_sets) == pathways$Pathways[i])
		if(length(ix) != 0) {
      # Retrieve genes
			genes <- gene_sets[[ix]]
			
			if(!is.null(probe_to_gene)) { 
				iy <- match(genes, probe_to_gene[, 1])
				genes <- probe_to_gene[iy, 2]
			}
			pathways$Genes[i] <- paste(genes, collapse = " ")
		}
	}
	return(pathways)
}