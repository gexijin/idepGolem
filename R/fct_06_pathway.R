#' fct_06_pathway.R This file holds all of the main data analysis functions
#' associated with eighth tab of the iDEP website.
#'
#'
#' @section fct_06_pathway.R functions:
#'
#'
#' @name fct_06_pathway.R
NULL

#' Pathway analysis with gage package
#'
#' Run pathway analysis with the gage package using the results
#' from the limma_value function.
#'
#' @param select_go String designating the section of the database to query for
#'   pathway analysis. See \code{\link{gmt_category}()} for choices.
#' @param select_contrast String designating the comparison from DEG analysis to
#'  filter for the significant genes. See the 'comparison' element from the list
#'  returned from \code{\link{limma_value}()} for options.
#' @param min_set_size Minimum gene set size for a pathway
#' @param max_set_size Maximum gene set size for a pathway
#' @param limma Results list from the \code{\link{limma_value}()}
#' @param gene_p_val_cutoff Significant p-value to filter
#'  the top genes fold change by
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database. See list returned from \code{\link{read_gene_sets}()}
#' @param absolute_fold TRUE/FALSE to use the absolute value of the fold
#'  change
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#'  result
#' @param expr_data Pre-processed gene expression data frame
#' @param sample_info Experiment design structure as a data frame
#' @param data_type Data type selected as integers; Fold Change = 1, 
#' Expression = 2
#'
#' @export
#' @return A data frame with the results of the pathway analysis.
#'  The data frame has five columns for the direction of the
#'  regulation, the pathway description, the stat value, the
#'  number of overlapping genes, and the p-value.
#'
#' @family pathway functions
gage_data <- function(expr_data,
                      sample_info = NULL,
                      data_type,
                      select_go,
                      select_contrast,
                      min_set_size,
                      max_set_size,
                      limma,
                      gene_p_val_cutoff,
                      gene_sets,
                      absolute_fold,
                      pathway_p_val_cutoff,
                      n_pathway_show) {
  if (select_go == "ID not recognized!") {
    return(as.data.frame("Gene ID not recognized."))
  }
  if (is.null(select_contrast)) {
    return(NULL)
  }
  my_range <- c(min_set_size, max_set_size)
  no_sig <- as.data.frame("No significant pathway found.")
  if (length(limma$top_genes) == 0) {
    return(no_sig)
  }
  
  if (length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]
  } else {
    top <- limma$top_genes
    ix <- match(select_contrast, names(top))
    if (is.na(ix)) {
      return(no_sig)
    }
    top_1 <- top[[ix]]
  }
  if (dim(top_1)[1] == 0) {
    return(no_sig)
  }
  colnames(top_1) <- c("Fold", "FDR")

  expr_data <- expr_data[rownames(top_1),]
  expr_data <- expr_data[which(top_1$FDR < gene_p_val_cutoff), ]
  top_1 <- top_1[which(top_1$FDR < gene_p_val_cutoff), ]

  gmt <- gene_sets
  if (length(gene_sets) == 0) {
    return(as.data.frame("No gene set found!"))
  }
  
  groups <- unlist(
    stringr::str_split(string = select_contrast, 
                       pattern = "-")
  )
  
  # Conditional statements for determining GAGE data
  # Expression data selected with experiment design/sample info
  if(!is.null(sample_info) && data_type == 2) {
    sample_info <- t(sample_info)
    
    samples <- sapply(groups, function(x){
      colnames(sample_info)[colSums(sample_info == x) > 0]
    })
    
    samples <- as.vector(samples)
    data <- expr_data[, samples]
    # Fold change selected
  } else if (data_type == 1){
    fold <- top_1[, 1]
    names(fold) <- rownames(top_1)
    if (absolute_fold) {
      fold <- abs(fold)
    } 
    data <- fold
    # Expression data selected without experiment design
  } else {
    col_groups <- detect_groups(colnames(expr_data), preserve_original = TRUE)
    
    cols <- sapply(groups, function(x){
      which(col_groups == x)
    })
    cols <- sort(as.vector(cols))
    
    data <- expr_data[, cols]
  }
  
  paths <- gage::gage(data, gsets = gmt, ref = NULL, samp = NULL)

  paths <- rbind(paths$greater, paths$less)

  if (dim(paths)[1] < 1 | dim(paths)[2] < 6) {
    return(no_sig)
  }
  top_1 <- paths[, c("stat.mean", "set.size", "q.val")]
  colnames(top_1) <- c("statistic", "Genes", "adj.Pval")
  top_1 <- top_1[order(top_1[, 3]), ]
  if (length(which(top_1[, 3] <= pathway_p_val_cutoff)) == 0) {
    return(no_sig)
  }
  top_1 <- top_1[which(top_1[, 3] <= pathway_p_val_cutoff), , drop = FALSE]
  if (dim(top_1)[1] > n_pathway_show) {
    top_1 <- top_1[1:n_pathway_show, , drop = FALSE]
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
  top_1 <- top_1[order(top_1$Direction, as.numeric(top_1$adj.Pval)), ]
  top_1[duplicated(top_1[, 1]), 1] <- ""
  rownames(top_1) <- seq(1, nrow(top_1), 1)

  return(top_1)
}

#' Pathway analysis with the PGSEA package
#'
#' Run pathway analysis with the PGSEA package using the results
#' from the limma_value function.
#'
#' @param processed_data Matrix of gene data that has been through
#'   \code{\link{pre_process}()}
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database. See returned list from \code{\link{read_gene_sets}()}
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of significant pathways to show
#'
#' @export
#' @return A list with a data frame and a numeric value that is used
#'  in the \code{\link{plot_pgsea}()} to create a heatmap.
#'
#' @family pathway functions
pgsea_data <- function(processed_data,
                       gene_sets,
                       my_range,
                       pathway_p_val_cutoff,
                       n_pathway_show) {
  subtype <- detect_groups(colnames(processed_data), preserve_original = TRUE)

  # Cut off to report in PGSEA. Otherwise NA
  p_value <- 0.01
  if (length(gene_sets) == 0) {
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
  if (dim(pg_results)[1] < 2) {
    return()
  }
  best <- max(abs(pg_results))

  if (length(subtype) < 4 || length(unique(subtype)) < 2 ||
    length(unique(subtype)) == dim(processed_data)[2]) {
    pg_results <- pg_results[order(-apply(pg_results, 1, sd)), ]
    return(list(pg_data = pg_results[1:n_pathway_show, ], best <- best))
  }

  cat("\nComputing P values using ANOVA\n")
  path_p_value <- function(k,
                           pg_results,
                           subtype) {
    return(summary(aov(pg_results[k, ] ~ subtype))[[1]][["Pr(>F)"]][1])
  }
  p_values <- sapply(1:dim(pg_results)[1], function(x) {
    path_p_value(
      k = x,
      pg_results = pg_results,
      subtype = subtype
    )
  })
  p_values <- stats::p.adjust(p_values, "fdr")


  if (sort(p_values)[2] > pathway_p_val_cutoff) {
    return(list(pg_data = NULL, best = best))
  } else {
    n_sig_t <- rowSums(pg$p.results < p_value)

    result <- cbind(as.matrix(p_values), n_sig_t, pg_results)
    result <- result[order(result[, 1]), ]
    result <- result[which(result[, 1] < pathway_p_val_cutoff), , drop = F]

    pg_results <- result[, -2]

    # When there is only 1 left in the matrix pg_results becomes a vector
    if (sum(p_values < pathway_p_val_cutoff) == 1) {
      pg_data <- t(as.matrix(pg_results))
      pg_data <- rbind(pg_data, pg_data)
    } else {
      if (dim(pg_results)[1] > n_pathway_show) {
        pg_data <- pg_results[1:n_pathway_show, ]
      } else {
        pg_data <- pg_results
      }
    }

    rownames(pg_data) <- sapply(rownames(pg_data), extract_under)
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

#' Heatmap of PGSEA pathway analysis
#'
#' Create a heatmap from the pathway analysis using the PGSEA
#' package. The heatmap shows the expression in each group for
#' each significantly enriched pathway.
#'
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param processed_data Matrix of gene data that has been through
#'   \code{\link{pre_process}()}
#' @param contrast_samples Sample columns that correspond to the
#'  selected comparison
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database. See list returned from \code{\link{read_gene_sets}()}
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of significant pathways to show
#' @param select_go pathway category.
#' @param show_pathway_id Whether to show pathway id for GO and KEGG pathways
#' @param plot_colors A vector of colors for activated/surpressed pathways
#' @export
#' @return A heatmap plot with the rows as the significant
#'  pathways and the columns corresponding to the samples.
#'
#' @family pathway functions
plot_pgsea <- function(my_range,
                       processed_data,
                       contrast_samples,
                       gene_sets,
                       pathway_p_val_cutoff,
                       n_pathway_show,
                       select_go,
                       show_pathway_id,
                       margin = c(3, 1, 13, 38),
                       plot_colors = NULL) {
  genes <- processed_data[, contrast_samples]
  if (length(gene_sets) == 0) {
    return(
      NULL
    )
  } else {
    subtype <- detect_groups(colnames(genes), preserve_original = TRUE)
    result <- pgsea_data(
      processed_data = genes,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show
    )

    if (is.null(result$pg_data)) {
      plot.new()
      text(0.5, 1, "No significant pathway found!")
    } else {

       # remove pathway ID if selected so
      if (!show_pathway_id) {
        row.names(result$pg_data) <- remove_pathway_id_second(
          strings = row.names(result$pg_data),
          select_go = select_go
        )
      }
      
      if(is.null(plot_colors)){
        color_vec <- PGSEA::.rwb
      }
      else{
        color_vec_1 <- colorRampPalette(c(plot_colors[[1]][1],"white"))(25)
        color_vec_2 <- colorRampPalette(c("white",plot_colors[[1]][2]))(25)
        color_vec <- c(color_vec_1,color_vec_2)
      }
      
      PGSEA::smcPlot(
        result$pg_data,
        factor(subtype),
        scale = c(-max(result$pg_data), max(result$pg_data)),
        show.grid = T,
        margins = margin,
        col = color_vec,
        cex.lab = 0.5
      )
    }
  }
}

#' Heatmap of GSVA pathway analysis
#'
#' Create a heatmap from the pathway analysis using the GSVA
#' package. The heatmap shows the expression in each group for
#' each significantly enriched pathway.
#'
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param processed_data Matrix of gene data that has been through
#'   \code{\link{pre_process}()}
#' @param contrast_samples Sample columns that correspond to the
#'  selected comparison
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database. See list returned from \code{\link{read_gene_sets}()}
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of significant pathways to show
#' @param select_go pathway category.
#' @param show_pathway_id Whether to show pathway id for GO and KEGG pathways
#'
#' @export
#' @return A heatmap plot with the rows as the significant
#'  pathways and the columns corresponding to the samples.
#'
#' @family pathway functions
plot_gsva <- function(my_range,
                       processed_data,
                       contrast_samples,
                       gene_sets,
                       pathway_p_val_cutoff,
                       n_pathway_show,
                       select_go,
                       show_pathway_id,
                       algorithm = "gsva") {
  genes <- processed_data[, contrast_samples]
  if (length(gene_sets) == 0) {
    plot.new()
    text(0.5, 1, "No significant pathway found!")
  } else {
    subtype <- detect_groups(colnames(genes), preserve_original = TRUE)
    result <- gsva_data(
      processed_data = genes,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show,
      algorithm = algorithm
    )

    # remove notificatoin from last analysis
    removeNotification("small_sample_size")
    if(ncol(genes) <= 10) {
      showNotification(
        ui = paste("Only ", ncol(genes), "samples! GSVA results are not reliable when sample sizes are small (<=10)"),
        id = "small_sample_size",
        duration = NULL,
        type = "error"
      )

    }

    if (is.null(result$pg_data)) {
      plot.new()
      text(0.5, 1, "No significant pathway found!")
    } else {

       # remove pathway ID if selected so
      if (!show_pathway_id) {
        row.names(result$pg_data) <- remove_pathway_id_second(
          strings = row.names(result$pg_data),
          select_go = select_go
        )
      }

      if(algorithm == "ssgsea") {
        result$pg_data <- result$pg_data - rowMeans(result$pg_data)        
      }

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


#' Pathway analysis with the GSVA package
#'
#' Run pathway analysis with the GSVA package using the results
#' from the limma_value function.
#'
#' @param processed_data Matrix of gene data that has been through
#'   \code{\link{pre_process}()}
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database. See returned list from \code{\link{read_gene_sets}()}
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of significant pathways to show
#' @param algorithm  Options for GSVA: plage, ssgsea, zscore or gsva
#'
#' @export
#' @return A list with a data frame and a numeric value that is used
#'  in the \code{\link{plot_gsva}()} to create a heatmap.
#'
#' @family pathway functions
gsva_data <- function(processed_data,
                       gene_sets,
                       my_range,
                       pathway_p_val_cutoff,
                       n_pathway_show,
                       algorithm = "gsva") {
  subtype <- detect_groups(colnames(processed_data), preserve_original = TRUE)

  # Cut off to report in PGSEA. Otherwise NA
  p_value <- 0.01
  if (length(gene_sets) == 0) {
    return(list(pg3 = NULL, best = 1))
  }
  
  # # Modern Syntax for current GSVA versions
  if (algorithm == "gsva"){
   
     param <- GSVA::gsvaParam(processed_data, gene_sets)
   
  } else if (algorithm == "ssgsea") {
   
    param <- GSVA::ssgseaParam(processed_data, gene_sets)
   
  } else if (algorithm == "plage") {
   
    param <- GSVA::plageParam(processed_data, gene_sets)
  }
   
  pg_results <- GSVA::gsva(param = param, verbose = FALSE)
   
 # Deprecated syntax for old versions
 #pg_results <- GSVA::gsva(processed_data, gene_sets, verbose = FALSE, method = algorithm)

  # Remove se/wrts with all missing(non-signficant)
  pg_results <- pg_results[rowSums(is.na(pg_results)) < ncol(pg_results), ]
  if (dim(pg_results)[1] < 2) {
    return()
  }
  best <- max(abs(pg_results))

  if (length(subtype) < 4 || length(unique(subtype)) < 2 ||
    length(unique(subtype)) == dim(processed_data)[2]) {
    pg_results <- pg_results[order(-apply(pg_results, 1, sd)), ]
    return(list(pg_data = pg_results[1:top, ], best <- best))
  }

  cat("\nComputing P values using ANOVA\n")
  path_p_value <- function(k,
                           pg_results,
                           subtype) {
    return(summary(aov(pg_results[k, ] ~ subtype))[[1]][["Pr(>F)"]][1])
  }
  p_values <- sapply(1:dim(pg_results)[1], function(x) {
    path_p_value(
      k = x,
      pg_results = pg_results,
      subtype = subtype
    )
  })
  p_values <- stats::p.adjust(p_values, "fdr")


  if (sort(p_values)[2] > pathway_p_val_cutoff) {
    return(list(pg_data = NULL, best = best))
  } else {
    result <- cbind(as.matrix(p_values), pg_results)
    result <- result[order(result[, 1]), ]
    result <- result[which(result[, 1] < pathway_p_val_cutoff), , drop = FALSE]

    pg_results <- result

    # When there is only 1 left in the matrix pg_results becomes a vector
    if (sum(p_values < pathway_p_val_cutoff) == 1) {
      pg_data <- t(as.matrix(pg_results))
      pg_data <- rbind(pg_data, pg_data)
    } else {
      if (dim(pg_results)[1] > n_pathway_show) {
        pg_data <- pg_results[1:n_pathway_show, ]
      } else {
        pg_data <- pg_results
      }
    }

    rownames(pg_data) <- sapply(rownames(pg_data), extract_under)
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

#' Data from GSVA plot
#'
#' Get the data matrix that is plotted in the heatmap created by
#' the \code{\link{plot_pgsea}()}.
#'
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param data Matrix of gene data that has been through
#'  \code{\link{pre_process}()}
#' @param select_contrast String designating the comparison from DEG analysis to
#'  filter for the significant genes. See the 'comparison' element from the list
#'  returned from \code{\link{limma_value}()} for options.
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database. See list returned from \code{\link{read_gene_sets}()}
#' @param sample_info Matrix of experiment design information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' @param select_model_comprions String designating selected comparisons to
#'  analyze in the DEG analysis. See \code{\link{list_model_comparisons_ui}()}
#'  for options
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#' @param algorithm  Options for GSVA: plage, ssgsea, zscore or gsva
#'  result
#'
#' @export
#' @return Data matrix with the rownames the descriptions of pathways
#'  and the matrix the returned expression calculation from the PGSEA
#'  package.
#'
#' @family pathway functions
get_gsva_plot_data <- function(my_range,
                                data,
                                select_contrast,
                                gene_sets,
                                sample_info,
                                select_factors_model,
                                select_model_comprions,
                                pathway_p_val_cutoff,
                                n_pathway_show,
                                algorithm = "gsva") {
  # Find sample related to the comparison
  iz <- match(
    detect_groups(colnames(data), preserve_original = TRUE),
    unlist(strsplit(select_contrast, "-"))
  )
  iz <- which(!is.na(iz))

  if (!is.null(sample_info) & !is.null(select_factors_model) & length(select_model_comprions) > 0) {
    # Strings like: "groups: mutant vs. control"
    comparisons <- gsub(".*: ", "", select_model_comprions)
    comparisons <- gsub(" vs\\. ", "-", comparisons)
    # Corresponding factors
    factors_vector <- gsub(":.*", "", select_model_comprions)
    # Selected contrast lookes like: "mutant-control"
    ik <- match(select_contrast, comparisons)
    if (is.na(ik)) {
      iz <- 1:(dim(data)[2])
    } else {
      # Interaction term, use all samples
      # Corresponding factors
      selected_factor <- factors_vector[ik]
      iz <- match(sample_info[, selected_factor], unlist(strsplit(select_contrast, "-")))
      iz <- which(!is.na(iz))
    }
  }

  if (grepl("I:", select_contrast)) {
    # If it is factor design use all samples
    iz <- 1:(dim(data)[2])
  }
  if (is.na(iz)[1] | length(iz) <= 1) {
    iz <- 1:(dim(data)[2])
  }

  genes <- data
  genes <- genes[, iz]

  subtype <- detect_groups(colnames(genes), preserve_original = TRUE)

  if (length(gene_sets) == 0) {
    return(as.data.frame("No significant pathway!"))
  } else {
    result <- gsva_data(
      processed_data = genes,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show,
      algorithm = algorithm
    )

    if (is.null(result$pg_data)) {
      return(as.data.frame("No significant pathway!"))
    } else {
      return(as.data.frame(result$pg_data))
    }
  }
}


#' Pathway analysis with the FGSEA package
#'
#' Run pathway analysis with the FGSEA package using the results
#' from the \code{\link{limma_value}()}.
#'
#' @param select_contrast String designating the comparison from DEG analysis to
#'  filter for the significant genes. See the 'comparison' element from the list
#'  returned from \code{\link{limma_value}()} for options.
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param limma Results list from the \code{\link{limma_value}()}
#' @param gene_p_val_cutoff Significant p-value to filter
#'  the top genes fold change by
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database. See results list from
#'  \code{\link{read_gene_sets}()}.
#' @param absolute_fold TRUE/FALSE to use the absolute value of the fold
#'  change
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#'  result
#'
#' @export
#' @return A data frame with the results of the pathway analysis.
#'  The data frame has five columns for the direction of the
#'  regulation, the pathway description, the stat value, the
#'  number of overlapping genes, and the p-value.
#'
#' @family pathway functions
fgsea_data <- function(select_contrast,
                       my_range,
                       limma,
                       gene_p_val_cutoff,
                       gene_sets,
                       absolute_fold,
                       pathway_p_val_cutoff,
                       n_pathway_show) {
  nPerm <- 10000 # number of permutations

  no_sig <- as.data.frame("No significant pathway found.")
  if (length(limma$top_genes) == 0) {
    return(no_sig)
  }
  if (length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]
  } else {
    top <- limma$top_genes
    ix <- match(select_contrast, names(top))
    if (is.na(ix)) {
      return(no_sig)
    }
    top_1 <- top[[ix]]
  }
  if (dim(top_1)[1] == 0) {
    return(no_sig)
  }
  colnames(top_1) <- c("Fold", "FDR")

  # Remove some genes
  top_1 <- top_1[which(top_1$FDR < gene_p_val_cutoff), ]

  if (length(gene_sets) == 0) {
    return(as.data.frame("No gene set found!"))
  }

  fold <- top_1[, 1]
  names(fold) <- rownames(top_1)

  # Use absolute value of fold change, disregard direction
  if (absolute_fold) {
    fold <- abs(fold)
  }

  # nproc = 0 lets fgsea auto-detect available cores. On Windows, BiocParallel
  # falls back to serial (1 core) automatically because SOCK clusters fail to
  # serialize large gene sets. On Linux/Mac, fork-based parallelism is used.
  paths <- fgsea::fgsea(
    pathways = gene_sets,
    stats = fold,
    minSize = my_range[1],
    maxSize = my_range[2],
    nPerm = nPerm,
    nproc = 0
  )

  if (dim(paths)[1] < 1) {
    return(no_sig)
  }
  paths <- as.data.frame(paths)
  # Sort by NES
  paths <- paths[order(-abs(paths[, 5])), ]
  top_1 <- paths[, c(1, 5, 7, 3)]
  colnames(top_1) <- c("Pathway", "NES", "Genes", "adj.Pval")

  if (length(which(top_1[, 4] <= pathway_p_val_cutoff)) == 0) {
    return(no_sig)
  }
  top_1 <- top_1[which(top_1[, 4] <= pathway_p_val_cutoff), , drop = FALSE]

  if (dim(top_1)[1] > n_pathway_show) {
    top_1 <- top_1[1:n_pathway_show, , drop = FALSE]
  }

  top_1 <- as.data.frame(top_1)
  top_1 <- cbind(rep(select_contrast, dim(top_1)[1]), top_1)
  top_1[, 4] <- as.character(round(as.numeric(top_1[, 4]), 4))
  top_1$adj.Pval <- sprintf("%-2.1e", as.numeric(top_1$adj.Pval))
  top_1[, 1] <- as.character(top_1[, 1])
  colnames(top_1)[1] <- "Direction"
  colnames(top_1)[2] <- paste("GSEA analysis:", gsub("-", " vs ", select_contrast))
  top_1[which(as.numeric(top_1[, 3]) > 0), 1] <- "Up"
  top_1[which(as.numeric(top_1[, 3]) < 0), 1] <- "Down"
  # sort by adj.Pval
  top_1 <- top_1[order(top_1[, 1], as.numeric(top_1$adj.Pval)), ]
  top_1[duplicated(top_1[, 1]), 1] <- ""
  top_1[, 3] <- as.character(round(as.numeric(top_1[, 3]), 4))

  return(top_1)
}

#' Pathway analysis with reactome package
#'
#' Run pathway analysis with the reactome package using the results
#' from the limma_value function.
#'
#' @param select_contrast String designating the comparison from DEG analysis to
#'  filter for the significant genes. See the 'comparison' element from the list
#'  returned from \code{\link{limma_value}()} for options.
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param limma Results list from the \code{\link{limma_value}()}
#' @param gene_p_val_cutoff Significant p-value to filter
#'  the top genes fold change by
#' @param converted Return value from the \code{\link{convert_id}()} function,
#'  contains information about the gene IDs for the matched species
#' @param idep_data List of data files from the database, returned from
#'  \code{\link{get_idep_data}()}
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#'  result
#' @param absolute_fold TRUE/FALSE to use the absolute value of the fold
#'  change
#'
#' @export
#' @return A data frame with the results of the pathway analysis.
#'  The data frame has five columns for the direction of the
#'  regulation, the pathway description, the stat value, the
#'  number of overlapping genes, and the p-value.
#'
#' @family pathway functions
#'
reactome_data <- function(select_contrast,
                          my_range,
                          limma,
                          gene_p_val_cutoff,
                          converted,
                          idep_data,
                          pathway_p_val_cutoff,
                          n_pathway_show,
                          absolute_fold) {
  
  ensembl_species <- c(
    "hsapiens_gene_ensembl", "rnorvegicus_gene_ensembl", "mmusculus_gene_ensembl",
    "celegans_gene_ensembl", "scerevisiae_gene_ensembl", "drerio_gene_ensembl",
    "dmelanogaster_gene_ensembl"
  )
  reactome_pa_species <- c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly")
  no_sig <- as.data.frame("No significant pathway found.")
  if (length(limma$top_genes) == 0) {
    return(no_sig)
  }
  if (length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]
  } else {
    top <- limma$top_genes
    ix <- match(select_contrast, names(top))
    if (is.na(ix)) {
      return(no_sig)
    }
    top_1 <- top[[ix]]
  }
  if (dim(top_1)[1] == 0) {
    return(no_sig)
  }
  colnames(top_1) <- c("Fold", "FDR")

  # Remove some genes
  top_1 <- top_1[which(top_1$FDR < gene_p_val_cutoff), ]

  fold <- top_1[, 1]
  names(fold) <- rownames(top_1)
  if (absolute_fold) {
    # Use absolute value of fold change, disregard direction
    fold <- abs(fold)
  }

  species <- converted$species[1, 1]
  ix <- match(species, ensembl_species)

  if (is.na(ix)) {
    return(as.data.frame("Species not coverted by ReactomePA package!"))
  }

  fold <- convert_ensembl_to_entrez(
    query = fold,
    species = species,
    org_info = idep_data$org_info,
    idep_data = idep_data
  )
  # Remove duplicate gene entrez IDs
  fold <- fold[!duplicated(names(fold))]
  
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

  if (is.null(paths)) {
    return(no_sig)
  }
  if (dim(paths)[1] == 0) {
    return(no_sig)
  }

  if (dim(paths)[1] < 1) {
    return(no_sig)
  }
  paths <- as.data.frame(paths)
  paths <- paths[order(-abs(paths[, 5])), ]

  top_1 <- paths[, c(2, 5, 3, 7)]

  colnames(top_1) <- c("Pathway", "NES", "Genes", "adj.Pval")

  if (length(which(top_1[, 4] <= pathway_p_val_cutoff)) == 0) {
    return(no_sig)
  }
  top_1 <- top_1[which(top_1[, 4] <= pathway_p_val_cutoff), , drop = FALSE]

  if (dim(top_1)[1] > n_pathway_show) {
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
    gsub("-", " vs ", select_contrast)
  )
  top_1[which(as.numeric(top_1[, 3]) > 0), 1] <- "Up"
  top_1[which(as.numeric(top_1[, 3]) < 0), 1] <- "Down"
  top_1 <- top_1[order(top_1[, 1], -abs(as.numeric(top_1[, 3]))), ]
  top_1[duplicated(top_1[, 1]), 1] <- ""
  top_1[, 3] <- as.character(round(as.numeric(top_1[, 3]), 4))

  return(top_1)
}

#' Pathway analysis with the PGSEA package on all samples
#'
#' Run pathway analysis with the PGSEA package using the results
#' from the \code{\link{limma_value}()} on all samples.
#'
#' @param go  String designating the section of the database to query for
#'   pathway analysis. See \code{\link{gmt_category}()} for choices.
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param data Matrix of gene data that has been through the
#'   \code{\link{pre_process}()}
#' @param select_contrast String designating the comparison from DEG analysis to
#'  filter for the significant genes. See the 'comparison' element from the list
#'  returned from \code{\link{limma_value}()} for options.
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database. See list returned from \code{\link{read_gene_sets}()}
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#'  result
#' @param select_go pathway category.
#' @param show_pathway_id Whether to show pathway id for GO and KEGG pathways
#' @param plot_colors A vector of colors for activated/surpressed pathways
#' 
#' @export
#' @return A data frame with the results of the pathway analysis.
#'  The data frame has five columns for the direction of the
#'  regulation, the pathway description, the stat value, the
#'  number of overlapping genes, and the p-value.
#'
#' @family pathway functions
pgsea_plot_all <- function(go,
                           my_range,
                           data,
                           select_contrast,
                           gene_sets,
                           pathway_p_val_cutoff,
                           n_pathway_show,
                           select_go,
                           show_pathway_id,
                           margin = c(3, 1, 13, 38),
                           plot_colors = NULL) {
  if (length(gene_sets) == 0) {
    plot.new()
    text(0, 1, "No gene sets!")
  } else {
    subtype <- detect_groups(colnames(data), preserve_original = TRUE)
    result <- pgsea_data(
      processed_data = data,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show
    )

    if (is.null(result$pg_data)) {
      plot.new()
      text(0.5, 1, "No significant pathway found!")
    } else {

       # remove pathway ID if selected so
      if (!show_pathway_id) {
        row.names(result$pg_data) <- remove_pathway_id_second(
          strings = row.names(result$pg_data), 
          select_go = select_go
        )
      }
      
      if(is.null(plot_colors)){
        color_vec <- PGSEA::.rwb
      }
      else{
        color_vec_1 <- colorRampPalette(c(plot_colors[[1]][1],"white"))(25)
        color_vec_2 <- colorRampPalette(c("white",plot_colors[[1]][2]))(25)
        color_vec <- c(color_vec_1,color_vec_2)
      }

      PGSEA::smcPlot(
        result$pg_data,
        factor(subtype),
        scale = c(-max(result$pg_data), max(result$pg_data)),
        show.grid = T,
        margins = margin,
        col = color_vec,
        cex.lab = 0.5
      )
    }
  }
}

#' Data from PGSEA plot
#'
#' Get the data matrix that is plotted in the heatmap created by
#' the \code{\link{plot_pgsea}()}.
#'
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param data Matrix of gene data that has been through
#'  \code{\link{pre_process}()}
#' @param contrast_samples Sample columns that correspond to the
#'  selected comparison
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database. See list returned from \code{\link{read_gene_sets}()}
#' @param sample_info Matrix of experiment design information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' @param select_model_comprions String designating selected comparisons to
#'  analyze in the DEG analysis. See \code{\link{list_model_comparisons_ui}()}
#'  for options
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#'  result
#'
#' @export
#' @return Data matrix with the rownames the descriptions of pathways
#'  and the matrix the returned expression calculation from the PGSEA
#'  package.
#'
#' @family pathway functions
get_pgsea_plot_data <- function(my_range,
                                data,
                                contrast_samples,
                                gene_sets,
                                sample_info,
                                select_factors_model,
                                select_model_comprions,
                                pathway_p_val_cutoff,
                                n_pathway_show) {
  genes <- data[, contrast_samples]

  subtype <- detect_groups(colnames(genes), preserve_original = TRUE)

  if (length(gene_sets) == 0) {
    return(as.data.frame("No significant pathway!"))
  } else {
    result <- pgsea_data(
      processed_data = genes,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show
    )

    if (is.null(result$pg_data)) {
      return(as.data.frame("No significant pathway!"))
    } else {
      return(as.data.frame(result$pg_data))
    }
  }
}

#' Data from PGSEA plot all samples
#'
#' Get the data matrix that is plotted in the heatmap created by
#' the pgsea_plot_all function.
#'
#' @param data Matrix of gene data that has been through
#'  \code{\link{pre_process}()}
#' @param select_contrast String designating the comparison from DEG analysis to
#'  filter for the significant genes. See the 'comparison' element from the list
#'  returned from \code{\link{limma_value}()} for options.
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database. See list returned from \code{\link{read_gene_sets}()}
#' @param my_range Vector of the (min_set_size, max_set_size)
#' @param pathway_p_val_cutoff Significant p-value to determine
#'  enriched pathways
#' @param n_pathway_show Number of pathways to return in final
#'  result
#'
#' @export
#' @return Data matrix with the rownames the descriptions of pathways
#'  and the matrix the returned expression calculation from the PGSEA
#'  package.
#'
#' @family pathway functions
get_pgsea_plot_all_samples_data <- function(data,
                                            select_contrast,
                                            gene_sets,
                                            my_range,
                                            pathway_p_val_cutoff,
                                            n_pathway_show) {
  genes <- data
  subtype <- detect_groups(colnames(genes), preserve_original = TRUE)

  if (length(gene_sets) == 0) {
    plot.new()
    text(0, 1, "No gene sets!")
  } else {
    result <- pgsea_data(
      processed_data = genes,
      gene_sets = gene_sets,
      my_range = my_range,
      pathway_p_val_cutoff = pathway_p_val_cutoff,
      n_pathway_show = n_pathway_show
    )

    if (is.null(result$pg_data)) {
      return(as.data.frame("No significant pathway!"))
    } else {
      return(as.data.frame(result$pg_data))
    }
  }
}

#' Transform Pathway Data
#' 
#' Transform data from various pathway methods for uniform download format
#'
#' @param data Pathway data from selected method
#' @param contrast Selected contrast from Stats tab
#' @param method  Selected pathway method
#' @param genes Gene data
#' @param org  Selected org from pre-process step
#' @param path_id Show pathway ID toggle
#' @param go Selected pathway database
#' @param deg Up/Down-regulation results from Stats
#'
#' @return Data frame with hyperlinks and urls for pathways
#' 
#' @export
#'          
pathway_data_transform <- function(data,
                                   contrast,
                                   method,
                                   genes,
                                   org,
                                   path_id,
                                   go,
                                   deg){
  
  if (data[1,1] == "No significant pathway found."){
    return(data)
  }
  
  if (method %in% c("PGSEA", "GSVA", "ssGSEA", "PLAGE")){
    rn <- rownames(data)
    
    data <- data.frame(
      adj.Pval = sub("^([0-9.eE+-]+)\\s+.*", "\\1", rn),
      pathway = sub("^[0-9.eE+-]+\\s+", "", rn),
      data,
      row.names = NULL,
      check.names = FALSE
    )
  }
  
  if(method == "ssGSEA") {
    data[,-c(1:3)] <- data[, -c(1:3)] - rowMeans(data[, -c(1:3)])        
  }
  
  # Name for hypertext column
  hypertext_name <- paste0(method," Analysis: ", contrast)
  # Rename 2nd column
  colnames(data)[2] <- paste(hypertext_name, "Pathways")
  
  if (ncol(data) > 1) {
    # add URL
    ix <- match(data[, 2], genes$pathway_info$description)
    
    paths <- genes$gene_lists[which(names(genes$gene_lists) %in% data[,2])]
    # Find gene matches in top pathways
    path_match <- lapply(rownames(deg), function(x) {
      groups <- names(paths)[sapply(paths, function(vec) x %in% vec)]
      if (length(groups) > 0) {
        data.frame(Gene = x, group = groups, stringsAsFactors = FALSE)
      } else {
        NULL  # skip this element if no groups matched
      }
    })
    
    # Combine into one data frame
    result_df <- do.call(rbind, path_match)

    # Extract the specific contrast column from deg results matrix
    if (is.matrix(deg) || is.data.frame(deg)) {
      # If deg has multiple columns (comparisons), select the one matching contrast
      if (!is.null(contrast) && contrast %in% colnames(deg)) {
        deg_contrast <- data.frame(Gene = rownames(deg),
                                   Expr = deg[, contrast],
                                   stringsAsFactors = FALSE)
      } else if (ncol(deg) == 1) {
        # If only one column, use it
        deg_contrast <- data.frame(Gene = rownames(deg),
                                   Expr = deg[, 1],
                                   stringsAsFactors = FALSE)
      } else {
        # Default to first column if contrast not found
        deg_contrast <- data.frame(Gene = rownames(deg),
                                   Expr = deg[, 1],
                                   stringsAsFactors = FALSE)
      }
    } else {
      # If deg is a vector
      deg_contrast <- data.frame(Gene = rownames(deg),
                                 Expr = deg,
                                 stringsAsFactors = FALSE)
    }

    counts <- merge(x = result_df, y = deg_contrast, by = "Gene", all.x = TRUE)
    # Find Up/Down Gene count
    counts <- dplyr::group_by(counts, group) |>
      dplyr::summarize(UpGenes = sum(Expr == 1, na.rm = TRUE),
                       DownGenes = sum(Expr == -1, na.rm = TRUE),
                       UnregGenes = sum(Expr == 0, na.rm = TRUE))
    
    data2 <- merge(x = data, y = counts, by.x = colnames(data)[2], by.y = "group")
    # sort by adj.Pval
    data2 <- data2[match(data[,2], data2[,1]), ]
    
    # remove pathway ID, but only in Ensembl species
    if (!path_id && org > 0) {
      data2[, 1] <- remove_pathway_id(data2[, 1], go)
    }
    
    # Add hypertext to the end of the data
    data2[hypertext_name] <- hyperText(
      data2[, 1],
      genes$pathway_info$memo[ix]
    )
    
    # create separate URL column for download
    data2 <- data.frame(data2[,c(2,1)],
                        URL = genes$pathway_info$memo[ix],
                        data2[,-c(1,2)],
                        check.names = FALSE
    )
    
    if (method %in% c("GSEA", "GAGE")){
      data2[, 7:9] <- lapply(data2[, 7:9], as.character)
    }
    
  }
  
  return(data2)
  
}

#' Get data from genes in selected pathway
#'
#' Return a data matrix that is a subset of the processed data and
#' only contains genes that are in the gene set of the desired
#' pathway.
#'
#' @param sig_pathways Description of the pathway for which to
#'  obtain the gene expression data
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database. See list returned from
#'  \code{\link{read_gene_sets}()}
#' @param contrast_samples Vector of sample columns that correspond to the
#'  selected comparison
#' @param data Matrix of gene data that has been through
#'  \code{\link{pre_process}()}
#' @param select_org String designating the organism being analyzed
#' @param all_gene_names Matrix of all the matched and converted
#'  gene IDs from \code{\link{get_all_gene_names}()}
#' @param deg data frame of Stats gene regulation results - i.e.-1, 0 , 1 for
#' every gene
#' @param select_contrast Selected comparison/contrast name to extract regulation
#' status from (optional, defaults to first column if deg has multiple columns)
#'
#' @export
#' @return Sub-data matrix from the processed data. Only contains
#'  genes from the selected pathway and samples that correspond to
#'  the comparison being analyzed.
#'
#' @family pathway functions
pathway_select_data <- function(sig_pathways,
                                gene_sets,
                                contrast_samples,
                                data,
                                select_org,
                                all_gene_names,
                                deg,
                                select_contrast = NULL) {
  if (sig_pathways == "All") {
    return(NULL)
  }

  # Find the gene set
  ix <- which(names(gene_sets) == sig_pathways)
  if (length(ix) == 0) {
    return(NULL)
  }
  # Retrieve genes
  genes <- gene_sets[[ix]]

  # Find related samples
  iz <- contrast_samples
  x <- data[which(rownames(data) %in% genes), iz]

  # Extract the specific contrast column from deg results matrix
  if (is.matrix(deg) || (is.data.frame(deg) && ncol(deg) > 1)) {
    # If deg has multiple columns (comparisons), select the one matching select_contrast
    if (!is.null(select_contrast) && select_contrast %in% colnames(deg)) {
      deg_col <- data.frame(Regulation = deg[, select_contrast], row.names = rownames(deg))
    } else {
      # Default to first column if select_contrast not found or not provided
      deg_col <- data.frame(Regulation = deg[, 1], row.names = rownames(deg))
    }
  } else {
    # If deg has only one column or is a vector
    deg_col <- data.frame(Regulation = deg[, 1], row.names = rownames(deg))
  }

  x <- merge(x = x, y = deg_col, by = "row.names")
  rownames(x) <- x$Row.names
  x <- dplyr::select(x,-1)
  x[,ncol(x)] <- dplyr::case_when(x[,ncol(x)] == 1 ~ "Up",
                                  x[,ncol(x)] == -1 ~ "Down",
                                  TRUE ~ "None")

  return(x)
}

#' Find list of genes in Reactome Pathway
#'
#' @param sig_pathway Pathway selection
#' @param data ReactomePA pathway dataset
#' @param gene_info Gene dataset
#' @param converted pre_process converted data
#'
#' @returns data frame of genes from the selected pathway
#' @export
#'
reactome_gene_list <- function(sig_pathway,
                               data,
                               gene_info,
                               converted){
  
  ensembl_species <- c(
    "hsapiens_gene_ensembl", "rnorvegicus_gene_ensembl", 
    "mmusculus_gene_ensembl", "celegans_gene_ensembl", 
    "scerevisiae_gene_ensembl", "drerio_gene_ensembl",
    "dmelanogaster_gene_ensembl"
  )
  reactome_pa_species <- c("human", "rat", "mouse", 
                           "celegans", "yeast", "zebrafish", 
                           "fly")
  # Find species of data
  species <- converted$species[1, 1]
  #Match to ensembl ID
  ix <- match(species, ensembl_species)
  
  # Find genes in pathway selected
  gene_search <- ReactomePA::viewPathway(sig_pathway,
                                         organism = reactome_pa_species[ix],
                                         keyType = "ENSEMBLID")
  gene_search <- as.vector(gene_search$data$name)
  
  # Match genes by symbol
  symbolID <- gene_info[, c("ensembl_gene_id", "symbol")]
  symbolID$symbol <- gsub(" ", "", symbolID$symbol)
  ix <- match(gene_search, symbolID$symbol)
  
  genes <- symbolID[ix, 1]
  
  # search for genes in submitted data
  x <- data[which(rownames(data) %in% genes),]
  
  return(x)
  
}

#' Create pathway table with gene sets
#'
#' Create a data frame of significant pathways and their analysis
#' values. Also add a column that contains the gene sets for the
#' pathway.
#'
#' @param pathway_method Integer indicating which pathway method to use. Should
#'  be one of 1 for "GAGE", 2 = "PGSEA", 3 for "GSEA (preranded fgsea)", 4
#'  for "PGSEA w/ all samples", and 5 for "ReactomePA".
#' @param gage_pathway_data Matrix returned from \code{\link{gage_data}()}
#' @param fgsea_pathway_data Matrix returned from \code{\link{fgsea_data}()}
#' @param pgsea_plot_data Matrix returned from
#'  \code{\link{get_pgsea_plot_data}()}
#' @param pgsea_plot_all_samples_data Matrix returned from
#'  \code{\link{get_pgsea_plot_all_samples_data}()}
#' @param gsva_plot_data Matrix returned from the \code{\link{gsva_plot_data}()}
#'  function
#' @param reactome_pa_pathway_data Data frame returned from
#'  \code{\link{reactome_data}()}
#' @param go String designating the section of the database to query for
#'   pathway analysis. See \code{\link{gmt_category}()} for choices.
#' @param select_org String designating which organism is being analyzed
#' @param gene_info Matrix returned from \code{\link{gene_info}()} function, all
#'  gene info from the database query with the User gene IDs
#' @param gene_sets List of vectors with each vector being the
#'  set of genes that correspond to a particular pathway in
#'  the database \code{\link{read_gene_sets function}()}
#' @param show_pathway_id whether to show pathway id or remove it
#'
#' @export
#' @return A data frame with the pathway analysis statistics and
#'  the gene sets for each significantly enriched pathway.
#'
#' @family pathway functions
get_pathway_list_data <- function(pathway_method,
                                  gage_pathway_data,
                                  fgsea_pathway_data,
                                  pgsea_plot_data,
                                  pgsea_plot_all_samples_data,
                                  gsva_plot_data,
                                  reactome_pa_pathway_data,
                                  go,
                                  select_org,
                                  gene_info,
                                  gene_sets,
                                  show_pathway_id) {
  pathways <- NULL
  if (pathway_method == 1) {
    if (!is.null(gage_pathway_data)) {
      if (dim(gage_pathway_data)[2] > 1) {
        pathways <- gage_pathway_data
        colnames(pathways)[2] <- "Pathways"
        colnames(pathways)[4] <- "nGenes"
      }
    }
  }
  if (pathway_method == 3) {
    if (!is.null(fgsea_pathway_data)) {
      if (dim(fgsea_pathway_data)[2] > 1) {
        pathways <- fgsea_pathway_data
        colnames(pathways)[2] <- "Pathways"
        colnames(pathways)[4] <- "nGenes"
      }
    }
  }
  if (pathway_method == 2) {
    if (!is.null(pgsea_plot_data)) {
      if (dim(pgsea_plot_data)[2] > 1) {
        pathways <- as.data.frame(pgsea_plot_data)
        pathways$Pathways <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
        pathways$adj.Pval <- gsub(" .*", "", rownames(pathways))
        pathways$Direction <- "Diff"
      }
    }
  }
  if (pathway_method == 4) {
    if (!is.null(pgsea_plot_all_samples_data)) {
      if (dim(pgsea_plot_all_samples_data)[2] > 1) {
        pathways <- as.data.frame(pgsea_plot_all_samples_data)
        pathways$Pathways <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
        pathways$adj.Pval <- gsub(" .*", "", rownames(pathways))
        pathways$Direction <- "Diff"
      }
    }
  }

  if (pathway_method >= 6 && pathway_method <= 8 ) {
    if (!is.null(gsva_plot_data)) {
      if (dim(gsva_plot_data)[2] > 1) {
        pathways <- as.data.frame(gsva_plot_data)
        pathways$Pathways <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
        pathways$adj.Pval <- gsub(" .*", "", rownames(pathways))
        pathways$Direction <- "Diff"
      }
    }
  }

  if (pathway_method == 5) {
    if (!is.null(reactome_pa_pathway_data)) {
      if (ncol(reactome_pa_pathway_data) >= 4) {
        pathways <- as.data.frame(reactome_pa_pathway_data)
        colnames(pathways)[2] <- "Pathways"
        colnames(pathways)[4] <- "nGenes"
      }
    }
  }

  if (is.null(pathways)) {
    return(NULL)
  }
  # if no gene set data, return pathway list
  if (is.null(gene_sets)) {
    return(pathways)
  }

  pathways$adj_p_val <- as.numeric(pathways$adj.Pval)
  pathways <- subset(pathways, select = -c(adj.Pval))
  pathways$adj_p_val <- as.character(pathways$adj_p_val)

  # Sometimes only one pathway is in the table
  if (nrow(pathways) > 1) {
    for (i in 2:nrow(pathways)) {
      if (nchar(pathways$Direction[i]) <= 1) {
        pathways$Direction[i] <- pathways$Direction[i - 1]
      }
    }
  }

  # Gene symbol matching symbols
  probe_to_gene <- NULL
  if (go != "ID not recognized!" & select_org != "NEW") {
    # If more than 50% genes has symbol
    if (sum(is.na(gene_info$symbol)) / dim(gene_info)[1] < .5) {
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

  pathways$Genes <- vector(mode = "list", length = nrow(pathways))
  # looking up genes for each pathway
  for (i in 1:nrow(pathways)) {
    # Find the gene set
    ix <- which(names(gene_sets) == pathways$Pathways[i])
    if (length(ix) != 0) {
      # Retrieve genes
      if (length(ix) > 1) {
        genes <- gene_sets[[ix[[1]]]]
      } else {
        genes <- gene_sets[[ix]]
      }

      if (!is.null(probe_to_gene)) {
        iy <- match(genes, probe_to_gene[, 1])
        genes <- probe_to_gene[iy, 2]
      }
      pathways$Genes[[i]] <- c(genes)
    }
  }

    # remove pathway ID if selected so
  if (!show_pathway_id) {
    pathways$Pathways <- remove_pathway_id(
      strings = pathways$Pathways,
      select_go = go
    )
  }

  return(pathways)
}

#' Use KEGG to create a pathway diagram
#'
#' In the database, use the KEGG information to create an image
#' that is a diagram of the pathway that is being enriched.
#'
#' @param go String designating the section of the database to query for
#'   pathway analysis. See \code{\link{gmt_category}()} for choices.
#' @param gage_pathway_data Matrix returned from \code{\link{gage_data}()}
#' @param sig_pathways Description of the pathway for which to
#'  obtain the gene expression data
#' @param select_contrast String designating the comparison from DEG analysis to
#'  filter for the significant genes. See the 'comparison' element from the list
#'  returned from \code{\link{limma_value}()} for options.
#' @param limma Results list from \code{\link{limma_value}()}
#' @param converted Return value from the \code{\link{convert_id}()}, contains
#'  information about the gene IDs for the matched species
#' @param idep_data Read data files from \code{\link{get_idep_data}()}
#' @param select_org The organism that the gene data is for
#' @param low_color String designating color value for the low-ly expressed
#'  genes
#' @param high_color String designating color value for the high-ly expressed
#'  genes
#'
#' @export
#' @return Make an image and return the path to the image to be
#'  rendered in the server.
#'
#' @family pathway functions
kegg_pathway <- function(go,
                         gage_pathway_data,
                         sig_pathways,
                         select_contrast,
                         limma,
                         converted,
                         idep_data,
                         select_org,
                         low_color = "green",
                         high_color = "red") {
  # First generate a blank image. Otherwise return(NULL) gives us errors.
  out_file <- tempfile(fileext = ".png")
  png(out_file, width = 400, height = 300)

  frame()
  dev.off()
  blank <- list(
    src = out_file,
    contentType = "image/png",
    width = 400,
    height = 300,
    alt = "Not downloaded."
  )


  if (is.null(go) || go != "KEGG") {
    return(blank)
  }
  if (is.null(gage_pathway_data)) {
    return(blank)
  }
  if (is.null(sig_pathways)) {
    return(blank)
  }

  if (is.null(select_contrast)) {
    return(blank)
  }

  if (sig_pathways == "All") {
    return(blank)
  }

  if (length(limma$top_genes) == 0) {
    return(blank)
  }

  # Get fold change
  if (length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]
  } else {
    top <- limma$top_genes
    ix <- match(select_contrast, names(top))
    if (is.na(ix)) {
      return(blank)
    }
    top_1 <- top[[ix]]
  }
  if (dim(top_1)[1] == 0) {
    return(blank)
  }


  colnames(top_1) <- c("Fold", "FDR")
  species <- converted$species[1, 1]

  fold <- top_1[, 1]
  names(fold) <- rownames(top_1)
  fold <- convert_ensembl_to_entrez(
    query = fold,
    species = species,
    org_info = idep_data$org_info,
    idep_data = idep_data
  )
  if (is.null(fold)) {
    warning("No Entrez IDs available for KEGG plotting; returning placeholder image.")
    return(blank)
  }
  fold <- fold[!is.na(fold)]
  if (!length(fold)) {
    warning("No fold-change values available for KEGG plotting; returning placeholder image.")
    return(blank)
  }


  kegg_species_id <- idep_data$kegg_species_id

  kegg_species <- as.character(
    kegg_species_id[which(kegg_species_id[, 1] == species), 3]
  )

  # look up KEGG species ID "hsa", "mmu"
  kegg_species <- as.character(
    idep_data$org_info$KEGG[which(idep_data$org_info$ensembl_dataset == species)]
  )

  if (nchar(kegg_species) <= 2) {
    return(blank)
  }

  # find pathway id
  # "Path:hsa04110 Cell cycle" --> "hsa04110"
  path_id <- gsub(" .*", "", sig_pathways)
  path_id <- gsub("Path:", "", path_id)

  if(0) {
    path_id <- kegg_pathway_id(
      sig_pathways,
      species,
      "KEGG",
      select_org,
      idep_data$gmt_files,
      idep_data$org_info,
      idep_data
    )
  }

  # Kegg pathway id not found.
  if (is.null(path_id)) {
    return(blank)
  }
  if (nchar(path_id) < 3) {
    return(blank)
  }
  random_string <- gsub(".*file", "", tempfile())
  temp_folder <- tempdir()
  out_file <- paste(
    temp_folder,
    "/",
    path_id,
    ".",
    random_string,
    ".png",
    sep = ""
  )

  pv.out <- mypathview(
    gene.data = fold,
    pathway.id = path_id,
    kegg.dir = temp_folder,
    out.suffix = random_string,
    species = kegg_species,
    kegg.native = TRUE,
    low = list(gene = low_color, cpd = "blue"),
    high = list(gene = high_color, cpd = "yellow")
  )

  # Return a list containing the filename
  list(
    src = out_file,
    contentType = "image/png",
    width = "100%",
    height = "100%",
    alt = "KEGG pathway image."
  )
}



#' Remove Pathway ID from pathway name 
#' Only for GO and KEGG pathways
#'
#' Path:hsa00270 Cysteine and methionine metabolism 
#'           --> Cysteine and methionine metabolism
#'
#' @param strings a vector of strings
#' @param select_go   GOBP, GOCC, GOMP or KEGG or something else
#'
#' @export
#' @return a vector of strings
#'
#' @family pathway functions
remove_pathway_id <- function(strings, select_go) {
    if (is.null(strings)) {
      return(NULL)
    } else {
      if (select_go %in% c("GOBP", "GOCC", "GOMF", "KEGG")) {
        strings <- sub(
          "^\\S+\\s",
          "",
          strings
        )
        strings <- proper(strings)
      }
      return(strings)
    }
}

#' Remove Pathway ID from pathway name in PGSEA
#' Only for GO and KEGG pathways
#'
#' 5.40e-05 Path:hsa04110 Cell cycle
#'     --> 5.40e-05 Cell cycle
#'
#' @param strings a vector of strings
#' @param select_go   GOBP, GOCC, GOMP or KEGG or something else
#'
#' @export
#' @return a vector of strings
#'
#' @family pathway functions
remove_pathway_id_second <- function(strings, select_go) {
    if (is.null(strings)) {
      return(NULL)
    } else {
      if (select_go %in% c("GOBP", "GOCC", "GOMF", "KEGG")) {
        FDRs <- gsub(" .*", "", strings)

        # remove FDR
        strings <- remove_pathway_id(
          strings = strings,
          select_go = select_go
        )

        # remove pathway ID
        strings <- remove_pathway_id(
          strings = strings,
          select_go = select_go
        )
        # add FDR back
        strings <- paste(FDRs, strings)
      }
      return(strings)
    }
}
