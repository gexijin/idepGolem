#' fct_06_deg1.R This file holds all of the main data analysis functions
#' associated with sixth tab of the iDEP website.
#'
#'
#' @section fct_06_deg1.R functions:
#' \code{change_names}
#'
#'
#' @name fct_06_deg1.R
NULL


# Change comparison names in limma from "KO_ko-WT_ko" to    "KO-WT_for_ko"
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param comparison DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
change_names <- function(comparison) {
  # check to see if work needs to be done
  # if no work needs to be done return input
  if (grepl(".*_.*-.*_.*", comparison)) {
    comparison_levels <- unlist(strsplit(comparison, "-"))
    comparison_levels <- unlist(strsplit(comparison_levels, "_"))
    if (length(comparison_levels) != 4) {
      comparison <- comparison_levels
    } else {
      comp_dup_index <- which(duplicated(comparison_levels))
      comparison <- paste0(
        comparison_levels[-c(comp_dup_index, comp_dup_index - 2)],
        collapse = "-"
      )
      comparison <- paste0(
        comparison,
        "_for_",
        comparison_levels[comp_dup_index]
      )
    }
  }
  return(comparison)
}

#' LIST FACTORS UI
list_factors_ui <- function(
  sample_info,
  data_file_format,
  counts_deg_method,
  id
) {
  ns <- NS(id)

  if (is.null(sample_info)) {
    return(
      HTML(
        "<font size = \"2\">A <a href=\"https://idepsite.wordpress.com/data-format/\">
        sample information file</a> can be uploaded to build a linear model according
        to experiment design. </font>"
      )
    ) 
  } else {
		factors <- colnames(sample_info)
    choices <- setNames(factors, factors)
		title <- "1. Select 1 or 2 main factors. Or leave it blank and just choose pairs
              of sample groups below."	
    if (data_file_format == 1 & counts_deg_method==3) {
      title <- "1. Select 6 or less main factors. Or skip this step and just choose 
                pairs of sample groups below."
    }
    return(
      checkboxGroupInput(
        inputId = ns("select_factors_model"), 
        h5(title), 
        choices = choices,
        selected = NULL
      )
    )	   	  
	}
}

#' LIST BLOCK FACTORS UI
list_block_factors_ui <- function(
  sample_info,
  select_factors_model,
  data_file_format,
  deg_method,
  id
) {
  ns <- NS(id)

  if (is.null(sample_info)) {
		return(NULL)		   
	} else { 
		factors <- colnames(sample_info)
		factors <- setdiff(factors, select_factors_model)

		if (length(factors) == 0 ) {
      return(NULL)
    }

		choices = setNames(factors, factors)
		title <- "Select a factor for batch effect or paired samples, if needed."	

    return(
      checkboxGroupInput(
        inputId = ns("select_block_factors_model"), 
        h5(title), 
        choices = choices,
        selected = NULL
      )
    )
	}
}

#' MODEL COMPARISONS UI
list_model_comparisons_ui <- function(
  sample_info,
  select_factors_model,
  processed_data,
  id
) {
  ns <- NS(id)

	if (is.null(sample_info) | is.null(select_factors_model)) { 
    factors <- as.character (
      detect_groups(
        colnames(processed_data)
      )
    )
		factors <- unique(factors)
		comparisons <- apply(
      t(combn(factors, 2)),
      1,
      function(x) paste(x, collapse = " vs. ")
    )
	  comparisons <- c(
      comparisons,
      apply(
        t(combn(rev(factors),2)),
        1,
        function(x) paste(x, collapse = " vs. ")
      )
    )	
		comparisons <- sort(comparisons)
		choices <- stats::setNames(gsub(" vs\\. ","-",comparisons), comparisons)

    return(
      checkboxGroupInput(
        inputId = ns("select_model_comprions"), 
			  label = h5("Select comparisons among sample groups:"),
        choices = choices,
        selected = choices[[1]]
      )
    )
  } else { 
		choices = list()
				
    for(selected_factors in select_factors_model) { 
			ix = match(selected_factors, colnames(sample_info))
				
      if(is.na(ix)) {
        next
      }

			factors <- unique(sample_info[, ix])
			comparisons <- apply(
        t(combn(factors, 2)),
        1,
        function(x) paste(x, collapse = " vs. ")
      )
				comparisons <- c(
          comparisons,
          apply(
            t(combn(rev(factors), 2)),
            1,
            function(x) paste(x, collapse = " vs. ")
          )
        )	
				comparisons <- sort(comparisons)
				comparisons <- paste0(selected_factors, ": ", comparisons)
				choices <- append(choices, stats::setNames(comparisons, comparisons))
		}
			
    if(length(choices) == 0) {
      return(NULL)
    } else {
      return(
        checkboxGroupInput(
          inputId = ns("select_model_comprions"), 
					label = h5("2. Select one or more comparisons:"), 
					choices = choices,
					selected = choices[[1]]
        )	   
      )
    }
	} 
}

#' LIST INTERACTION TERMS
list_interaction_terms_ui <- function(
  sample_info,
  select_factors_model,
  id
) {
  ns <- NS(id)
	if (is.null(sample_info) | is.null(select_factors_model)) {
    return(NULL)
	}	 else { 
		selected_factors = select_factors_model[!grepl(":", select_factors_model)]
    
    if(length(selected_factors) <= 1) {
      return(NULL) 
    }
		interactions <- apply(
      t(combn(selected_factors, 2)),
      1,
      function(x) paste(x, collapse = ":")
    )
		# choices <- setNames(interactions, interactions)
    return(
      checkboxGroupInput(
        inputId = ns("select_interactions"), 
				label = h5(
          "Interaction terms between factors(e.g. genotypes repond differently
          to treatment?):"
        ),
				choices = interactions,
        selected = NULL
      )
    )
  }
}

#' EXPERIMENT DESIGN TEXT
experiment_design_txt <- function(
  sample_info,
  select_factors_model,
  select_block_factors_model,
  select_interactions
) {
  if (is.null(sample_info) | is.null(select_factors_model)) {
    return(NULL)
	} else {
    model <- paste(
      "Model: expression ~ ",
      paste(select_factors_model, collapse = " + ")
    )
		if(!is.null(select_block_factors_model)) {
      model <- paste0(
        model,
        " + ",
        paste(select_block_factors_model, collapse = " + " )
      )
    }
    if(!is.null(select_interactions)) {
      model <- paste0(
        model,
        " + ",
        paste(select_interactions, collapse = " + ")
      )
    }
    
    return(model)									
	}
}

#' SELECT REFERENCE LEVELS
select_reference_levels_ui <- function(
  sample_info,
  select_factors_model,
  data_file_format,
  counts_deg_method,
  id
) {
  ns <- NS(id)

  if (is.null(sample_info) | is.null(select_factors_model)) {
    return(NULL)
  }	else {
    selected_factors <- select_factors_model[!grepl(":", select_factors_model)]
    
    if(length(selected_factors) == 0) {
      return(NULL)
    }
    if (!(data_file_format == 1 & counts_deg_method == 3)) {
      return(NULL)
    }
    select_choices <- c()
    for(i in selected_factors){
      if(is.na(match(i, colnames(sample_info)))) {
        select_choices[[i]] <- NULL
      } else {
        select_choices[[i]] <- unique(
          sample_info[, i]
        )      
      }
    }

    return(lapply(names(select_choices), function(x) {
        tagList(
          column(
            width = 4,
            selectInput(
              inputId = ns(
                paste0(
                  "reference_level_factor_",
                  which(names(select_choices) == x)
                )
              ), 
							label = h5(paste0("Reference/baseline level for ", x)),
							choices= setNames(
                as.list(
                  paste0(
                    x,
                    ":",
                    select_choices[[x]]
                  )
                ),
                select_choices[[x]]
              )
            )
          )
        )
      })
    )	
	}
}

#' LIMMA REACTIVE VALUE
limma_value <- function(
  data_file_format,
  counts_deg_method,
  raw_counts,
  limma_p_val,
  limma_fc,
  select_model_comprions,
  sample_info,
  select_factors_model,
  select_interactions,
  select_block_factors_model,
  factor_reference_levels,
  processed_data,
  counts_log_start,
  p_vals
) {
  if(data_file_format == 1) {
    if(counts_deg_method == 3) {
      return(
        deg_deseq2(
          raw_counts = raw_counts,
          max_p_limma = limma_p_val,
          min_fc_limma = limma_fc,
				  selected_comparisons = select_model_comprions,
          sample_info = sample_info,
					model_factors = c(select_factors_model, select_interactions), 
				  block_factor = select_block_factors_model,
          reference_levels = factor_reference_levels
        )
      )
		} else if(counts_deg_method < 3) {
      return(
        deg_limma(
          processed_data = processed_data,
          max_p_limma = limma_p_val,
          min_fc_limma = limma_fc,
					raw_counts = raw_counts,
          counts_deg_method = counts_deg_method,
					prior_counts = counts_log_start,
          data_file_format = data_file_format,
					selected_comparisons = select_model_comprions,
          sample_info = sample_info,
					model_factors = c(select_factors_model, select_interactions),
					block_factor = select_block_factors_model
        )
      )
    }
  } else if(data_file_format == 2) {
	  return(
      deg_limma(
        processed_data = processed_data,
        max_p_limma = limma_p_val,
        min_fc_limma = limma_fc,
				raw_counts = raw_counts,
        counts_deg_method = counts_deg_method,
				prior_counts = counts_log_start,
        data_file_format = data_file_format,
				selected_comparisons = select_model_comprions,
        sample_info = sample_info,
				model_factors = c(select_factors_model, select_interactions),
				block_factor = select_block_factors_model
      )
    )
	} else {
		if(!is.null(p_vals)) {
		  ix <- match(rownames(processed_data), rownames(p_vals))
		  p_vals <- p_vals[ix, ]
		}

		# Looks like ratio data, take log2
		if(sum(
      round(apply(processed_data, 2, median) + .2) == 1
    ) == dim(x)[2] & min(x) > 0) {
      processed_data <- log2(processed_data)
    }
		
		exp_type <- "None standard data without replicates."
		all_calls <- processed_data
		for(i in 1:dim(all_calls)[2]) { 
			tem <- all_calls[, i]
			all_calls[which(
        tem <= log2(limma_fc) & tem >= -log2(limma_fc)
      ) , i] <- 0			
			all_calls[which(tem > log2(limma_fc)), i] <- 1
			all_calls[which(tem < -log2(limma_fc)), i] <- -1		
			if(!is.null(p_vals)) {
        all_calls[which(p_vals[, i] > limma_p_val), i] <- 0
      } 
		}
		comparisons <- colnames(all_calls)
		extract_column <- function(
      i,
      processed_data,
      p_vals,
      top_genes
    ) {
			top_genes <- as.data.frame(processed_data[, i, drop = FALSE])
			if(is.null(p_vals)) {
        top_genes$FDR <- 0
      } else {
        top_genes$FDR <- p_vals[, i]
      }
			colnames(top_genes) <- c("Fold","FDR")
			return(top_genes)	
		} 
		top_genes <- lapply(1:dim(processed_data)[2], function(x) {
      extract_column(
        i = x,
        processed_data = processed_data,
        p_vals = p_vals,
        top_genes = top_genes
      )
    })
		top_genes <- setNames(top_genes, colnames(processed_data))
		
		return(list(
      results = all_calls,
      comparisons = colnames(processed_data),
      exp_type = exp_type,
      top_genes = top_genes
    ))
  }
}

# Differential expression using DESeq2
deg_deseq2 <- function(
  raw_counts,
  max_p_limma = .05,
  min_fc_limma = 2,
  selected_comparisons = NULL,
  sample_info = NULL,
  model_factors = NULL,
  block_factor = NULL,
  reference_levels = NULL
){
  max_samples <- 100
  max_comparisons <- 20

	groups <- as.character(detect_groups(colnames(raw_counts),sample_info))
	unique_groups <- unique(groups)	
	
	# Check for replicates, removes samples without replicates
  # Number of replicates per biological sample
	reps <- as.matrix(table(groups))
  # Less than 2 samples with replicates
	if (sum(reps[, 1] >= 2) < 2) {
    return(list(
      results= NULL,
      comparisons = NULL,
      exp_type = 
        "Failed to parse sample names to define groups. Cannot perform DEGs
         and pathway analysis. Please double check column names! Use 
         WT_Rep1, WT_Rep2 etc. ",
      top_genes = NULL
    ))
  }
	# Remove samples without replicates
	unique_groups <- rownames(reps)[which(reps[, 1] > 1)]
	ix <- which(groups %in% unique_groups)  
	groups <- groups[ix]   
	raw_counts <- raw_counts[, ix] 
	exp_type <- paste(length(unique_groups)," sample groups detected.")	
	
	# Too many samples 
	if(ncol(raw_counts)  > max_samples) { 
		return(list(
      results = NULL,
      comparisons = NULL, 
			exp_type = paste(
        exp_type,
        "Too many samples for DESeq2. Please choose limma-voom or 
         limma-trend."
      ),
			top_genes = NULL
    ))
	}	
		
	# All pair-wise comparisons
	comparisons <- ""
	for(i in 1:(length(unique_groups) - 1)) {
    for (j in (i + 1):length(unique_groups)){
      comparisons <- c(
        comparisons,
        paste(unique_groups[j], "-", unique_groups[i], sep = "")
      )
    }
  }
	comparisons <- comparisons[-1]

   # Too many comparisons 
	if(length(comparisons)  > max_comparisons) { 
		exp_type = paste(
      exp_type,
      " Too many comparisons. Only first",
      max_comparisons,
      "of the ",
      length(comparisons), 
			"comparisons calculated. Please choose comparisons."
    )
		comparisons <- comparisons[1:max_comparisons]
	}	
	
	col_data = cbind(colnames(raw_counts), groups)

	# No sample file, but user selected comparisons using column names
	if(is.null(model_factors) & length(selected_comparisons) > 0) {
		comparisons = selected_comparisons
  }

	comparison_names <- comparisons
	# Set up the DESeqDataSet Object and run the DESeq pipeline
	dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = col_data,
    design = ~groups
  )								

	if(is.null(model_factors)) {
    dds = DESeq2::DESeq(dds)
  } else {
    # Using selected factors and comparisons ----------
		# Build model
    # Block factor is just added in
		model_factors = c(model_factors, block_factor)  
    
    # Selected factors and interactions:
    # c( "strain", "treatment", "strain:treatment")
		factors <- model_factors  
    # Non-interaction terms
		factors <- factors[!grepl(":", factors)]
		# Interaction terms like strain:treatment
		interactions <- model_factors[grepl(":", model_factors)]
		
		col_data <- sample_info  
    # Factors are encoded as "A", "B", "C"; Avoids illegal letters
		factors_coded <- toupper(letters)[1: dim(col_data)[2]]
    # For look up; each column of sample_info
		names(factors_coded) <- colnames(col_data)
    # All columns named A B C D  
		colnames(col_data) <- factors_coded  

		col_data = as.data.frame(col_data)
		
		# Set reference levels for factors
    # c("genotype:wt", "treatment:control")
		if(! is.null(reference_levels) ) {
			# First factor
			for(refs in reference_levels) {
        if(!is.null(refs)) {
					# Corresponding column id for factor
          ix <- match(
            gsub(":.*", "", refs),
            colnames(sample_info)
          )
					col_data[, ix] <- as.factor(col_data[, ix])
					col_data[, ix] <- relevel(
            col_data[,ix],
            gsub(".*:", "", refs)
          )
				}
      }
		}
				
		# Base model
    deseq2_object <- paste(
      "dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts,
      colData = col_data, design = ~ ", 
			paste(factors_coded[factors], collapse = "+")
    )	
		exp_type = paste(
      "Model: ~",
      paste(model_factors, collapse = " + ")
    )

		
		# Create model
		if(length(interactions) > 0) {
			for(interaction_terms in interactions) {
				# Split strain:mutant as "strain" and "mutant"
        interacting_factors <- unlist(
          strsplit(interaction_terms, ":")
        )
        # Convert "strain:mutant" to "A:B"
				tem <- paste(
          factors_coded[interacting_factors],
          collapse = ":"
        )   
				deseq2_object <- paste(deseq2_object, " + ", tem)
			}			
		}

    # End the model
		deseq2_object <- paste(deseq2_object, ")") 

		eval(parse(text = deseq2_object) )
	
		dds = DESeq2::DESeq(dds)  # main function		


		# Comparisons 
		# "group: control vs. mutant"
		comparisons <- gsub(".*: ", "", selected_comparisons)
		comparisons <- gsub(" vs\\. ", "-", comparisons)
    # Corresponding factors for each comparison
		factors_vector <- gsub(":.*", "", selected_comparisons)
		
		# comparison_names holds names for display with real factor names
		# comparisons is used in calculation it is A, B, C for factors
		comparison_names <- comparisons

		# Note that with interaction terms, not all meaningful comparisons is
    # listed for selection. This is complex. Only under reference level.
		
		# Comparisons due to interaction terms
		if(length(interactions) > 0) {
			interaction_comparisons <- DESeq2::resultsNames(dds)
			interaction_comparisons <- interaction_comparisons[grepl(
        "\\.", interaction_comparisons
      )]
	
			comparisons <- c(comparisons, interaction_comparisons)
		
			# Translate comparisons generated in interaction terms back to factor names
			interaction_comparison_names <- interaction_comparisons
			for(i in 1:length(interaction_comparison_names)) {
				tem <- unlist(strsplit(interaction_comparison_names[i], "\\."))
				tem_factors <- substr(tem, 1, 1) 
				
        # Get the first letter and translate into real factor names
				tem_factors[1] <- names(factors_coded)[factors_coded == tem_factors[1]]
        # Get the 2nd letters and translate into real factor names  
				tem_factors[2] <- names(factors_coded)[factors_coded == tem_factors[2]]  

				interaction_comparisons[i] <- paste0(
          "I:",
          tem_factors[1],
          "_",
          substr(tem[1], 2, nchar(tem[1])),
          ".",
          tem_factors[2],
          "_",
          substr(tem[2], 2, nchar(tem[2])) 
			  )				
			}
			comparison_names = c(comparison_names,interaction_comparison_names)
		}
	}
	
	# Extract contrasts according to comparisons defined above
	result_first <- NULL
  all_calls <- NULL
	top_genes <- list()
  # Counter
  pk <- 1
  # First results
	pp <- 0 

	for(kk in 1:length(comparisons)) {
		tem = unlist(strsplit(comparisons[kk], "-"))
		
    # Group comparison using sample names
		if(is.null(model_factors)) {
      selected <- DESeq2::results(dds, contrast = c("groups", tem[1], tem[2]) )
    } else {
      # Not interaction term: they contain . interaction term
			if(!grepl("\\.", comparisons[kk])) {
        selected <- DESeq2::results(
          dds,
          contrast = c(factors_coded[factors_vector[kk]], tem[1], tem[2])
        )
      # either A, B, C ...
      } else {
        # Interaction term
        selected <- DESeq2::results(dds, name = comparisons[kk])
      }
		}

		selected$calls <- 0

		selected$calls[which(
      selected$log2FoldChange > log2(min_fc_limma) & selected$padj < max_p_limma
      )]  <-  1
		
    selected$calls[which(
      selected$log2FoldChange < -log2(min_fc_limma) & selected$padj < max_p_limma
    )] <- -1
		
    colnames(selected) <- paste(
      as.character(comparison_names[kk]),
      "___",
      colnames(selected),
      sep = ""
    )
		selected <- as.data.frame(selected)
    # First one with significant genes, collect gene list and Pval+ fold
		if(pp == 0) {
      result_first <- selected
      pp <- 1
		  top_genes[[1]] <- selected[, c(2, 6)] 
		  names(top_genes)[1] <- comparison_names[kk]
    } else {
      result_first <- merge(result_first, selected, by = "row.names") 
			rownames(result_first) <- result_first[, 1]
      result_first <- result_first[, -1]
      pk <- pk + 1
      top_genes[[pk]] <- selected[, c(2,6)]
      # Assign name to comparison
			names(top_genes)[pk] <- comparison_names[kk]
		}
	}

	interactions <- c()
	if(!is.null(model_factors)) {
    interactions <- model_factors[grepl(":", model_factors )]
  }
		
  # Add comprisons for non-reference levels. It adds to the results_first object.	
	if(length(interactions) > 0) {
    # Factor whose values are factors and names are factor and level combination
    # conditionTreated, genotypeWT
    factor_lookup <- c()
		level_lookup <- c()
		
		for(i in 1:dim(sample_info)[2]) {
      unique_sample_info <- unique(sampleInfo)
			tem <- rep(toupper(letters)[i], dim(unique_sample_info)[1])
			names(tem) <- paste0(toupper(letters)[i], unique_sample_info[, i])
			factor_lookup <- c(factor_lookup, tem)  
			
			tem <- as.character(unique_sample_info[, i])
			names(tem) <- paste0(toupper(letters)[i], unique_sample_info[, i])
			level_lookup <- c(level_lookup, tem)
		}
		
		# None interaction terms 
		none_inter_terms <- DESeq2::resultsNames(dds)[!grepl(
      "\\.", DESeq2::resultsNames(dds)
    )]
		none_inter_terms <- none_inter_terms[-1]
		all_interaction_terms <- DESeq2::resultsNames(dds)[grepl(
      "\\.", DESeq2::resultsNames(dds)
    )]

    # Each none interaction term
		for(kk in 1:length(none_inter_terms)) { 
			# Not just group comparison using sample names
      if(!is.null(model_factors)) {
				# Current factor
				c_factor <- gsub("_.*", "", none_inter_terms[kk])

				for(interaction_term in all_interaction_terms) {
          # 4 components
					splits <- split_interaction_terms(
            interaction_term,
            factor_lookup = factor_lookup,
            level_lookup = level_lookup
          )

					if (c_factor != splits[1] & c_factor != splits[3]) {
            next
          }

					selected <- DESeq2::results(
            dds,
            list(c(none_inter_terms[kk], interaction_term))
          ) 
					comparison_name <- paste0(
            none_inter_terms[kk],
            "__",
            gsub("\\.", "", interaction_term)
          )
						
					if(c_factor == splits[1]) {
            other_level <- splits[4]
          } else {
            other_level = splited[2]
          }
							
					comparison_name = paste0(
            gsub(
              "_vs_",
              "-",
              substr(none_inter_terms[kk], 3, nchar(none_inter_terms[kk]))
            ), 
						"_for_",
            other_level
          )
					comparison_names <- c(comparison_names, comparison_name)
					selected$calls <- 0   
					selected$calls[which(
            selected$log2FoldChange > log2(min_fc_limma) & selected$padj < max_p_limma
          )] <- 1
					selected$calls[which(
            selected$log2FoldChange <  -log2(min_fc_limma) & selected$padj < max_p_limma
          )] <- -1
					colnames(selected) <- paste(comparison_name, "___", colnames(selected), sep = "")
					selected <- as.data.frame(selected)
          # First one with significant genes, collect gene list and Pval+ fold
					if(pp == 0) {
						result_first <- selected
            pp = 1 
						top_genes[[1]] <- selected[, c(2, 6)] 
						names(top_genes)[1] <- comparison_name
          } else {
            result_first = merge(result_first, selected, by = "row.names") 
						rownames(result_first) <- result_first[, 1] 
						result_first <- result_first[, -1]
						pk <- pk + 1 
						top_genes[[pk]] <- selected[, c(2, 6)]
            names(top_genes)[pk] <- comparison_name
					}
				}			
			} 
		}
	}

	if(!is.null(result_first)) { 
		# Note that when you only select 1 column from a data frame it automatically 
    # converts to a vector. drop = FALSE prevents that.
		all_calls <- as.matrix(
      result_first[, grep("calls", colnames(result_first)), drop = FALSE]
    )
		colnames(all_calls) <- gsub("___.*", "", colnames(all_calls))
    # Note that samples names should have no "."
		colnames(all_calls) <- gsub("\\.", "-", colnames(all_calls))
		colnames(all_calls) <- gsub("^I-", "I:", colnames(all_calls))
	}

	return(list(
    results= all_calls,
    comparisons = comparison_names,
    exp_type = exp_type,
    top_genes = top_genes
  )) 
}

#' SPLIT INTERACION TERMS
#' Split  genotypeI.conditionTrt --> c("genotype","I","conditoin","Trt")
split_interaction_terms <- function(
  term,
  factor_lookup,
  level_lookup
) {
	if(!grepl("\\.", term)) {
    return(NULL)
  }
	terms_split <- unlist(strsplit(term, "\\."))
	# factor1, level1, factor2, level2
	return(
    c(
      factor_lookup[terms_split[1]],
      level_lookup[terms_split[1]],
      factor_lookup[terms_split[2]],
      level_lookup[terms_split[2]]
    )
  )
}

deg_limma <- function(
  processed_data,
  max_p_limma = .1,
  min_fc_limma = 2,
  raw_counts,
  counts_deg_method,
  prior_counts,
  data_file_format,
  selected_comparisons = NULL,
  sample_info = NULL,
  model_factors = NULL,
  block_factor = NULL
){
  # Many different situations:
  # 1. Just use sample names
  # 2. Just one factor
  # 3. Two factors no interaction
	# 4. Two factors with interaction
  # 5. Block factor 

	top_genes <- list()
  limma_trend <- FALSE
	if(data_file_format == 2) {
    eset <- methods::new("ExpressionSet", exprs = as.matrix(processed_data))
  } else {
    # Limma-trend method selected for counts data
		if(counts_deg_method == 1) {
      # Use transformed data for limma  
			eset <- methods::new("ExpressionSet", exprs = as.matrix(processed_data))
			limma_trend = TRUE
		}
	}

	groups <- colnames(processed_data)
	groups <-  detect_groups(groups, sample_info) 
	unique_groups <- unique(groups)  
	
	# Check for replicates, removes samples without replicates
  # Number of replicates per biological sample
	reps <- as.matrix(table(groups))
  # Less than 2 samples with replicates
	if(sum(reps[, 1] >= 2) < 2) {
    return(
      list(
        results = NULL,
        comparisons = NULL,
        exp_type = 
          "Failed to parse sample names to define groups. Cannot perform
          DEGs and pathway analysis. Please double check column names! Use
          WT_Rep1, WT_Rep2 etc. ",
        topGenes=NULL
      )
    )
  }
	 
	# Remove samples without replicates
	unique_groups <- rownames(reps)[which(reps[, 1] > 1)]
	ix <- which(groups %in% unique_groups)  
	groups <- groups[ix]   
	processed_data <- processed_data[, ix]
  raw_counts <- raw_counts[, ix] 
	
  # Just two groups
	if(length(unique_groups) == 2) {  
		unique_groups <- unique(groups)
    # "Mutant-WT"
		comparisons <-  paste(unique_groups[2], "-", unique_groups[1], sep = "")
		
		# No sample file, but user selected comparisons using column names
		if(is.null(model_factors) & length(selected_comparisons) > 0) {
      comparisons <- selected_comparisons
    }

    # Set reference level based on the order in which the levels appear
    # The first appearing level is set as reference; otherwise, we get
    # up and down-regulation reversed.
    groups <- factor(groups, levels = unique_groups) 

		design <- model.matrix(~0 + groups)
		colnames(design) <- unique_groups
		
    # Voom
		if(!is.null(raw_counts) && counts_deg_method == 2) {
			dge <- edgeR::DGEList(counts = raw_counts)
      # Normalization
			dge <- edgeR::calcNormFactors(dge, method = "TMM")
			voom_results <- limma::voom(dge, design)
      fit <- limma::lmFit(v, design)
    } else {
      # Regular limma
			fit <- limma::lmFit(eset, design)
    }

    cont_matrix <- limma::makeContrasts(contrasts = comparisons, levels = design)
		contrasts_fit <- limma::contrasts.fit(fit, cont_matrix)
		fit <- limma::eBayes(contrasts_fit, trend = limma_trend)

		# Calls differential gene expression 1 for up, -1 for down
		results <- limma::decideTests(
      fit,
      p.value = max_p_limma,
      lfc = log2(min_fc_limma)
    )

		top_genes_table <- limma::topTable(fit, number = 1e12, sort.by = "M")
		if(dim(top_genes_table)[1] != 0) {
      top_genes_table <- top_genes_table[, c('logFC', 'adj.P.Val')]
      top_genes[[1]] <- top_genes_table
    }

		# Log fold change is actually substract of means. So if the data is natral log
    # transformed, it should be natral log.
		exp_type = "2 sample groups."
  }
	
  # More than two sample groups
	if(length(unique_groups) > 2) {
	  # Set reference level based on the order in which the levels appear
    # The first appearing level is set as reference; otherwise, we get up and 
    # down-regulation reversed.
    groups <- factor(groups, levels = unique_groups)

		design <- model.matrix(~ 0 + groups)
		colnames(design) <- gsub("^groups", "", colnames(design))

		if(!is.null(raw_counts) && counts_deg_method == 2) {
			limma_voom <- limma::voom(raw_counts, design) 
			fit <- limma::lmFit(limma_voom, design) 
		} else {
      fit <- limma::lmFit(eset, design)
    }
		
		fit <- limma::eBayes(fit, trend = limma_trend)
		
		comparisons <- ""
		for(i in 1:(length(unique_groups) - 1)) {
      for (j in (i + 1):length(unique_groups)) {
        comparisons <- c(
          comparisons,
          paste(unique_groups[j], "-", unique_groups[i], sep = "")
        )
      }
    }
		comparisons <- comparisons[-1]

		# No sample file, but user selected comparisons using column names
		if(is.null(model_factors) & length(selected_comparisons) > 0) {
      comparisons <- selected_comparisons
    }
		
		make_contrast <- limma::makeContrasts(contrasts = comparisons[1], levels = design)
		if(length(comparisons) > 1) {
      for(kk in 2:length(comparisons) ) {
        make_contrast <- cbind(
          make_contrast,
          limma::makeContrasts(contrasts = comparisons[kk], levels = design)
        )
      }
    }
		exp_type = paste(length(unique_groups), " sample groups detected.")	 
		
		# Factorial design 2x2, 2x3, 3x5 etc.
		# All samples must be something like WT_control_rep1
		if(sum(sapply(strsplit(unique_groups, "_"), length) == 2) == length(unique_groups)) {
			comparisons <- ""
			for(i in 1:(length(unique_groups) - 1)) {
        for (j in (i + 1):length(unique_groups)) {
          # Only compare WT_control vs. WT_treatment
          if(strsplit(unique_groups[i], "_")[[1]][1] == strsplit(unique_groups[j], "_")[[1]][1] |
             strsplit(unique_groups[i], "_")[[1]][2] == strsplit(unique_groups[j], "_")[[1]][2]) {
            comparisons <- c(comparisons, paste(unique_groups[j], "-", unique_groups[i], sep = ""))
          }
        }
      }
			comparisons <- comparisons[-1]

			# Factors genotype treatment levels
			extract_treatment <- function(x) {
        paste(gsub(".*_", "", unlist(strsplit(x, "-"))), collapse = "-")
      }
			extract_genotype <- function (x) {
        gsub("_.*", "", unlist(strsplit(x, "-")))[1]
      }
			extract_treatment_counting <- unique(gsub(".*_", "", unlist(strsplit(unique_groups, "-"))))
			treatments <- sapply(comparisons, extract_treatment)
			genotypes <- sapply(comparisons, extract_genotype)
			exp_type <- paste(
        exp_type,
        "\nFactorial design:",
        length(unique(genotypes)),
        "X",
        length(extract_treatment_counting),
        sep = ""
      )

			# Pairwise contrasts
			make_contrast <- limma::makeContrasts(
        contrasts = comparisons[1],
        levels = design
      )
			for(kk in 2:length(comparisons)) {
        make_contrast <- cbind(
          make_contrast,
          limma::makeContrasts(contrasts = comparisons[kk], levels = design)
        )
      }
			contrast_names = colnames(make_contrast)

			# Interaction contrasts
			for (kk in 1:(length(comparisons) - 1)) {
        for(kp in (kk + 1):length(comparisons)) {
          if(treatments[kp] == treatments[kk]) {  
					  make_contrast <- cbind(
              make_contrast,
              make_contrast[, kp] - make_contrast[, kk]
            )
					  contrast_names <- c(
              contrast_names,
              paste(
                "I:", 
                genotypes[kp],
                "-",
                genotypes[kk],
                "(",
                gsub("-", ".vs.", treatments[kp]),
                ")",
                sep = ""
              )
            )
				  }
        }   
			}

			colnames(make_contrast) <- contrast_names
			comparisons <- contrast_names
		}
		

		# Sample information is uploaded and user selected factors and comparisons
		if(!is.null(model_factors) & length(selected_comparisons) > 0) {
			exp_type <- paste("Model: ~", paste(model_factors, collapse = " + "))
      # Default value to be re-write if needed
			interaction_term <- FALSE
			
			# Model factors that does not contain interaction terms
			# model_factors "genotype", "condition", "genotype:condition"
			key_model_factors <- model_factors[!grepl(":", model_factors)]
			
			# "genotype: control vs. mutant"
      # Corresponding factors for each comparison
			factors_vector <- gsub(":.*", "", selected_comparisons) 
			# Remove factors not used in comparison, these are batch effects/pairs/blocks, 	
			
			# A factor is selected both in block and main factors, then use it as block factor
			key_model_factors <- key_model_factors[
        is.na(match(key_model_factors, block_factor))
      ]	
		
      # Design matrix ----------
      # Remove factors not used.
			sample_info_filter <- sample_info[, key_model_factors, drop = F]
			groups <- apply(sample_info_filter, 1, function(x) paste(x, collapse = "_"))
			unique_groups <- unique(groups)  

      groups <- factor(groups, levels = unique_groups)

		  design <- stats::model.matrix(~ 0 + groups)  
		  colnames(design) <- gsub("^groups", "", colnames(design))
            	
			if(!is.null(raw_counts) && counts_deg_method == 2) {
				voom_results <- limma::voom(raw_counts, design) 
				fit <- limma::lmFit(voom_results, design) 
			} else {
        fit <- limma::lmFit(eset, design)
      }
		
			fit <- limma::eBayes(fit, trend = limma_trend)	
			
			# Making comaprisons-----------
      # Only one factor, or more than two then use all pairwise comparisons
			if(length(key_model_factors) != 2 | length(block_factor) > 1)  {
				comparisons <- gsub(".*: ", "", selected_comparisons)
				comparisons <- gsub(" vs\\. ", "-", comparisons)
      # Two key factors
			} else if(length(key_model_factors) == 2){ 
				if(sum(grepl(":", model_factors) > 0)) { 
					interaction_term <- TRUE
					# All pairwise comparisons
					comparisons <- ""
					for(i in 1:(length(unique_groups) - 1)) {
            for(j in (i + 1):length(unique_groups)) {
              # Only compare WT_control vs. WT_treatment
              if(strsplit(unique_groups[i], "_")[[1]][1] == strsplit(unique_groups[j], "_")[[1]][1] |
                 strsplit(unique_groups[i], "_")[[1]][2] == strsplit(unique_groups[j], "_")[[1]][2]) {
                comparisons <- c(comparisons, paste(unique_groups[j], "-", unique_groups[i], sep = ""))
              }
            }
          }
					comparisons <- comparisons[-1]
					
					# Pairwise contrasts
				  make_contrast <- limma::makeContrasts(
            contrasts = comparisons[1],
            levels = design
          )
					if(length(comparisons) > 1) {
            for(kk in 2:length(comparisons)) {
              make_contrast <- cbind(
                make_contrast,
                limma::makeContrasts(contrasts = comparisons[kk], levels = design)
              )
            }
          }
					
						
					contrast_names <- colnames(make_contrast)		
				
					# All possible interactions
					# Interaction contrasts
					
					contrast_compare <- NULL
					contrast_names <- ""
					for (kk in 1:(dim(make_contrast)[2] - 1)) {
					  for(kp in (kk+1):dim(make_contrast)[2]) {
              if(is.null(contrast_compare)) {
                contrast_compare <- make_contrast[, kp] - make_contrast[,kk]
              } else {
                contrast_compare <- cbind(
                  contrast_compare,
                  make_contrast[, kp] - make_contrast[, kk]
                )
              }
							contrast_names <- c(
                contrast_names,
                paste0(
                  "I:",
                  colnames(make_contrast)[kp],
                  ".vs.",
                  colnames(make_contrast)[kk]
                )
              )
						}   
					}
					colnames(contrast_compare) <- contrast_names[-1]
					
					# Remove nonsense contrasts from interactions
					contrast_compare <- contrast_compare[, which(
            apply(abs(contrast_compare), 2, max) == 1), drop = F
          ]
					contrast_compare <- contrast_comapre[, which(
            apply(abs(contrast_compare), 2, sum) == 4), drop = F
          ]
          # Remove duplicate columns
					contrast_compare <- t(unique(t(contrast_compare)))		
					
					# Remove unwanted contrasts involving more than three levels in 
          # either factor
					keep <- c()
					for(i in 1:dim(contrast_compare)[2]) {
						tem <- rownames(contrast_compare)[contrast_compare[, i] != 0]
						tem1 <- unique(unlist(gsub("_.*", "", tem)))
						tem2 <- unique(unlist(gsub(".*_", "", tem)))
						if(length(tem1) == 2 & length(tem2) == 2) {
              keep <- c(keep, colnames(contrast_compare)[i])
            }
					}
					contrast_comapre <- contrast_compare[, keep, drop = F]
					comparison_names = colnames(contrast_compare) 				 
				}

				# "stage: MN vs. EN"  -->  c("MN_AB-EN_AB", "EN_Nodule-EN_AB") 
				#  comparisons in all levels of the other factor 
				transform_comparisons <- function(
          comparison,
          key_model_factors
        ) {
					tem <- gsub(".*: ", "", comparison)
          # control  mutant
					tem <- unlist(strsplit(tem, " vs\\. ")) 							
					factor <- gsub(":.*", "", comparison)
          
          # 1: first factor, 2: 2nd factor
					ix <- match(factor, key_model_factors)
          # 3-1 = 2; 3-1=1
					other_factor <- key_model_factors[3 - ix]
					other_factor_levels <- unique(sample_info_filter[, other_factor])				
					comparisons <- c()
					
					for(factor_levels in other_factor_levels) {
						if(ix == 1) {
							comparisons <- c(
                comparisons,
                paste(paste0(tem, "_", factor_levels), collapse = "-")
              )
						} else {
							comparisons <- c(
                comparisons,
                paste(paste0(factor_levels, "_", tem), collapse = "-")
              )
						}
					}
					return(comparisons)		
				}	
				
				comparisons <- unlist(sapply(selected_comparisons, transform_comparisons))
				comparisons <- as.vector(comparisons)
			}
			
			# make contrasts
			make_contrast <- limma::makeContrasts(contrasts = comparisons[1], levels = design)
			if(length(comparisons) > 1) {
        for(kk in 2:length(comparisons)) {
          make_contrast <- cbind(
            make_contrast,
            limma::makeContrasts(contrasts = comparisons[kk], levels = design)
          )
        }
      }
				
			if(interaction_term) {
				make_contrast <- cbind(make_contrast, contrast_compare)
				contrast_names <- c(colnames(make_contrast), colnames(contrast_compare))
				comparisons <- c(comparisons, comparison_names)
			}
			
      # Factor is selected as block
			if(length(block_factor) >= 1) { 
				if(length(block_factor) >= 1) {
          # If multiple use the first one
          block_factor <- block_factor[1] 
        }

				block <- sample_info[, block_factor]
			
				if(!is.null(raw_counts) && counts_deg_method == 2) {
					voom_results <- limma::voom(raw_counts, design)
					corfit <- limma::duplicateCorrelation(voom_results, design, block = block)			
					fit <- limma::lmFit(
            voom_results,
            design,
            block = block,
            correlation = corfit$consensus
          ) 
				} else {
					corfit <- limma::duplicateCorrelation(eset, design, block = block)
					fit <- limma::lmFit(
            eset,
            design,
            block = block,
            correlation = corfit$consensus
          )
				}
				fit <- limma::eBayes(fit, trend = limma_trend)	
			}
		}
		
		fit_contrast <- limma::contrasts.fit(fit, make_contrast)
		fit_contrast <- limma::eBayes(fit_contrast, trend = limma_trend)
		results <- limma::decideTests(
      fit_contrast,
      p.value = max_p_limma,
      lfc = log2(min_fc_limma )
    )
		# Extract fold change for each comparison
		# There is issues with direction of foldchange. Sometimes opposite
		top <- function(
      comp,
      fit_contrast,
      processed_data
    ) {
			tem <- limma::topTable(
        fit_contrast,
        number = 1e12,
        coef = comp,
        sort.by = "M"
      ) 
			if(dim(tem)[1] == 0) {
        return(1) 
			} else { 			
				# Compute fold change for the first gene (ranked by absolute value)
				tem2 <- as.numeric(
          processed_data[which(rownames(processed_data) == rownames(tem)[1]), ]
        )
				names(tem2) <- colnames(processed_data) 
					
				return(tem[, c(1,5)]) 
			}											
		}
		
		
		top_genes <- lapply(comparisons, function(x) {
      top(
        comp = x,
        fit_contrast = fit_contrast,
        processed_data = processed_data
      )
    })
		top_genes <- setNames(top_genes, comparisons)

		ix <- which(unlist(lapply(top_genes, class)) == "numeric")
		if(length(ix) > 0) {
      top_genes <- top_genes[-ix]
    } 
	}
    
	return(list(
    results = results,
    comparisons = comparisons,
    exp_type = exp_type,
    top_genes = top_genes
  )) 
}

sig_genes_plot <- function(
  results
) {
	Up <-  apply(results, 2, function(x) sum(x == 1))
	Down <- apply(results, 2, function(x) sum(x == -1)) 
	stats <- rbind(Up, Down)
				 
	gg <- reshape2::melt(stats)

	colnames(gg) <- c("Regulation","Comparisons","Genes")
		 
	plot_bar <- ggplot2::ggplot(
    gg,
    ggplot2::aes(x = Comparisons, y = Genes, fill = Regulation)
  ) +
  ggplot2::geom_bar(position = "dodge", stat = "identity") +
  ggplot2::coord_flip() +
  ggplot2::theme(
    legend.position = "top",
    axis.title.y = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 12)
  ) +
  ggplot2::theme_light() +
  ggplot2::ylab("Number of differntially expressed genes") +
	ggplot2::geom_text(
    ggplot2::aes(label = Genes),
    position = ggplot2::position_dodge(width = 0.9),
    vjust = 0.5,
    hjust = 0
  )

	return(plot_bar)
}

#' SIG GENES STAT TABLE
genes_stat_table <- function(
  limma
) {
  results <- limma$results
  
  # If only one comparison
  if(dim(results)[2] == 1) { 
    Up <- sum(results == 1)
    Down <- sum(results == -1)
    stats <- c(colnames(results), Up, Down)
    stats <- t(as.data.frame(stats))
    row.names(stats) <- colnames(results)
    colnames(stats) <- c("Comparison","Up", "Down")
  
  # More than one comparisons
  } else {  
		Up <-  apply(results, 2, function(x) sum(x == 1))
		Down <- apply(results, 2, function(x) sum(x == -1)) 
		stats <- rbind(Up, Down)
		stats <- t(stats)
		stats <- cbind(rownames(stats), stats)
		colnames(stats)[1] <- "Comparisons"
    # Reverse row order, to be the same with plot
		stats <- stats[dim(stats)[1]:1, ] 
  }
		 
  return(as.data.frame(stats))
}

#' COMPARISON LIST FOR VENN DIAGRAM
list_comp_venn <- function(
  limma,
  up_down_regulated,
  id
) {
  ns <- NS(id)
  if(is.null(limma$comparisons)) {
    return(
      selectInput(
        inputId = ns("select_comparisons_venn"),
        label = NULL,
        choices = list("All" = "All"),
        selected = "All"
      )
    )  
	}	else {
    choices <- setNames(limma$comparisons, limma$comparisons)
		
    if(up_down_regulated) {
      tem <- c(
        paste0("Up_", limma$comparisons),
        paste0("Down_", limma$comparisons)
      )
			choices <- setNames(tem, tem)
		}
				
		choices_first_three <- choices
		
    if(length(choices_first_three) > 3) {
      # By default only 3 are selected
      choices_first_three <- choices[1:3]
    }
		return(
      checkboxGroupInput(
        inputId = ns("select_comparisons_venn"), 
			  label = h4("Select up to 5 comparisons"), 
			  choices = choices,
			  selected = choices_first_three
      )
    )	
	} 
}

#' VENN DIAGRAM PLOT FUNCTION
plot_venn <- function(
  limma,
  up_down_regulated,
  select_comparisons_venn
) {
  results <- limma$results

	# Split by up or down regulation
	if(up_down_regulated) {
    result_up <- results 
		result_up[result_up < 0] <- 0
		colnames(result_up) <- paste0("Up_", colnames(result_up))
		result_down <- results 
		result_down[result_down > 0] <- 0
		colnames(result_down) <- paste0("Down_", colnames(result_down))				
		results <- cbind(result_up, result_down)
	}			
	
  ixa <- c()
	for(comps in select_comparisons_venn) {
    # If not interaction term
		if(!grepl("^I:|^I-|^Up_I:|^Up_I-|^Down_I:|^Down_I-", comps) ) {  
			ix <- match(comps, colnames(results)) 
		} else {
			# Mismatch in comparison names for interaction terms for DESeq2
			# I:water_Wet.genetic_Hy 	in the selected Contrast
			# Diff-water_Wet-genetic_Hy  in column names
			tem <- gsub("^I-", "I:", colnames(results))
			tem <- gsub("-", "\\.", tem)
			ix <- match(comps, tem) 
      
      # This is for limma package
			if(is.na(ix)) {
        ix <- match(comps, colnames(results))
      } 						
		}
		ixa <- c(ixa,ix)
	}
  # Only use selected comparisons
	results <- results[, ixa, drop = FALSE]
	if(dim(results)[2] > 5) {
    results <- results[, 1:5]
  }
  colnames(results) <- gsub("^I-", "I:", colnames(results))	
	
  return(
    limma::vennDiagram(
      results,
      circle.col = rainbow(5),
      cex = c(1., 1, 0.7)
    )
  )
}

#' HEATMAP DATA FOR DEG
deg_heat_data <- function(
  limma,
  select_contrast,
  converted_data,
  sample_info,
  select_factors_model,
  select_model_comprions,
  factor_reference_levels,
  counts_deg_method,
  data_file_format
) {
	genes <- limma$results
	
  if(is.null(genes)) {
    return(NULL)
  }
	
  # If not interaction term
  if(!grepl("I:", select_contrast)) {
		ix <- match(input$selectContrast, colnames(genes)) 
	} else {
		# Mismatch in comparison names for interaction terms for DESeq2
		# I:water_Wet.genetic_Hy in the selected Contrast
		# Diff-water_Wet-genetic_Hy in column names
	  tem <- gsub("I-", "I:", colnames(genes))
		tem <- gsub("-", "\\.", tem)
		ix <- match(select_contrast, tem) 
		
    # This is for limma package
		if(is.na(ix)) {
      ix <- match(select_contrast, colnames(genes)) 	
    }		
			
	}
	  
	if(is.null(ix) || is.na(ix)) {
    return(NULL)
  }
  # No significant genes for this comparison
	if(sum(abs(genes[, ix])) <= 1) {
    return(NULL) 
  }
  if(dim(genes)[2] < ix) {
    return(NULL)
  }	
  query <- rownames(genes)[which(genes[, ix] != 0)]
	if(length(query) == 0) {
    return(NULL)
  }
  iy <- match(query, rownames(converted_data))
		  

	iz <- find_contrast_samples(
    select_contrast, 
		colnames(converted_data),
		sample_info,
		select_factors_model,
		select_model_comprions, 
		factor_reference_levels,
		counts_deg_method,
		data_file_format
	)
	
		# color bar
		 bar = as.vector( genes[,ix]  ); # new R versions stopped autoconvert single column data frames to vectors.
		 names(bar) = row.names( genes[,ix] )
		 bar = bar[bar!=0]

		 # retreive related data		 
		 genes = convertedData()[iy,iz,drop=FALSE]
		 
		 genes = genes[order(bar),,drop=FALSE] # needs to be sorted because myheatmap2 does not reorder genes
		 bar = sort(bar)


		 return(list(genes=genes, bar=bar ))

	
	 })
	})
}