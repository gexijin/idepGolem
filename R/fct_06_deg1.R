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
    # choices <- setNames(factors, factors)
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
        choices = factors,
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

		# choices = setNames(factors, factors)
		title <- "Select a factor for batch effect or paired samples, if needed."	

    return(
      checkboxGroupInput(
        inputId = ns("select_block_factors_model"), 
        h5(title), 
        choices = factors,
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
        inputId = ns("select_model_comps"), 
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
  counts_deg_method
) {
  if(data_file_format == 1) {
    if(counts_deg_method == 3) {
      return(
        DEG.DESeq2( convertedCounts(),input$limmaPval, input$limmaFC,
									input$selectModelComprions, readSampleInfo(),
									c(input$selectFactorsModel,input$selectInteractions), 
									input$selectBlockFactorsModel, factorReferenceLevels() )  )
			}
    }
  }
}
  if(input$dataFileFormat == 1 ) {  # if count data
		 if(input$CountsDEGMethod == 3 ) {    # if DESeq2 method
				# rawCounts = read.csv("exampleData/airway_GSE52778.csv", row.names=1)
				# res =DEG.DESeq2(rawCounts, .05, 2) 
				# res1 =DEG.limma(rawCounts, .1, 1.5,rawCounts, 2,3) 
					
			return(   DEG.DESeq2( convertedCounts(),input$limmaPval, input$limmaFC,
									input$selectModelComprions, readSampleInfo(),
									c(input$selectFactorsModel,input$selectInteractions), 
									input$selectBlockFactorsModel, factorReferenceLevels() )  )
			}
			if(input$CountsDEGMethod < 3 )    # voom or limma-trend
				return( DEG.limma(convertedData(), input$limmaPval, input$limmaFC,
									convertedCounts(), input$CountsDEGMethod,
									priorCounts=input$countsLogStart,input$dataFileFormat,
									input$selectModelComprions, readSampleInfo(),
									c(input$selectFactorsModel,input$selectInteractions),
									input$selectBlockFactorsModel) )
	} else if (input$dataFileFormat == 2 ){ # normalized data
	 return( DEG.limma(convertedData(), input$limmaPval, input$limmaFC,
						convertedCounts(), input$CountsDEGMethod,
						priorCounts=input$countsLogStart,input$dataFileFormat,
						input$selectModelComprions, readSampleInfo(),
						c(input$selectFactorsModel,input$selectInteractions),
						input$selectBlockFactorsModel) )
	} else {   # dataFileFormat == 3 user just uploaded fold change matrix
	
		x = convertedData()
		
		pvals = convertedPvals()
		if(!is.null(pvals) ) {
		  ix = match(rownames(x), rownames(pvals))
		  pvals = pvals[ix,]
		}


		# looks like ratio data, take log2
		if( sum(round(apply(x,2, median) + .2) == 1 ) == dim(x)[2] & min(x) > 0) 
			x = log2(x)
		
		Exp.type = "None standard data without replicates."
		all.Calls = x # fake calls
		for( i in 1: dim(all.Calls)[2]) { 
			tem <- all.Calls[,i]
			all.Calls[which( tem <= log2(input$limmaFC) & tem >=  -log2(input$limmaFC) ) ,i] = 0			
			all.Calls[which( tem >  log2(input$limmaFC) ) ,i] = 1
			all.Calls[which( tem < -log2(input$limmaFC) ) ,i] = -1		
			if(!is.null(pvals) ) 
				all.Calls[ which( pvals[,i] > input$limmaPval),i] = 0
		}
		comparisons = colnames(all.Calls)
		extractColumn <- function (i) {
			topGenes = as.data.frame( convertedData()[,i,drop=FALSE])
			if(is.null(pvals) ) topGenes$FDR = 0 else 
				topGenes$FDR = pvals[,i]# fake fdr
				
			colnames(topGenes) = c("Fold","FDR")
			return(topGenes)	
		} 
		topGenes = lapply( 1:dim( x )[2], extractColumn )
		topGenes <- setNames(topGenes, colnames(x ) )
		
		return( list(results= all.Calls, comparisons = colnames(x ), Exp.type=Exp.type, topGenes=topGenes) )
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
	if(ncol(rawCounts)  > max_samples) { 
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
        paste(g[j], "-", g[i], sep = "")
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
      "dds <- DESeqDataSetFromMatrix(countData = raw_counts,
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
	
		dds = DESeq(dds)  # main function		


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