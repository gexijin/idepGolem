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
          inputId = ns("selectModelComprions"), 
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