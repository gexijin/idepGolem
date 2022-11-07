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


#' Change comparison names
#'
#' Change comparison names in limma from "KO_ko-WT_ko" to    "KO-WT_for_ko"
#'
#' @param comparison Comparison to change the name for
#'
#' @export
#' @return The altered comparison string, see the description for example.
#' 
#' @family DEG functions
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

#' Create choices and title for model factors
#'
#' Create a list of options for comparisons from the
#' sample_info. Also create a title to be used in a
#' UI checkbox based off the date file format and the
#' DEG method. The list of options will become
#' checkboxes for the user to create their experiment
#' design from.
#'
#' @param sample_info Matrix of experiment design information for grouping
#' @param data_file_format Integer indicating the data format. This should be
#'   one of 1 for read counts data, 2 for normalized expression, or 3 for
#'   fold changes and adjusted P-values
#' @param counts_deg_method Integer indicating method of DEG analysis being 
#'   performed. This should be one of 1 for limma-trend, 2 for limma-voom, and 
#'   3 for DESeq2
#'  
#' @export
#' @return A list containing a string title and a vector of
#'  comparisons to choose from.
#'  
#' @family DEG functions 
list_factors_ui <- function(sample_info,
                            data_file_format = c(1, 2, 3),
                            counts_deg_method = c(1, 2, 3)) {
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
    title <- "1. Select main factors (3+ factors are not well tested). Or leave it blank and just choose pairs
              of sample groups below."
    if (data_file_format == 1 & counts_deg_method == 3) {
      title <- "1. Select 6 or less main factors. Or skip this step and just choose
                pairs of sample groups below."
    }
    return(list(
      title = title,
      choices = choices
    ))
  }
}

#' Get the block factor choices
#'
#' This function uses the sample_info file and the selected
#' factors for the model to create a selection for the batch
#' effect. Returns a vector that turns into a checkbox for
#' the User.
#'
#' @param sample_info Matrix of experiment design information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#'
#' @export
#' @return This function returns a vector of choices for a batch
#'  effect or paired samples.
#'  
#' @family DEG functions 
list_block_factors_ui <- function(sample_info,
                                  select_factors_model) {
  if (is.null(sample_info)) {
    return(NULL)
  } else {
    factors <- colnames(sample_info)
    factors <- setdiff(factors, select_factors_model)

    if (length(factors) == 0) {
      return(NULL)
    }

    choices <- setNames(factors, factors)

    return(choices)
  }
}

#' Create model comparisons choices
#'
#' This function uses the sample_info file and the selected
#' factors to create a list of options for model comparisons.
#' Changes with the input of select_factors_model. If there
#' is no selected factor then it defaults to comparisons that
#' can be created from the processed data.
#'
#' @param sample_info Matrix of experiment design information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' @param processed_data Matrix of gene data that has been through the 
#'  \code{\link{pre_process}()}
#'
#' @export
#' @return Returns a list containing a vector of choices and a
#'  title for the UI element.
#'
#' @family DEG functions 
list_model_comparisons_ui <- function(sample_info,
                                      select_factors_model,
                                      processed_data) {
  if (is.null(sample_info) | is.null(select_factors_model)) {
    factors <- as.character(
      detect_groups(
        colnames(processed_data)
      )
    )

    # group is not appropriate,
    # all in one group or too many groups
    factors <- unique(factors)
    if (length(factors) == 1 |
      length(factors) >= nrow(processed_data)
    ) {
      return(list(
        choices = NULL,
        title = "No groups!"
      ))
    }

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
    choices <- stats::setNames(gsub(" vs\\. ", "-", comparisons), comparisons)
    title <- "Select comparisons among sample groups:"

    return(list(
      choices = choices,
      title = title
    ))
  } else {
    choices <- list()

    for (selected_factors in select_factors_model) {
      ix <- match(selected_factors, colnames(sample_info))

      if (is.na(ix)) {
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
      title <- "2. Select one or more comparisons:"
    }

    if (length(choices) == 0) {
      return(NULL)
    } else {
      return(list(
        choices = choices,
        title = title
      ))
    }
  }
}

#' List model interaction terms
#'
#' This functions uses the sample info file and the selected
#' model factors to create interaction terms to be used in the
#' DEG process.
#'
#' @param sample_info Matrix of experiment design information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#'
#' @export
#' @return Returns a character string of an interaction term
#'  between the selected factors. Used in a checkbox for the
#'  User to create a model expression
#'  
#' @family DEG functions
list_interaction_terms_ui <- function(sample_info,
                                      select_factors_model) {
  if (is.null(sample_info) | is.null(select_factors_model)) {
    return(NULL)
  } else {
    selected_factors <- select_factors_model[!grepl(":", select_factors_model)]

    if (length(selected_factors) <= 1) {
      return(NULL)
    }
    interactions <- apply(
      t(combn(selected_factors, 2)),
      1,
      function(x) paste(x, collapse = ":")
    )
    return(interactions)
  }
}

#' Create a string of the model design
#'
#' Use the model design selections to create a string of the
#' model design being used for the DEG analysis.
#'
#' @param sample_info Matrix of experiment design information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' @param select_block_factors_model The selected factors for
#'  batch effect
#' @param select_interactions The interaction terms being used in
#'  the model design
#'
#' @export
#' @return Returns a string of the model design being used for
#'  the DEG analysis
#'  
#' @family DEG functions 
experiment_design_txt <- function(sample_info,
                                  select_factors_model,
                                  select_block_factors_model,
                                  select_interactions) {
  if (is.null(sample_info) | is.null(select_factors_model)) {
    return(NULL)
  } else {
    model <- paste(
      "Model: expression ~ ",
      paste(select_factors_model, collapse = " + ")
    )
    if (!is.null(select_block_factors_model)) {
      model <- paste0(
        model,
        " + ",
        paste(select_block_factors_model, collapse = " + ")
      )
    }
    if (!is.null(select_interactions)) {
      model <- paste0(
        model,
        " + ",
        paste(select_interactions, collapse = " + ")
      )
    }

    return(model)
  }
}

#' Refernce levels for selected factor
#'
#' This function uses a vector of selected factors to create
#' choices for the reference level to use for the factor in
#' the DEG analysis.
#'
#' @param sample_info Matrix of experiment design information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' @param data_file_format Integer indicating the data format. This should be
#'   one of 1 for read counts data, 2 for normalized expression, or 3 for
#'   fold changes and adjusted P-values
#' @param counts_deg_method Integer indicating method of DEG analysis being 
#'   performed. This should be one of 1 for limma-trend, 2 for limma-voom, and 
#'   3 for DESeq2
#'
#' @export
#' @return A list the same length as the vector of selected factors.
#'  Each entry in the list corresponds to the choice of group to
#'  use for the reference level.
select_reference_levels_ui <- function(sample_info,
                                       select_factors_model,
                                       data_file_format = c(1, 2, 3),
                                       counts_deg_method = c(1, 2, 3)) {
  if (is.null(sample_info) | is.null(select_factors_model)) {
    return(NULL)
  } else {
    selected_factors <- select_factors_model[!grepl(":", select_factors_model)]

    if (length(selected_factors) == 0) {
      return(NULL)
    }
    if (!(data_file_format == 1 & counts_deg_method == 3)) {
      return(NULL)
    }
    select_choices <- c()
    for (i in selected_factors) {
      if (is.na(match(i, colnames(sample_info)))) {
        select_choices[[i]] <- NULL
      } else {
        select_choices[[i]] <- unique(
          sample_info[, i]
        )
      }
    }
    return(select_choices)
  }
}

#' DEG analysis function
#'
#' Use the limma or DESeq2 package to perform DEG analysis
#' with the specified model design. Core function for the
#' DEG panel of iDEP.
#'
#' @param data_file_format Integer indicating the data format. This should be
#'   one of 1 for read counts data, 2 for normalized expression, or 3 for
#'   fold changes and adjusted P-values
#' @param counts_deg_method Integer indicating method of DEG analysis being 
#'   performed. This should be one of 1 for limma-trend, 2 for limma-voom, and 
#'   3 for DESeq2
#' @param raw_counts Matrix of  raw counts before processing for
#'  gene expression data
#' @param limma_p_val Significant p-value to use for expressed
#'  genes
#' @param limma_fc Minimum fold-change cutoff for the DEG
#'  analysis
#' @param select_model_comprions Selected comparisons to analyze
#'  in the DEG analysis
#' @param sample_info Experiment file information for grouping
#' @param select_factors_model The selected factors for the model
#'  expression
#' @param select_interactions The interaction terms being used in
#'  the model design
#' @param select_block_factors_model The selected factors for
#'  batch effect
#' @param factor_reference_levels Vector of reference levels to
#'  use for the selected factors
#' @param processed_data Matrix of gene data that has been through the 
#'  \code{\link{pre_process}()}
#' @param counts_log_start The constant added to the log transformation
#'  from pre-processing
#' @param p_vals The vector of p-vals calculated in pre-process for
#'  significant expression
#' @param threshold_wald_test TRUE/FALSE to use threshold-based Wald test
#' to test null hypothesis that the absolute value of fold-change is
#' bigger than a value. Default is FALSE 
#' @param independent_filtering TRUE/FALSE to conduct independent
#' filtering in DESeq2 results function. Default is true.
#'
#' @export
#' @return List with the results of the DEG analysis. When the function
#'  is successful there are four entries in the list. "results" is a
#'  matrix with the same dimensions as the processed data. The entries
#'  in "results" are c(-1, 0, 1) for (negative fold change, no significant
#'  change, positive fold change) respectively. The second entry is
#'  "comparisons" and is a character vector of the different comparisons
#'  that were analyzed in the function. Third is "exp_type" and details
#'  the model expression that was used for the DEG analysis. Lastly is
#'  "top_genes" which is itself a list. The "top_genes" list has an entry
#'  for each comparison. Each entry is a data frame with two columns. One
#'  column is the calculated fold change for the comparison and the other
#'  is the adjusted p-value for the fold change calculation.
limma_value <- function(data_file_format = c(1, 2, 3),
                        counts_deg_method = c(1, 2, 3),
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
                        p_vals,
                        threshold_wald_test = FALSE,
                        independent_filtering = TRUE) {
  # read counts data -----------------------------------------------------------
  if (data_file_format == 1) {

    # DESeq2----------------------------
    if (counts_deg_method == 3) {
      return(
        deg_deseq2(
          raw_counts = raw_counts,
          max_p_limma = limma_p_val,
          min_fc_limma = limma_fc,
          selected_comparisons = select_model_comprions,
          sample_info = sample_info,
          model_factors = c(select_factors_model, select_interactions),
          block_factor = select_block_factors_model,
          reference_levels = factor_reference_levels,
          threshold_wald_test = threshold_wald_test,
          independent_filtering = independent_filtering
        )
      )
      # limma-voom 2 or limma-trend 1 --------------------
    } else if (counts_deg_method < 3) {
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

    # normalized data -----------------------------------------------------------
  } else if (data_file_format == 2) {
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

    # P value, LFC data ---------------------------------------------------------
  } else {
    if (!is.null(p_vals)) {
      ix <- match(rownames(processed_data), rownames(p_vals))
      p_vals <- p_vals[ix, ]
    }

    # Looks like ratio data, take log2
    if (
      sum(
        round(
          apply(processed_data, 2, median) + .2
        ) == 1
      ) == dim(processed_data)[2] & min(processed_data) > 0
    ) {
      processed_data <- log2(processed_data)
    }

    exp_type <- "None standard data without replicates."
    all_calls <- processed_data
    for (i in 1:dim(all_calls)[2]) {
      tem <- all_calls[, i]
      all_calls[which(
        tem <= log2(limma_fc) & tem >= -log2(limma_fc)
      ), i] <- 0
      all_calls[which(tem > log2(limma_fc)), i] <- 1
      all_calls[which(tem < -log2(limma_fc)), i] <- -1
      if (!is.null(p_vals)) {
        all_calls[which(p_vals[, i] > limma_p_val), i] <- 0
      }
    }
    comparisons <- colnames(all_calls)
    extract_column <- function(i,
                               processed_data,
                               p_vals,
                               top_genes) {
      top_genes <- as.data.frame(processed_data[, i, drop = FALSE])
      if (is.null(p_vals)) {
        top_genes$FDR <- 0
      } else {
        top_genes$FDR <- p_vals[, i]
      }
      colnames(top_genes) <- c("Fold", "FDR")
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
      expr = NULL,
      top_genes = top_genes
    ))
  }
}

#' Differential expression using DESeq2 package
#'
#' Used in the limma_value function to perform DEG analysis using the
#' DESeq2 package. It is not recommended to use this function on its own.
#'
#' @param raw_counts Matrix of raw counts before processing for
#'  gene expression data
#' @param max_p_limma Significant p-value to use for the fold-change
#'  values
#' @param min_fc_limma Minimum fold-change to include in the results
#' @param selected_comparisons Comparisons being analyzed in the DEG
#'  analysis
#' @param sample_info Experiment file information for grouping
#' @param model_factors Vector of selected factors and interaction terms
#'  from the model design
#' @param block_factor The selected factors for batch effect
#' @param reference_levels Vector of reference levels to use for the
#'  selected factors
#' @param threshold_wald_test TRUE/FALSE to use threshold-based Wald test
#'  to test null hypothesis that the absolute value of fold-change is bigger than 
#'  a value.. Default is FALSE. 
#' @param independent_filtering TRUE/FALSE to conduct independent filtering. 
#'  Default is TRUE. 
#'
#' @export
#' @return The return value is the results of the DEG analysis. These
#'  results are filtered and formatted by the limma_value function.
#'   results, a data frame with up or down regulated genes for all comparisons
#'   comparisons, a vectors holding comparison_names,
#'   exp_type, a character holding experimental design or error messages.
#'   top_genes, a list, each elements hold the lfc & FDR for a comparison
deg_deseq2 <- function(raw_counts,
                       max_p_limma = .05,
                       min_fc_limma = 2,
                       selected_comparisons = NULL,
                       sample_info = NULL,
                       model_factors = NULL,
                       block_factor = NULL,
                       reference_levels = NULL,
                       threshold_wald_test = FALSE,
                       independent_filtering = TRUE) {
  # Local parameters---------------------------------------------
  max_samples <- 500
  max_comparisons <- 500

  # Define groups------------------------------------------------
  # if factors are not selected, ignore the design matrix
  # this solve the error cased when design matrix is available but
  # factors are not selected.
  if (is.null(model_factors)) {
    sample_info <- NULL
  }
  groups <- as.character(
    detect_groups(
      colnames(raw_counts),
      sample_info
    )
  )
  unique_groups <- unique(groups)


  # sample preprocess-------------------------------------------------
  # Check for replicates, removes samples without replicates
  # Number of replicates per biological sample
  reps <- as.matrix(table(groups))
  # Less than 2 samples with replicates
  if (sum(reps[, 1] >= 2) < 2) {
    return(list(
      results = NULL,
      comparisons = NULL,
      exp_type =
        "Failed to parse sample names to define groups. Cannot perform DEGs
         and pathway analysis. Please double check column names! Use
         WT_Rep1, WT_Rep2 etc. ",
      top_genes = NULL
    ))
  }

  # Remove samples without replicates (Disabled)
  # unique_groups now only holds groups with replicates
  # unique_groups <- rownames(reps)[which(reps[, 1] > 1)]
  # ix <- which(groups %in% unique_groups)
  # groups <- groups[ix]
  # raw_counts <- raw_counts[, ix]

  exp_type <- paste(length(unique_groups), " sample groups detected.")

  # Too many samples
  if (ncol(raw_counts) > max_samples) {
    return(list(
      results = NULL,
      comparisons = NULL,
      exp_type =
        "Too many samples for DESeq2. Please choose limma-voom or
         limma-trend.",
      top_genes = NULL
    ))
  }

  # If factors are selected, but comparisons are not
  if (!is.null(model_factors) & is.null(selected_comparisons)) {
    return(list(
      results = NULL,
      comparisons = NULL,
      exp_type =
        "Please select comparisons.",
      top_genes = NULL
    ))
  }

  # list comparisons----------------------------------------------------
  # All pair-wise comparisons
  comparisons <- ""
  for (i in 1:(length(unique_groups) - 1)) {
    for (j in (i + 1):length(unique_groups)) {
      comparisons <- c(
        comparisons,
        paste(unique_groups[j], "-", unique_groups[i], sep = "")
      )
    }
  }
  comparisons <- comparisons[-1]

  # Too many comparisons
  if (length(comparisons) > max_comparisons) {
    exp_type <- paste(
      " Too many comparisons. Only the first",
      max_comparisons,
      "of the ",
      length(comparisons),
      "comparisons calculated. Please choose comparisons."
    )
    comparisons <- comparisons[1:max_comparisons]
  }

  col_data <- cbind("sample" = colnames(raw_counts), groups)


  # Build DESeq2 commands
  expr <- paste0(
    "# DESeq2 script generated by iDEP on ", date(), "\n",
    "# Please cite https://doi.org/10.1186/s12859-018-2486-6\n\n",
    "# ", version$version.string, "\n",
    "if (!require(\"BiocManager\", quietly = TRUE))\n",
    "  install.packages(\"BiocManager\")  # v. ", packageVersion("BiocManager"), "\n",
    "if (!require(\"DESeq2\", quietly = TRUE))\n",
    "  BiocManager::install(\"DESeq2\")  # v. ", packageVersion("DESeq2"), "\n",
    "library(DESeq2) # D.E.G.\n",
    "FC <- ", min_fc_limma, " # Fold-change cutoff\n",
    "FDR <- ", max_p_limma, " # FDR cutoff\n"
  )

  # if padj cutoff is bigger than the default 0.1
  # alpha for results function needs to be that.
  if (max_p_limma <= 0.1) {
    alpha <- 0.1
    expr <- paste0(expr, "alpha <- 0.1 # independent filtering, default\n")
  } else {
    alpha <- max_p_limma
    expr <- paste0(expr, "alpha <- FDR # use FDR cutoff for independent filtering\n")
  }

  expr <- paste0(
    expr,
    "\n#  Prepare data --------------------\n",
    "# Use the \"Converted counts\" button in the Pre-Process tab\n",
    "# to download the filtered counts file with gene IDs converted to Ensembl.\n",
    "raw_counts = read.csv(\"converted_counts_data.csv\")\n",
    "row.names(raw_counts) <- raw_counts$User_ID\n",
    "raw_counts <- raw_counts[, -(1:3)] # delete 3 columns of IDs\n",
    "str(raw_counts)\n\n"
  )

  # No design file, but user selected comparisons using column names
  if (is.null(model_factors) && length(selected_comparisons) > 0) {
    comparisons <- selected_comparisons
  }

  comparison_names <- comparisons


  # Run DESeq2 -----------------------------------------------------------
  # Set up the DESeqDataSet Object and run the DESeq pipeline
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = col_data,
    design = ~groups
  )

  # no factors selected
  if (is.null(model_factors)) {
    dds <- DESeq2::DESeq(dds)

    expr <- paste0(
      expr,
      print_dataframe(
        df = col_data,
        df_name = "col_data"
      )
    )

    expr <- paste0(
      expr,
      "\n# Run DESeq2--------------------\n"
    )
    expr <- paste0(
      expr,
      "dds <- DESeq2::DESeqDataSetFromMatrix(\n",
      "  countData = raw_counts,\n",
      "  colData = col_data,\n",
      "  design = ~groups\n",
      ")\n",
      "dds <- DESeq2::DESeq(dds)\n"
    )
  } else { # factors selected
    # Using selected factors and comparisons ----------
    # Build model
    # Block factor is just added in
    model_factors <- c(model_factors, block_factor)

    # Selected factors and interactions:
    # c( "strain", "treatment", "strain:treatment")
    factors <- model_factors
    # only non-interaction terms
    factors <- factors[!grepl(":", factors)]
    # Interaction terms like strain:treatment
    interactions <- model_factors[grepl(":", model_factors)]

    col_data <- sample_info
    # Factors are encoded as "A", "B", "C"; Avoids illegal letters
    factors_coded <- toupper(letters)[1:dim(col_data)[2]]
    # For look up; each column of sample_info
    names(factors_coded) <- colnames(col_data)
    # All columns named A B C D
    colnames(col_data) <- factors_coded

    col_data <- as.data.frame(col_data)

    expr <- paste0(
      expr,
      "# Factors coded: ",
      paste(
        paste(names(factors_coded), "-->", factors_coded),
        collapse = ", "
      ),
      "\n"
    )

    expr <- paste0(
      expr,
      print_dataframe(
        df = col_data,
        df_name = "col_data"
      )
    )
    expr <- paste0(
      expr,
      "row.names(col_data) <- colnames(raw_counts)\n",
      "col_data\n"
    )

    # Set reference levels for factors
    # c("genotype:wt", "treatment:control")
    if (!is.null(reference_levels)) {
      expr <- paste0(
        expr,
        "\n#Set reference level \n"
      )
      # First factor
      for (refs in reference_levels) {
        if (!is.null(refs)) {
          # Corresponding column id for factor
          ix <- match(
            gsub(":.*", "", refs),
            colnames(sample_info)
          )
          col_data[, ix] <- as.factor(col_data[, ix])
          col_data[, ix] <- relevel(
            col_data[, ix],
            gsub(".*:", "", refs)
          )

          expr <- paste0(
            expr,
            "col_data[, ", ix, "]",
            " <- as.factor(col_data[, ", ix, "])\n",
            "col_data[, ", ix, "] <- relevel(col_data[, ", ix, "], \"",
            gsub(".*:", "", refs),
            "\")\n"
          )
        }
      }
    }

    # Base model
    deseq2_object <- paste(
      "dds <- DESeq2::DESeqDataSetFromMatrix(\n",
      "  countData = raw_counts,\n",
      "  colData = col_data,\n",
      "  design = ~ ",
      paste(factors_coded[factors], collapse = " + ")
    )
    exp_type <- paste(
      "Model: ~",
      paste(model_factors, collapse = " + ")
    )

    # Create model, add interaction terms if selected
    if (length(interactions) > 0) {
      for (interaction_terms in interactions) {
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
    deseq2_object <- paste(deseq2_object, "\n)")

    eval(parse(text = deseq2_object))

    dds <- DESeq2::DESeq(dds) # main function

    expr <- paste0(
      expr,
      "\n# Run DESeq2--------------------\n",
      deseq2_object, "\n",
      "dds = DESeq2::DESeq(dds) \n"
    )
    #------------------------------------------------------------------
    # list selected comparisons
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
    if (length(interactions) > 0) {
      interaction_comparisons <- DESeq2::resultsNames(dds)
      interaction_comparisons <- interaction_comparisons[grepl(
        "\\.", interaction_comparisons
      )]

      comparisons <- c(comparisons, interaction_comparisons)

      # Translate comparisons generated in interaction terms to factor names
      interaction_comparison_names <- interaction_comparisons
      for (i in 1:length(interaction_comparison_names)) {
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
      } # for all interactions
      comparison_names <- c(comparison_names, interaction_comparisons)
    } # if interaction terms
  } # if factors selected


  #-------------------------------------------------------------------------
  # Extract contrasts according to comparisons defined above
  result_first <- NULL
  all_calls <- NULL
  top_genes <- list()
  # Counter
  pk <- 1
  # First results
  pp <- 0
  expr <- paste0(expr, "\n# Extract results--------------------\n")

  for (kk in 1:length(comparisons)) {
    tem <- unlist(strsplit(comparisons[kk], "-"))

    expr <- paste0(
      expr, "\n# Comparison ", kk,
      " of ", length(comparisons), ":  ",
      comparisons[kk], "\n"
    )

    # Group comparison using sample names
    if (is.null(model_factors)) {

      # whether testing the null hypothesis FC = 0, or |FC| > 2
      if (!threshold_wald_test) {
        selected <- DESeq2::results(dds,
          contrast = c("groups", tem[1], tem[2]),
          independentFiltering = independent_filtering,
          alpha = alpha
        )
        expr <- paste0(
          expr,
          "res <- DESeq2::results(dds, \n  contrast = c(\"groups\", \"",
          tem[1], "\", \"", tem[2], "\"),\n",
          "  independentFiltering = ", independent_filtering, ",\n",
          "  alpha = alpha\n",
          ")\n"
        )
      } else {
        selected <- DESeq2::results(dds,
          contrast = c("groups", tem[1], tem[2]),
          lfcThreshold = log2(min_fc_limma),
          altHypothesis = "greaterAbs",
          independentFiltering = independent_filtering,
          alpha = alpha
        )
        expr <- paste0(
          expr,
          "res <- DESeq2::results(dds, \n  contrast = c(\"groups\", \"",
          tem[1], "\", \"", tem[2], "\"),\n",
          "  lfcThreshold = log2(FC),\n",
          "  altHypothesis = \"greaterAbs\",\n",
          "  independentFiltering = ", independent_filtering, ",\n",
          "  alpha = alpha\n",
          ")\n"
        )
      }
    } else { # factors selected

      # Not interaction term: they contain . interaction term
      if (!grepl("\\.", comparisons[kk])) {
        if (threshold_wald_test) {
          selected <- DESeq2::results(
            dds,
            contrast = c(factors_coded[factors_vector[kk]], tem[1], tem[2]),
            lfcThreshold = log2(min_fc_limma),
            altHypothesis = "greaterAbs",
            independentFiltering = independent_filtering,
            alpha = alpha
          )
          expr <- paste0(
            expr,
            "res <- DESeq2::results(dds,\n  contrast = c(\"",
            paste(c(factors_coded[factors_vector[kk]], tem[1], tem[2]),
              collapse = "\", \""
            ), "\"),\n",
            "  lfcThreshold = log2(FC),\n",
            "  altHypothesis = \"greaterAbs\",\n",
            "  independentFiltering = ", independent_filtering, ",\n",
            "  alpha = alpha\n",
            ")\n"
          )
        } else { # do threshold-based wald test
          selected <- DESeq2::results(
            dds,
            contrast = c(factors_coded[factors_vector[kk]], tem[1], tem[2]),
            independentFiltering = independent_filtering,
            alpha = alpha
          )
          expr <- paste0(
            expr,
            "res <- DESeq2::results(dds,\n  contrast = c(\"",
            paste(c(factors_coded[factors_vector[kk]], tem[1], tem[2]),
              collapse = "\", \""
            ),
            "\"),\n",
            "  independentFiltering = ", independent_filtering, ",\n",
            "  alpha = alpha\n",
            ")\n"
          )
        }
      } else { # Interaction term--------
        expr <- paste0(expr, "# Interaction term\n")
        if (threshold_wald_test) { # Wald test
          selected <- DESeq2::results(dds,
            name = comparisons[kk],
            lfcThreshold = log2(min_fc_limma),
            altHypothesis = "greaterAbs",
            independentFiltering = independent_filtering,
            alpha = alpha
          )
          expr <- paste0(
            expr,
            "res <- DESeq2::results(dds, \n  name = \"",
            comparisons[kk],
            "\",\n",
            "  lfcThreshold = log2(FC),\n",
            "  altHypothesis = \"greaterAbs\",\n",
            "  independentFiltering = ", independent_filtering, ",\n",
            "  alpha = alpha\n",
            ")\n"
          )
        } else {
          selected <- DESeq2::results(dds,
            name = comparisons[kk],
            independentFiltering = independent_filtering,
            alpha = alpha
          )
          expr <- paste0(
            expr,
            "res <- DESeq2::results(dds, name = \"",
            comparisons[kk],
            "\",\n",
            "  independentFiltering = ", independent_filtering, ",\n",
            "  alpha = alpha\n",
            ")\n"
          )
        }
      }
    }

    selected$calls <- 0

    # upregulated genes marked as 1
    selected$calls[which(
      selected$log2FoldChange > log2(min_fc_limma) &
        selected$padj < max_p_limma
    )] <- 1

    # downregulated genes marked as -1
    selected$calls[which(
      selected$log2FoldChange < -log2(min_fc_limma) &
        selected$padj < max_p_limma
    )] <- -1

    colnames(selected) <- paste(
      as.character(comparison_names[kk]),
      "___",
      colnames(selected),
      sep = ""
    )
    selected <- as.data.frame(selected)

    expr <- paste0(
      expr,
      "# Examine results \n",
      "summary(res)\n",
      "plotMA(res)\n",
      "plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[1])\n",
      "res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select\n",
      "table(sign(res$log2FoldChange)) # N. of genes Down, Up\n",
      "res <- res[order(-res$log2FoldChange), ] #sort\n",
      "head(res) #top upregulated\n",
      "tail(res) #top downregulated\n"
    )

    # First one with significant genes, collect gene list and Pval+ fold
    # top_genes is a list, whose elements are data frames lfc, FDR
    # result_first is a large data frame, collecting all results
    if (pp == 0) {
      result_first <- selected
      pp <- 1
      top_genes[[1]] <- selected[, c(2, 6)] # fold and FDR in columns 2 & 6
      names(top_genes)[1] <- comparison_names[kk]
    } else {
      result_first <- merge(result_first, selected, by = "row.names")
      rownames(result_first) <- result_first[, 1]
      result_first <- result_first[, -1]
      pk <- pk + 1
      top_genes[[pk]] <- selected[, c(2, 6)]
      # Assign name to comparison
      names(top_genes)[pk] <- comparison_names[kk]
    }
  }

  #-----------------------------------------------------------------------------
  # add comparisons for non-reference levels, deleted 7/30/2022
  #-----------------------------------------------------------------------------
  # collect up or down calls for all comparisons
  if (!is.null(result_first)) {
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
    results = all_calls,
    comparisons = comparison_names,
    exp_type = exp_type,
    expr = expr,
    top_genes = top_genes
  ))
}



#' Differential expression using limma package
#'
#' Used in the limma_value function to perform DEG analysis using the
#' limma package. It is not recommended to use this function on its own.
#'
#' @param processed_data Data that has been through the pre-processing
#' @param max_p_limma Significant p-value to use for the fold-change
#'  values
#' @param min_fc_limma Minimum fold-change to include in the results
#' @param raw_counts The matrix of counts before processing for
#'  gene expression data
#' @param counts_deg_method The method or package being used for
#'  the DEG analysis
#' @param prior_counts The constant added to the log transformation
#'  from pre-processing
#' @param data_file_format Type of gene data being examined
#' @param selected_comparisons Selected comparisons to analyze
#'  in the DEG analysis
#' @param sample_info Experiment file information for grouping
#' @param model_factors Vector of selected factors and interaction terms
#'  from the model design
#' @param block_factor The selected factors for batch effect
#'
#' @export
#' @return The return value is the results of the DEG analysis. These
#'  results are filtered and formatted by the limma_value function.
#'         processed_data = processed_data,
deg_limma <- function(processed_data,
                      max_p_limma = .1,
                      min_fc_limma = 2,
                      raw_counts,
                      counts_deg_method,
                      prior_counts,
                      data_file_format,
                      selected_comparisons = NULL,
                      sample_info = NULL,
                      model_factors = NULL,
                      block_factor = NULL) {
  # Many different situations:
  # 1. Just use sample names
  # 2. Just one factor
  # 3. Two factors 2x2 with or without interaction
  # 4. Two factors 3x2 with or without interaction
  # 5. Block factor
  # 6. Three factors no interaction
  # 7. Three factors with interaction
  # 8. Four or more factors. (not tested!!!!)

  top_genes <- list()
  limma_trend <- FALSE

  # Build DESeq2 commands
  expr <- paste0(
    "# R script for differential expression using the limma package\n",
    "# generated by iDEP on ", date(), "\n",
    "# Please cite https://doi.org/10.1186/s12859-018-2486-6\n\n",
    "# ", version$version.string, "\n",
    "if (!require(\"BiocManager\", quietly = TRUE))\n",
    "  install.packages(\"BiocManager\")  # v. ", packageVersion("BiocManager"), "\n",
    "if (!require(\"limma\", quietly = TRUE))\n",
    "  BiocManager::install(\"limma\")  # v. ", packageVersion("limma"), "\n",
    "if (!require(\"edgeR\", quietly = TRUE))\n",
    "  BiocManager::install(\"edgeR\")  # v. ", packageVersion("edgeR"), "\n",
    "if (!require(\"Biobase\", quietly = TRUE))\n",
    "  BiocManager::install(\"Biobase\")  # v. ", packageVersion("Biobase"), "\n",
    "library(limma) # D.E.G.\n",
    "library(edgeR) # Data transformation\n",
    "library(Biobase) # ExpressionSet object\n",
    "FC <- ", min_fc_limma, " # Fold-change cutoff\n",
    "FDR <- ", max_p_limma, " # FDR cutoff\n",
    "limma_trend <- FALSE \n",
    "\n# Prepare data--------------------------------\n"
  )

  # Data prep -----------------------------------------------------------------
  # if normalized data, construct ExpressionSet directly
  if (data_file_format == 2 | counts_deg_method == 1) {
    eset <- methods::new("ExpressionSet", exprs = as.matrix(processed_data))
    expr <- paste0(
      expr,
      "# Use normalized data downloaded by clicking the \"Processed Data\"\n",
      "# on the Pre-Process tab.\n",
      "df = read.csv(\"processed_data.csv\")\n",
      "row.names(df) <- df$User_ID\n",
      "df <- df[, -(1:3)] # delete 3 columns of IDs\n",
      "str(df)\n",
      "eset <- methods::new(\"ExpressionSet\", exprs = as.matrix(df))\n"
    )

    # "DESeq2" = 3,
    # "limma-voom" = 2,
    # "limma-trend" = 1
    # Use transformed data for limma-trend
    if (counts_deg_method == 1) {
      limma_trend <- TRUE

      expr <- paste0(
        expr,
        "# Using limma-trend\n",
        "limma_trend <- TRUE\n"
      )
    }
    # read counts data, use limma-voom
  } else if (!is.null(raw_counts) &&
    counts_deg_method == 2) {
    eset <- edgeR::DGEList(counts = raw_counts)
    # Normalization
    eset <- edgeR::calcNormFactors(eset, method = "TMM")
    expr <- paste0(
      expr,
      "# Use the \"Converted counts\" button in the Pre-Process tab\n",
      "# to download the filtered counts file with gene IDs converted to Ensembl.\n",
      "raw_counts = read.csv(\"converted_counts_data.csv\")\n",
      "row.names(raw_counts) <- raw_counts$User_ID\n",
      "raw_counts <- raw_counts[, -(1:3)] # delete 3 columns of IDs\n",
      "str(raw_counts)\n",
      "# Use Limma-voom\n",
      "eset <- edgeR::DGEList(counts = raw_counts)\n",
      "eset <- edgeR::calcNormFactors(eset, method = \"TMM\")\n\n"
    )
  } else { # exceptions
    return(list(
      results = NULL,
      comparisons = NULL,
      exp_type = NULL,
      expr = NULL,
      top_genes = NULL
    ))
  }


  # sample groups -------------------------------------------------------------
  sample_info_effective <- NULL
  selected_factors <- c(model_factors, block_factor)

  # if factors are selected, only use the selected factors to define groups
  if (!is.null(selected_factors)) {
    # remove interaction terms
    selected_factors <- selected_factors[
      !grepl(":", selected_factors)
    ]
    sample_info_effective <- sample_info[, selected_factors, drop = FALSE]
  }
  groups <- detect_groups(
    colnames(processed_data),
    sample_info_effective
  )
  unique_groups <- unique(groups)

  # Check for replicates, removes samples without replicates
  # Number of replicates per biological sample
  reps <- as.matrix(table(groups))
  # Less than 2 samples with replicates
  if (sum(reps[, 1] >= 2) < 2) {
    return(
      list(
        results = NULL,
        comparisons = NULL,
        exp_type =
          "Failed to parse sample names to define groups. Cannot perform
          DEGs and pathway analysis. Please double check column names! Use
          WT_Rep1, WT_Rep2 etc. ",
        expr = NULL,
        topGenes = NULL
      )
    )
  }

  # Remove samples without replicates; disabled 7/31/2022
  # unique_groups <- rownames(reps)[which(reps[, 1] > 1)]
  # ix <- which(groups %in% unique_groups)
  # groups <- groups[ix]
  # processed_data <- processed_data[, ix]
  # raw_counts <- raw_counts[, ix]


  # Two groups------------------------------------------------------------------
  # Just two groups, no design matrix uploaded
  if (is.null(model_factors) && length(unique_groups) == 2) {
    unique_groups <- unique(groups)

    expr <- paste0(
      expr,
      print_vector(groups, "groups"),
      "unique_groups <- unique(groups)\n"
    )
    # No sample file, but user selected comparisons using column names
    if (length(selected_comparisons) > 0) {
      comparisons <- selected_comparisons[1] # only use one if two selected
      expr <- paste0(
        expr,
        print_vector(selected_comparisons, "comparisons")
      )
    } else {
      # "Mutant-WT"
      comparisons <- paste(unique_groups[2], "-", unique_groups[1], sep = "")
      expr <- paste0(
        expr,
        "comparisons <-  \"",
        paste(unique_groups[2], "-", unique_groups[1], sep = ""),
        "\"\n"
      )
    }

    # Set reference level based on the order in which the levels appear
    # The first appearing level is set as reference; otherwise, we get
    # up and down-regulation reversed.
    groups <- factor(groups, levels = unique_groups)
    design <- model.matrix(~ 0 + groups)
    colnames(design) <- unique_groups

    expr <- paste0(
      expr,
      "\n# Build model----------------------------\n",
      "groups <- factor(groups, levels = unique_groups)\n",
      "design <- model.matrix(~0 + groups)\n",
      "colnames(design) <- unique_groups\n"
    )

    # Voom--------------------------
    if (!is.null(raw_counts) && counts_deg_method == 2) {
      eset <- limma::voom(eset, design)
      fit <- limma::lmFit(eset, design)

      expr <- paste0(
        expr,
        "eset <- limma::voom(eset, design)\n",
        "fit <- limma::lmFit(eset, design)\n"
      )

      # Regular limma ----------------
    } else {
      fit <- limma::lmFit(eset, design)
      expr <- paste0(expr, "fit <- limma::lmFit(eset, design)\n")
    }

    cont_matrix <- limma::makeContrasts(
      contrasts = comparisons,
      levels = design
    )
    contrasts_fit <- limma::contrasts.fit(fit, cont_matrix)
    fit <- limma::eBayes(contrasts_fit, trend = limma_trend)

    expr <- paste0(
      expr,
      "cont_matrix <- limma::makeContrasts(\n",
      "  contrasts = comparisons,\n",
      "  levels = design\n",
      ")\n",
      "contrasts_fit <- limma::contrasts.fit(fit, cont_matrix)\n",
      "fit <- limma::eBayes(contrasts_fit, trend = limma_trend)\n"
    )

    # Calls differential gene expression 1 for up, -1 for down
    results <- limma::decideTests(
      fit,
      p.value = max_p_limma,
      lfc = log2(min_fc_limma)
    )

    expr <- paste0(
      expr,
      "\n# Examine results ----------------------------\n",
      "results <- limma::decideTests(\n",
      "  fit,\n",
      "  p.value = FDR,\n",
      "  lfc = log2(FC)\n",
      ")\n",
      "summary(results)\n",
      "head(results)\n"
    )

    # LFC and FDR
    top_genes_table <- limma::topTable(fit, number = Inf, sort.by = "M")

    expr <- paste0(
      expr,
      "stats <- limma::topTable(fit, number = Inf, sort.by = \"logFC\")\n",
      "head(stats)\n",
      "table(sign(subset(stats, adj.P.Val < FDR & abs(logFC) > log2(FC))$logFC))\n",
      "plotMDS(eset) # MDS plot of original data\n",
      "plotMD(fit, status = results)\n",
      "volcanoplot(fit, names = row.names(eset), highlight = 10)\n"
    )

    # only keep FC and FDR columns
    if (dim(top_genes_table)[1] != 0) { # have rows
      top_genes_table <- top_genes_table[, c("logFC", "adj.P.Val")]
      top_genes[[1]] <- top_genes_table
    }

    # Log fold change is actually substract of means. So if the data is natral log
    # transformed, it should be natral log.
    exp_type <- "2 sample groups."
  } else {
    # 3+ groups -----------------------------------------------------------------
    # More than two sample groups, or sample design uploaded
    #  if(!is.null(model_factors) || length(unique_groups) > 2) {

    # no design file, or no comparisons selected
    if (is.null(model_factors) # no factors selected
      #|| length(selected_comparisons) == 0
    ) {
      # Set reference level based on the order in which the levels appear
      # The first appearing level is set as reference; otherwise, we get up and
      # down-regulation reversed.
      groups <- factor(groups, levels = unique_groups)

      design <- model.matrix(~ 0 + groups)
      colnames(design) <- gsub("^groups", "", colnames(design))

      expr <- paste0(
        expr,
        print_vector(groups, "groups"),
        "unique_groups <- unique(groups)\n",
        "\n# Build model----------------------------\n",
        "groups <- factor(groups, levels = unique_groups)\n",
        "design <- model.matrix(~0 + groups)\n",
        "colnames(design) <- gsub(\"^groups\", \"\", colnames(design))\n"
      )

      # Voom----------------------------
      if (!is.null(raw_counts) && counts_deg_method == 2) {
        eset <- limma::voom(eset, design)
        fit <- limma::lmFit(eset, design)
        expr <- paste0(
          expr,
          "eset <- limma::voom(eset, design)\n",
          "fit <- limma::lmFit(eset, design)\n"
        )
        # regular limma---------------------
      } else {
        fit <- limma::lmFit(eset, design)
        expr <- paste0(expr, "fit <- limma::lmFit(eset, design)\n")
      }

      fit <- limma::eBayes(fit, trend = limma_trend)
      expr <- paste0(
        expr,
        "fit <- limma::eBayes(contrasts_fit, trend = limma_trend)\n"
      )
      # all comparisons ------------------------------
      comparisons <- ""
      for (i in 1:(length(unique_groups) - 1)) {
        for (j in (i + 1):length(unique_groups)) {
          comparisons <- c(
            comparisons,
            paste(unique_groups[j], "-", unique_groups[i], sep = "")
          )
        }
      }
      comparisons <- comparisons[-1]

      # No design file-----------------------------------------------------------
      # No design file, but user selected comparisons using column names
      if (is.null(model_factors) && length(selected_comparisons) > 0) {
        comparisons <- selected_comparisons
      }


      expr <- paste0(expr, print_vector(comparisons, "comparisons"))

      make_contrast <- limma::makeContrasts(contrasts = comparisons[1], levels = design)

      expr <- paste0(
        expr,
        "make_contrast <- limma::makeContrasts(\n",
        "  contrasts = comparisons[1], \n",
        "  levels = design\n)\n"
      )

      if (length(comparisons) > 1) {
        for (kk in 2:length(comparisons)) {
          make_contrast <- cbind(
            make_contrast,
            limma::makeContrasts(contrasts = comparisons[kk], levels = design)
          )
        }
        expr <- paste0(
          expr,
          "for(kk in 2:length(comparisons) ) {\n",
          "  make_contrast <- cbind(\n",
          "    make_contrast,\n",
          "    limma::makeContrasts(contrasts = comparisons[kk], levels = design)\n",
          "  )\n",
          "}\n"
        )
      }
      exp_type <- paste(length(unique_groups), " sample groups detected.")

      # Design file uploaded ----------------------------------------------------
      # Sample information is uploaded and user selected factors and comparisons
    } else {
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
      # groups <- apply(sample_info_filter, 1, function(x) paste(x, collapse = "_"))
      groups <- detect_groups(colnames(processed_data), sample_info_filter)
      unique_groups <- unique(groups)

      groups <- factor(groups, levels = unique_groups)

      design <- stats::model.matrix(~ 0 + groups)
      colnames(design) <- gsub("^groups", "", colnames(design))

      expr <- paste0(
        expr,
        print_dataframe(sample_info, "sample_info"),
        print_vector(model_factors, "model_factors"),
        print_vector(key_model_factors, "key_model_factors"),
        "sample_info_filter <- sample_info[, key_model_factors, drop = F]\n",
        print_vector(groups, "groups"),
        "unique_groups <- unique(groups)\n",
        "\n# Build model----------------------------\n",
        "groups <- factor(groups, levels = unique_groups)\n",
        "design <- model.matrix(~0 + groups)\n",
        "colnames(design) <- gsub(\"^groups\", \"\", colnames(design))\n"
      )


      # Making comaprisons------------------------------------------------

      # not well tested for 3 factors and beyond
      # if interaction term exists, make contrast
      if (sum(grepl(":", model_factors) > 0)) {
        interaction_term <- TRUE

        # All pairwise comparisons
        comparisons <- ""
        for (i in 1:(length(unique_groups) - 1)) {
          for (j in (i + 1):length(unique_groups)) {
            # number of factors sharing the same level
            # wt_IR vs p53_mock --> 0
            # wt_IR vs wt_mock --> 1
            n_same_levels <- length(
              intersect(
                unlist(strsplit(unique_groups[i], "_")),
                unlist(strsplit(unique_groups[j], "_"))
              )
            )
            # only one factor has different levels
            if (n_same_levels == length(key_model_factors) - 1) {
              comparisons <- c(
                comparisons,
                paste(unique_groups[j],
                  "-",
                  unique_groups[i],
                  sep = ""
                )
              )
            }
          }
        }
        comparisons <- comparisons[-1]

        # Pairwise contrasts
        make_contrast <- limma::makeContrasts(
          contrasts = comparisons[1],
          levels = design
        )
        if (length(comparisons) > 1) {
          for (kk in 2:length(comparisons)) {
            make_contrast <- cbind(
              make_contrast,
              limma::makeContrasts(contrasts = comparisons[kk], levels = design)
            )
          }
        }

        contrast_names <- colnames(make_contrast)

        # All possible interactions------------------------------------
        # Interaction contrasts as the the difference between two ordinary contrasts

        contrast_interact <- NULL
        contrast_names <- ""
        for (kk in 1:(dim(make_contrast)[2] - 1)) {
          for (kp in (kk + 1):dim(make_contrast)[2]) {
            if (is.null(contrast_interact)) {
              contrast_interact <- make_contrast[, kp] - make_contrast[, kk]
            } else {
              contrast_interact <- cbind(
                contrast_interact,
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
        colnames(contrast_interact) <- contrast_names[-1]

        # Remove nonsense contrasts from interactions
        # only keep columns with 0, 1, or -1
        contrast_interact <- contrast_interact[, which(
          apply(abs(contrast_interact), 2, max) == 1
        ), drop = F]
        # only keep columns with the sum of absolute values to 4
        contrast_interact <- contrast_interact[, which(
          apply(abs(contrast_interact), 2, sum) == 4
        ), drop = F]
        # Remove duplicate columns
        contrast_interact <- t(unique(t(contrast_interact)))


        # Remove unwanted contrasts involving different factors
        keep <- c()
        for (i in 1:dim(contrast_interact)[2]) {
          # I:null_IR_yes-null_mock_yes.vs.wt_IR_no-wt_IR_yes
          tem <- gsub("I:", "", colnames(contrast_interact)[i])

          comparion_1 <- gsub("\\.vs\\..*", "", tem)
          group1 <- gsub("-.*", "", comparion_1)
          group1 <- unlist(strsplit(group1, "_"))
          group2 <- gsub(".*-", "", comparion_1)
          group2 <- unlist(strsplit(group2, "_"))
          factor1 <- which(!(group1 %in% group2))
          factor1 <- colnames(sample_info_filter)[factor1]
          # contrast factor in comparison 1

          comparion_2 <- gsub(".*\\.vs\\.", "", tem)
          group3 <- gsub("-.*", "", comparion_2)
          group3 <- unlist(strsplit(group3, "_"))
          group4 <- gsub(".*-", "", comparion_2)
          group4 <- unlist(strsplit(group4, "_"))
          factor2 <- which(!(group3 %in% group4))
          factor2 <- colnames(sample_info_filter)[factor2]

          # Filtering using selected interaction terms
          level_in_first <- setdiff(
            intersect(group1, group2), # unchanged in comparison 1
            intersect(group3, group4)
          )

          # find the factor (s) that stays in the same in each comparison
          # but differs between cthe two comparisons
          ix <- match(level_in_first, group1)
          control_factors <- colnames(sample_info_filter)[ix]

          # selected interacting factors
          interacting_terms <- model_factors[grepl(":", model_factors)]
          selected <- FALSE
          # if user selects multiple interaction terms
          for (interacting_term in interacting_terms) {
            interacting_factors <- unlist(strsplit(interacting_term, ":"))
            for (control_factor in control_factors) {
              if (setequal(
                unique(c(control_factor, factor1, factor2)),
                interacting_factors
              )) {
                selected <- TRUE
              }
            }
          }

          if (selected & factor1 == factor2) { # same factor
            keep <- c(keep, colnames(contrast_interact)[i])
          }
        }
        # drop
        contrast_interact <- contrast_interact[, keep, drop = F]

        # Remove unwanted contrasts involving more than three levels in
        # 2-factor, or more than 3 levels in 3-factor models.
        # Not well tested  for more complex models!!!!!!!!!
        keep <- c()
        for (i in 1:dim(contrast_interact)[2]) {
          tem <- rownames(contrast_interact)[contrast_interact[, i] != 0]
          # split all groups in terms of factor levels
          list1 <- lapply(tem, function(x) unlist(strsplit(x, "_")))
          df <- t(as.data.frame(list1))
          # this data frame look like. Each column a factor
          # "wt"  "mock"
          # "wt2" "mock"
          # "wt"  "IR"
          # "wt2" "IR"
          unique_levels <- apply(df, 2, function(x) length(unique(x)))

          if (sum(unique_levels <= length(unique_levels))
          == length(unique_levels)) {
            keep <- c(keep, colnames(contrast_interact)[i])
          }
        }

        # contrast_interact stores contrast due to interactions
        contrast_interact <- contrast_interact[, keep, drop = F]
      } # has interaction ?


      # formulate comparison---------------------------
      # "stage: MN vs. EN"  -->  c("MN_AB-EN_AB", "EN_Nodule-EN_AB")
      comparisons <- unlist(
        sapply(
          selected_comparisons,
          function(x) {
            transform_comparisons(
              comparison = x,
              key_model_factors = key_model_factors,
              sample_info = sample_info_filter
            )
          }
        )
      )
      comparisons <- as.vector(comparisons)

      # some comparisons does not make sense as the combination is not present
      # in design matrix

      validate_comparison <- rep(FALSE, length(comparisons))
      for (i in 1:length(comparisons)) {
        sample_1 <- gsub("-.*", "", comparisons[i])
        sample_2 <- gsub(".*-", "", comparisons[i])
        if (sum(grepl(sample_1, colnames(design)) > 0) &
          sum(grepl(sample_2, colnames(design)) > 0)
        ) {
          validate_comparison[i] <- TRUE
        }
      }
      comparisons <- comparisons[validate_comparison]

      expr <- paste0(expr, print_vector(comparisons, "comparisons"))

      # make contrasts -------------------------------------------------
      make_contrast <- limma::makeContrasts(contrasts = comparisons[1], levels = design)

      expr <- paste0(
        expr,
        "make_contrast <- limma::makeContrasts(\n",
        "  contrasts = comparisons[1], \n",
        "  levels = design\n)\n"
      )

      if (length(comparisons) > 1) {
        for (kk in 2:length(comparisons)) {
          make_contrast <- cbind(
            make_contrast,
            limma::makeContrasts(contrasts = comparisons[kk], levels = design)
          )
        }
        expr <- paste0(
          expr,
          "for(kk in 2:length(comparisons) ) {\n",
          "  make_contrast <- cbind(\n",
          "    make_contrast,\n",
          "    limma::makeContrasts(contrasts = comparisons[kk], levels = design)\n",
          "  )\n",
          "}\n"
        )
      }

      # add contrast due to interaction term ---------------------------
      if (interaction_term) {
        make_contrast <- cbind(make_contrast, contrast_interact)
        contrast_names <- c(colnames(make_contrast), colnames(contrast_interact))
        comparisons <- c(comparisons, colnames(contrast_interact))

        expr <- paste0(
          expr,
          print_dataframe(
            df = contrast_interact,
            df_name = "contrast_interact",
            convert_matrix = TRUE
          ),
          "make_contrast <- cbind(make_contrast, contrast_interact)\n",
          "contrast_names <- c(colnames(make_contrast), colnames(contrast_interact))\n",
          "comparisons <- c(comparisons, colnames(contrast_interact))\n"
        )
      }

      # No block factors ---------------------------------------------
      if (length(block_factor) == 0) {
        # voom
        if (!is.null(raw_counts) && counts_deg_method == 2) {
          eset <- limma::voom(eset, design)
          fit <- limma::lmFit(eset, design)
          expr <- paste0(
            expr,
            "eset <- limma::voom(eset, design)\n",
            "fit <- limma::lmFit(eset, design)\n"
          )
        } else {
          # limma-trend
          fit <- limma::lmFit(eset, design)
          expr <- paste0(expr, "fit <- limma::lmFit(eset, design)\n")
        }

        fit <- limma::eBayes(fit, trend = limma_trend)
        expr <- paste0(
          expr,
          "fit <- limma::eBayes(fit, trend = limma_trend)\n"
        )
      } else {
        # default use the first block
        block <- sample_info[, block_factor[1]]

        # if more than one block factor, paste them together to make one vector
        if (length(block_factor) > 1) {
          for (i in 2:length(block_factor)) {
            block <- paste(block, block_factor[i])
          }
        }

        expr <- paste0(
          expr,
          "\n# Add block factor(s)\n",
          print_vector(block, "block")
        )

        if (!is.null(raw_counts) && counts_deg_method == 2) {
          eset <- limma::voom(eset, design)
          corfit <- limma::duplicateCorrelation(eset, design, block = block)
          fit <- limma::lmFit(
            eset,
            design,
            block = block,
            correlation = corfit$consensus
          )
          expr <- paste0(
            expr,
            "eset <- limma::voom(eset, design)\n",
            "corfit <- limma::duplicateCorrelation(eset, design, block = block)\n",
            "fit <- limma::lmFit(\n",
            "  eset,\n",
            "  design,\n",
            "  block = block,\n",
            "  correlation = corfit$consensus\n",
            ")\n"
          )
        } else {
          corfit <- limma::duplicateCorrelation(eset, design, block = block)
          fit <- limma::lmFit(
            eset,
            design,
            block = block,
            correlation = corfit$consensus
          )
          expr <- paste0(
            expr,
            "corfit <- limma::duplicateCorrelation(eset, design, block = block)\n",
            "fit <- limma::lmFit(\n",
            "  eset,\n",
            "  design,\n",
            "  block = block,\n",
            "  correlation = corfit$consensus\n",
            ")\n"
          )
        }
        fit <- limma::eBayes(fit, trend = limma_trend)
        expr <- paste0(
          expr,
          "fit <- limma::eBayes(fit, trend = limma_trend)\n"
        )
      }
    }

    # Extract contrasts ---------------------------------------------------
    fit_contrast <- limma::contrasts.fit(fit, make_contrast)
    fit_contrast <- limma::eBayes(fit_contrast, trend = limma_trend)

    expr <- paste0(
      expr,
      "fit_contrast <- limma::contrasts.fit(fit, make_contrast)\n",
      "fit_contrast <- limma::eBayes(fit_contrast, trend = limma_trend)\n"
    )

    results <- limma::decideTests(
      fit_contrast,
      p.value = max_p_limma,
      lfc = log2(min_fc_limma)
    )
    expr <- paste0(
      expr,
      "\n# Examine results ----------------------------\n",
      "results <- limma::decideTests(\n",
      "  fit_contrast,\n",
      "  p.value = FDR,\n",
      "  lfc = log2(FC)\n",
      ")\n",
      "summary(results)\n",
      "head(results)\n",
      "plotMDS(eset) # MDS plot of original data\n"
    )

    top_genes <- lapply(comparisons, function(x) {
      extract_fcfdr(
        comp = x,
        fit_contrast = fit_contrast
      )
    })

    for (i in 1:length(comparisons)) {
      expr <- paste0(
        expr,
        "\n# Comparison: ", i, " of ", length(comparisons), ":", comparisons[i], "\n",
        "stats <- limma::topTable(fit_contrast, coef = ", i, ", number = Inf, sort.by = \"logFC\")\n",
        "table(sign(subset(stats, adj.P.Val < FDR & abs(logFC) > log2(FC))$logFC))\n",
        "plotMD(fit_contrast, coef = ", i, ", status = results)\n",
        "volcanoplot(fit_contrast, coef = ", i, ", names = row.names(eset), highlight = 10)\n"
      )
    }

    top_genes <- setNames(top_genes, comparisons)

    ix <- which(unlist(lapply(top_genes, class)) == "numeric")
    if (length(ix) > 0) {
      top_genes <- top_genes[-ix]
    }
  }

  return(list(
    results = results,
    comparisons = comparisons,
    exp_type = exp_type,
    expr = expr,
    top_genes = top_genes
  ))
}


#' Extract fold change and FDR for each comparison from limma
#'
#'
#' @param comp Comparison
#' @param fit_contrast fitted contrust from limma
#'  returned list
#'
#' @export
#' @return a data frame of LFC and FDR columns
extract_fcfdr <- function(comp, # comparison
                          fit_contrast # contrast
) {
  tem <- limma::topTable(
    fit_contrast,
    number = Inf,
    coef = comp,
    sort.by = "M"
  )
  if (dim(tem)[1] == 0) {
    return(1)
  } else {
    return(tem[, c(1, 5)])
  }
}


#' comparisons in all levels of the other factor
#' "stage: MN vs. EN"  -->  c("MN_AB-EN_AB", "EN_Nodule-EN_AB")
#'
#' @param comparison Comparison
#' @param key_model_factors model factors
#' @param sample_info a matrix of experimental design
#'  returned list
#'
#' @export
#' @return a list of comparisons
#
transform_comparisons <- function(comparison,
                                  key_model_factors,
                                  sample_info) {
  levels <- gsub(".*: ", "", comparison)
  # two levels of contrast: IR mock
  levels <- unlist(strsplit(levels, " vs\\. "))
  current_factor <- gsub(":.*", "", comparison)

  comparisons <- c()
  for (factor in key_model_factors) {
    if (factor == current_factor) {
      if (factor == key_model_factors[1]) { # if it is the first factor
        comparisons <- paste0( # IR-mock_
          comparisons,
          paste0(levels, collapse = "-"),
          "_"
        )
      } else { # if it is not the first: wt_  --> "wt_IR-wt_mock_"
        comparisons <- paste0(
          comparisons,
          levels[1],
          "-",
          comparisons,
          levels[2],
          "_"
        )
      }
    } else { #  not current factor
      other_factor_levels <- unique(sample_info[, factor])
      comparison0 <- c()
      for (other_factor_level in other_factor_levels) {
        # "wt-null_mock_"
        comparison1 <- paste0(comparisons, other_factor_level, "_")
        # "wt_mock-null_mock_"
        comparison1 <- gsub("-", paste0("_", other_factor_level, "-"), comparison1)
        # collect for levels
        comparison0 <- c(comparison0, comparison1)
      }
      # update for factor
      comparisons <- comparison0
    }
  }
  # remove the last "_"
  comparisons <- gsub("_$", "", comparisons)

  return(comparisons)
}

#' Print out a matrix's content as R code to regenerate the matrix
#' This is used for generating R code from DESeq2 and limma commands
#' Used mostly for design matrix in limma
#'
#' @param df a matrix or data frame of numbers or characters
#' @param df_name a string, specifying the name of the matrix to be print out
#' @param convert_matrix, whether to convert to a matrix
#'
#' @export
#' @return a string. R code that can reconstruct the matrix
print_dataframe <- function(df, df_name = "ma", convert_matrix = FALSE) {
  expr <- paste0(df_name, " <- data.frame(\n")
  for (c1 in 1:ncol(df)) {
    # if character, use quotation marks  \"
    q <- ""
    if (!is.numeric(df[, c1])) {
      q <- "\""
    }

    # export a column
    expr <- paste0(
      expr,
      "  \"", colnames(df)[c1], "\" = c(", q,
      paste0(
        df[, c1],
        collapse = paste0(q, ", ", q)
      ),
      q,
      ")"
    )
    if (c1 != ncol(df)) { # not last
      expr <- paste0(expr, ",\n")
    } else { # last one
      expr <- paste0(expr, "\n)\n")
    }
  }

  # convert to matrix
  if (convert_matrix) {
    expr <- paste0(
      expr,
      df_name,
      " <- as.matrix(",
      df_name,
      ")\n"
    )
  }
  return(expr)
}


#' print the content of a vector as R code for the reconstruc of the same vector
#'
#'
#' @param x vector
#' @param name name for the vector
#'  returned
#'
#' @export
#' @return a string object containing the R code
print_vector <- function(x, # comparison
                         name # contrast
) {
  # if character, use quotation marks  \"
  q <- ""
  if (!is.numeric(x)) {
    q <- "\""
  }
  expr <- paste0(
    name,
    " <- c(",
    q,
    paste0(
      x,
      collapse = paste0(q, ", ", q)
    ),
    q,
    ")\n"
  )
  return(expr)
}





#' Significant genes bar plot
#'
#' Create a bar plot of the number of genes with a significant fold
#' change. The plot will break down all of the comparisons that
#' were analyzed and the count of significant genes for both up
#' and down changes.
#'
#' @param results Results matrix from the limma_value function
#'  returned list
#'
#' @export
#' @return Formatted gg barplot of the significantly expressed
#'  genes.
sig_genes_plot <- function(results) {
  Up <- apply(results, 2, function(x) sum(x == 1))
  Down <- apply(results, 2, function(x) sum(x == -1))
  stats <- rbind(Up, Down)

  gg <- reshape2::melt(stats)

  colnames(gg) <- c("Regulation", "Comparisons", "Genes")

  plot_bar <- ggplot2::ggplot(
    gg,
    ggplot2::aes(x = Comparisons, y = Genes, fill = Regulation)
  ) +
    ggplot2::geom_bar(position = "dodge", stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "top",
      axis.title.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
    ) +
    ggplot2::ylab("Number of differntially expressed genes") +
    ggplot2::geom_text(
      ggplot2::aes(label = Genes),
      position = ggplot2::position_dodge(width = 0.9),
      vjust = 0.5,
      hjust = 0
    ) +
    ggplot2::ylim(0, max(gg$Genes) * 1.1)

  return(plot_bar)
}

#' Create a table from DEG results
#'
#' Using the limma_value return list, create a table of the
#' number of significantly expressed genes for each analyzed
#' comparison.
#'
#' @param limma Return list from the limma_value function
genes_stat_table <- function(limma) {
  results <- limma$results
  if (is.null(results)) {
    return(NULL)
  }
  # If only one comparison
  if (dim(results)[2] == 1) {
    Up <- sum(results == 1)
    Down <- sum(results == -1)
    stats <- c(colnames(results), Up, Down)
    stats <- t(as.data.frame(stats))
    row.names(stats) <- colnames(results)
    colnames(stats) <- c("Comparison", "Up", "Down")

    # More than one comparisons
  } else {
    Up <- apply(results, 2, function(x) sum(x == 1))
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

#' List comparisons for venn diagram plot
#'
#' Create a list of the comparisons that were detected and
#' analyzed by the limma_value function. These comparisons
#' can be used in the plot_venn function to find the number
#' of overlapping significantly expressed genes.
#'
#' @param limma Returned list of results from the limma_value
#'  function
#' @param up_down_regulated Split the comparisons into either
#'  up or down regulated
#'
#' @export
#' @return A character vector of the comparisons that were used
#'  in the DEG analysis and can be plotted with the venn_plot
#'  function.
list_comp_venn <- function(limma,
                           up_down_regulated) {
  if (is.null(limma$comparisons)) {
    return(NULL)
  } else {
    choices <- stats::setNames(limma$comparisons, limma$comparisons)

    if (up_down_regulated) {
      tem <- c(
        paste0("Up_", limma$comparisons),
        paste0("Down_", limma$comparisons)
      )
      choices <- stats::setNames(tem, tem)
    }

    choices_first_three <- choices

    if (length(choices_first_three) > 3) {
      # By default only 3 are selected
      choices_first_three <- choices[1:3]
    }
    return(list(
      choices = choices,
      choices_first_three = choices_first_three
    ))
  }
}

#' Prep data for venn diagram
#'
#' @description This function prepares the data from the DEG results
#' with some data transformations and filtering based on specified
#' parameters. The data is then passed into the \code{plot_venn} and
#' \code{plot_upset} functions.
#'
#' @param limma Returned list of results from the limma_value
#'  function
#' @param up_down_regulated Split the comparisons into either
#'  up or down regulated
#' @param select_comparisons_venn The comparisons to plot on the
#'  venn diagram
#'
#' @return A data frame of formatted data for use in a venn diagram or
#' upset plot
#' @export
prep_venn <- function(limma,
                      up_down_regulated,
                      select_comparisons_venn) {
  results <- limma$results

  # Split by up or down regulation
  if (up_down_regulated) {
    result_up <- results
    result_up[result_up < 0] <- 0
    colnames(result_up) <- paste0("Up_", colnames(result_up))
    result_down <- results
    result_down[result_down > 0] <- 0
    colnames(result_down) <- paste0("Down_", colnames(result_down))
    results <- cbind(result_up, result_down)
  }

  ixa <- c()
  for (comps in select_comparisons_venn) {
    # If not interaction term
    if (!grepl("^I:|^I-|^Up_I:|^Up_I-|^Down_I:|^Down_I-", comps)) {
      ix <- match(comps, colnames(results))
    } else {
      # Mismatch in comparison names for interaction terms for DESeq2
      # I:water_Wet.genetic_Hy   in the selected Contrast
      # Diff-water_Wet-genetic_Hy  in column names
      tem <- gsub("^I-", "I:", colnames(results))
      tem <- gsub("-", "\\.", tem)
      ix <- match(comps, tem)

      # This is for limma package
      if (is.na(ix)) {
        ix <- match(comps, colnames(results))
      }
    }
    ixa <- c(ixa, ix)
  }
  # Only use selected comparisons
  results <- results[, ixa, drop = FALSE]

  colnames(results) <- gsub("^I-", "I:", colnames(results))

  return(results)
}

#' Create a venn diagram plot
#'
#' Plot a venn diagram that illustrates the number of significantly
#' expressed genes that overlap for multiple comparisons.
#'
#' @param limma limma Returned list of results from the limma_value
#'  function
#' @param up_down_regulated Split the comparisons into either
#'  up or down regulated
#' @param select_comparisons_venn The comparisons to plot on the
#'  venn diagram
#'
#' @export
#' @return A formatted venn diagram plot of the selected comparisons.
plot_venn <- function(results) {
  if (dim(results)[2] > 5) {
    results <- results[, 1:5]
  }

  return(
    limma::vennDiagram(
      results,
      circle.col = rainbow(5),
      cex = c(1., 1, 0.7)
    )
  )
}

#' Plot upset graph
#'
#' @param results Dataframe from \code{venn_data()}
#'
#' @return A \code{ggplot2} object
#'
#' @export
#'
plot_upset <- function(results) {

  # get the groups of data by category
  data <- results |>
    tidyr::as_tibble() |>
    dplyr::mutate(id = 1:dplyr::n()) |>
    tidyr::gather(Group, GroupMember, dplyr::all_of(colnames(results))) |>
    dplyr::group_by_at(dplyr::vars(-c(Group, GroupMember))) |>
    tidyr::nest() |>
    dplyr::mutate(
      Group = purrr::map(
        data,
        ~ if (all(.x$GroupMember == 0)) {
          NA
        } else {
          .x$Group[.x$GroupMember != 0]
        }
      )
    ) |>
    dplyr::filter(!is.na(Group))

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = Group)) +
    ggplot2::geom_bar(fill = "grey") +
    ggplot2::geom_text(
      stat = "count",
      ggplot2::aes(label = ggplot2::after_stat(count)),
      vjust = -1,
      size = 5
    ) +
    ggplot2::theme_light() +
    ggplot2::scale_y_continuous(limits = function(x) {
      x * 1.1
    }) +
    ggupset::scale_x_upset() +
    ggupset::theme_combmatrix(
      combmatrix.panel.line.size = 0,
      combmatrix.panel.point.size = 4,
      combmatrix.label.text = ggplot2::element_text(size = 17),
      combmatrix.label.extra_spacing = 10
    ) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(NULL)

  return(plot)
}


#' Find data for the heatmap
#'
#' Filter the processed data into a submatrix that only contains
#' the genes that had a significant fold change. The samples will
#' also be subsetted to only contain samples that pertain to the
#' selected comparison or contrast.
#'
#' @param limma Return results list from the limma_value function
#' @param select_contrast Comparison from DEG analysis to filter
#'  for the significant genes
#' @param processed_data Data that has been through the pre-processing
#' @param contrast_samples Columns that are in the group of the
#'  selected comparison
#'
#' @export
#' @return Submatrix of the processed data with only the significantly
#'  expressed genes and the columns that are in the selected contrast
#'  group.
deg_heat_data <- function(limma,
                          select_contrast,
                          processed_data,
                          contrast_samples) {
  genes <- limma$results

  if (is.null(genes)) {
    return(NULL)
  }

  # If not interaction term
  if (!grepl("I:", select_contrast)) {
    ix <- match(select_contrast, colnames(genes))
  } else {
    # Mismatch in comparison names for interaction terms for DESeq2
    # I:water_Wet.genetic_Hy in the selected Contrast
    # Diff-water_Wet-genetic_Hy in column names
    tem <- gsub("I-", "I:", colnames(genes))
    tem <- gsub("-", "\\.", tem)
    ix <- match(select_contrast, tem)

    # This is for limma package
    if (is.na(ix)) {
      ix <- match(select_contrast, colnames(genes))
    }
  }

  if (is.null(ix) || is.na(ix)) {
    return(NULL)
  }
  # No significant genes for this comparison
  if (sum(abs(genes[, ix])) <= 1) {
    return(NULL)
  }
  if (dim(genes)[2] < ix) {
    return(NULL)
  }
  query <- rownames(genes)[which(genes[, ix] != 0)]
  if (length(query) == 0) {
    return(NULL)
  }
  iy <- match(query, rownames(processed_data))


  iz <- contrast_samples

  # Color bar
  bar <- as.vector(genes[, ix])
  names(bar) <- row.names(genes)
  bar <- bar[bar != 0]

  # Retreive related data
  genes <- processed_data[iy, iz, drop = FALSE]

  genes <- genes[order(bar), , drop = FALSE]
  bar <- sort(bar)

  return(list(
    genes = genes,
    bar = bar
  ))
}

#' Data processing for volcano and ma plot
#'
#' @param select_contrast Comparison from DEG analysis to filter
#'  for the significant genes
#' @param comparisons The comparisons vector from the results list
#'  of the limma_value function
#' @param top_genes top_genes list from results list of the
#'  limma_value function
#' @param limma_p_val Significant p-value to use to in determining
#'  the expressed genes
#' @param limma_fc Minimum fold change value to use in determining
#'  the expressed genes
#' @param processed_data Data matrix that has gone through
#'  pre-processing
#' @param contrast_samples Samples that are included in the selected
#'  comparison
#'
#' @return A list containing processed data for volcano plot in
#'  \code{plot_volcano()} & \code{plot_ma()} and list of differently expressed
#'  genes for labeling
#' @export
#'
volcano_data <- function(select_contrast,
                         comparisons,
                         top_genes,
                         limma_p_val,
                         limma_fc,
                         processed_data,
                         contrast_samples) {
  if (is.null(select_contrast) || is.null(comparisons) ||
    length(top_genes) == 0) {
    return(NULL)
  }
  if (length(comparisons) == 1) {
    top_1 <- top_genes[[1]]
  } else {
    top <- top_genes
    ix <- match(select_contrast, names(top))
    if (is.na(ix)) {
      return(NULL)
    }
    top_1 <- top[[ix]]
  }
  if (dim(top_1)[1] == 0) {
    grid::grid.newpage()
    return(grid::grid.text("Not available.", 0.5, 0.5))
  }
  colnames(top_1) <- c("Fold", "FDR")
  # Convert to data frame
  top_1 <- as.data.frame(top_1)
  # Remove NA's
  top_1 <- top_1[which(!(is.na(top_1$Fold) | is.na(top_1$FDR))), ]
  top_1$upOrDown <- "None"
  top_1$upOrDown[which(
    top_1$FDR <= limma_p_val & top_1$Fold >= log2(limma_fc)
  )] <- "Up"
  top_1$upOrDown[which(
    top_1$FDR <= limma_p_val & top_1$Fold <= -log2(limma_fc)
  )] <- "Down"

  iz <- contrast_samples

  average_data <- as.data.frame(apply(processed_data[, iz], 1, mean))
  colnames(average_data) <- "Average"
  rownames(average_data) <- rownames(processed_data)

  genes <- merge(average_data, top_1, by = "row.names")

  anotate_genes <- genes |>
    dplyr::filter(upOrDown != "None")

  return(list(
    data = genes,
    anotate_genes = anotate_genes
  ))
}


#' Volcano DEG plot
#'
#' Use the results from limma-value to create a volcano plot to
#' illustrate the significantly expressed genes. Up and down
#' regulated genes are colored on the ggplot.
#'
#' @param plot_colors List containing three colors to differentiate between
#'  the up-regulated, down-regulated, and other genes
#' @param anotate_genes List of gene names to anotate. Default is NULL.
#' @param data Dataframe of processed data from \code{volcano_data()}.
#'
#' @export
#' @return ggplot with the fold value as the X-axis and the log 10
#'  value of the adjusted p-value as the Y-axis.
plot_volcano <- function(data,
                         anotate_genes = NULL,
                         plot_colors) {
  anotate_data <- data |>
    dplyr::filter(Row.names %in% anotate_genes)

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = Fold, y = -log10(FDR))) +
    ggplot2::geom_point(ggplot2::aes(color = upOrDown)) +
    ggplot2::scale_color_manual(values = plot_colors) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        size = 16
      ),
      axis.text.y = ggplot2::element_text(
        size = 16
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      )
    ) +
    ggplot2::labs(
      title = "Fold Change vs. Adjusted p-Value",
      y = "-log10(Adjusted p-Val)",
      x = "Fold Change",
      color = "Regulated"
    ) +
    ggrepel::geom_text_repel(
      data = anotate_data,
      ggplot2::aes(label = Row.names),
      size = 3,
      min.segment.length = 0,
      max.time = 1,
      max.overlaps = 25,
      direction = "both",
      nudge_x = .5,
      nudge_y = 2
    )

  return(plot)
}

#' Plot mean expression and fold change
#'
#' Draw a ggplot of the overal mean expression for each gene
#' and the calculated fold change for the selected comparison.
#'
#' @param data Dataframe of processed data from \code{volcano_data()}
#' @param anotate_genes List containing the gene names to be labeled on the
#'   plot. Default is NULL.
#' @param plot_colors List containing three colors to differentiate between
#'   the up-regulated, down-regulated, and other genes
#'
#' @export
#' @return A ggplot with the X-axis the mean expression value and
#'  the Y-axis the calculated fold-change from the DEG analysis.
plot_ma <- function(data,
                    plot_colors,
                    anotate_genes = NULL) {
  anotate_data <- data |>
    dplyr::filter(Row.names %in% anotate_genes)

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = Average, y = Fold)) +
    ggplot2::geom_point(ggplot2::aes(color = upOrDown)) +
    ggplot2::scale_color_manual(values = plot_colors) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        size = 16
      ),
      axis.text.y = ggplot2::element_text(
        size = 16
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      )
    ) +
    ggplot2::labs(
      title = "Average Expression vs. Log2 Fold Change",
      y = "Log2 Fold Change",
      x = "Average Expression",
      color = "Regulated"
    ) +
    ggrepel::geom_text_repel(
      data = anotate_data,
      ggplot2::aes(label = Row.names),
      size = 3,
      min.segment.length = 0,
      max.time = 2,
      max.overlaps = 25,
      direction = "both",
      nudge_x = 0.5,
      nudge_y = 2
    )

  return(
    plot
  )
}

#' Scatter plot of the comparison groups
#'
#' Create a scatter plot of the expression value for each gene
#' in the two groups for the selected contrast. For the selected
#' contrast, the mean expression is calculated for a gene in both
#' group of samples in the contrast and plotted in a scatter plot.
#'
#' @param select_contrast Comparison from DEG analysis to filter
#'  for the significant genes
#' @param comparisons The comparisons vector from the results list
#'  of the limma_value function
#' @param top_genes top_genes list from results list of the
#'  limma_value function
#' @param limma_p_val Significant p-value to use to in determining
#'  the expressed genes
#' @param limma_fc Minimum fold change value to use in determining
#'  the expressed genes
#' @param contrast_samples Samples that are included in the selected
#'  comparison
#' @param processed_data Data matrix that has gone through
#'  pre-processing
#' @param sample_info Experiment file information for grouping
#'
#' @export
#' @return A formatted ggplot with the X-axis as the mean expression
#'  of one contrast group and the Y-axis as the mean expression of
#'  the other contrast group.
plot_deg_scatter <- function(select_contrast,
                             comparisons,
                             top_genes,
                             limma_p_val,
                             limma_fc,
                             contrast_samples,
                             processed_data,
                             sample_info) {
  if (grepl("I:", select_contrast)) {
    grid::grid.newpage()
    return(
      grid::grid.text("Not available for interaction terms.", 0.5, 0.5)
    )
  }

  if (length(comparisons) == 1) {
    top_1 <- top_genes[[1]]
  } else {
    top <- top_genes
    ix <- match(select_contrast, names(top))
    if (is.na(ix)) {
      grid::grid.newpage()
      return(grid::grid.text("Not available.", 0.5, 0.5))
    }
    top_1 <- top[[ix]]
  }

  if (dim(top_1)[1] == 0) {
    grid::grid.newpage()
    return(grid::grid.text("Not available.", 0.5, 0.5))
  }
  colnames(top_1) <- c("Fold", "FDR")
  # Convert to data frame
  top_1 <- as.data.frame(top_1)
  # Remove NA's
  top_1 <- top_1[which(!(is.na(top_1$Fold) | is.na(top_1$FDR))), ]
  top_1$upOrDown <- "None"
  top_1$upOrDown[which(
    top_1$FDR <= limma_p_val & top_1$Fold >= log2(limma_fc)
  )] <- "Up"
  top_1$upOrDown[which(
    top_1$FDR <= limma_p_val & top_1$Fold <= -log2(limma_fc)
  )] <- "Down"

  iz <- contrast_samples

  genes <- processed_data[, iz]

  g <- detect_groups(colnames(genes), sample_info)

  if (length(unique(g)) > 2) {
    plot.new()
    text(0.5, 0.5, "Not available for more than two groups.")
  } else {
    average_1 <- apply(genes[, which(g == unique(g)[1])], 1, mean)

    average_2 <- apply(genes[, which(g == unique(g)[2])], 1, mean)

    genes_1 <- cbind(average_1, average_2)
    rownames(genes_1) <- rownames(genes)
    genes_1 <- merge(genes_1, top_1, by = "row.names")

    colors <- gg_color_hue(2)

    return(
      ggplot2::ggplot(genes_1, ggplot2::aes(x = average_1, y = average_2)) +
        ggplot2::geom_point(ggplot2::aes(color = upOrDown)) +
        ggplot2::scale_color_manual(values = c(colors[1], "grey45", colors[2])) +
        ggplot2::theme_light() +
        ggplot2::theme(
          legend.position = "right",
          axis.title.x = ggplot2::element_text(
            color = "black",
            size = 14
          ),
          axis.title.y = ggplot2::element_text(
            color = "black",
            size = 14
          ),
          axis.text.x = ggplot2::element_text(
            size = 16
          ),
          axis.text.y = ggplot2::element_text(
            size = 16
          ),
          plot.title = ggplot2::element_text(
            color = "black",
            size = 16,
            face = "bold",
            hjust = .5
          )
        ) +
        ggplot2::labs(
          title = "Average Expression in Group",
          y = paste0("Average Expression: ", unique(g)[2]),
          x = paste0("Average Expression: ", unique(g)[1]),
          color = "Regulated"
        )
    )
  }
}


#' returns results from DEG analysis
#'
#'
#' @param limma_value Results from DESeq2 or limma
#' @param gene_names Gene IDs
#' @param processed_data  Processed data
#' @param no_id_conversion If true, gene ids will not be mapped to Ensembl
#'
#' @export
#' @return A list. deg_data and limma$results
deg_information <- function(limma_value,
                            gene_names,
                            processed_data,
                            no_id_conversion = FALSE) {
  if (no_id_conversion) {

    # get the first comparison level
    degs_data <- limma_value$top_genes[[1]]
    degs_data$User_ID <- rownames(degs_data)

    # get the additional comparison levels if they exists
    if (length(names(limma_value$top_genes)) > 1) {
      for (i in 2:length(names(limma_value$top_genes))) {
        temp <- limma_value$top_genes[[i]]
        temp$User_ID <- rownames(temp)
        degs_data <- dplyr::inner_join(degs_data, temp, by = "User_ID")
      }
    }

    # connect to gene symbols and original user id
    processed_data <- as.data.frame(processed_data)
    processed_data$User_ID <- rownames(processed_data)

    degs_data <- dplyr::full_join(degs_data, gene_names, by = "User_ID")
    degs_data <- dplyr::full_join(degs_data, processed_data, by = "User_ID")


    degs_data <- degs_data |>
      dplyr::relocate(User_ID)
  } else {

    # get the first comparison level
    degs_data <- limma_value$top_genes[[1]]
    colnames(degs_data) <- c(
      (paste(limma_value$comparisons[[1]], "logFC", sep = "_")),
      (paste(limma_value$comparisons[[1]], "adjPval", sep = "_"))
    )
    degs_data$ensembl_ID <- rownames(degs_data)

    # get the additional comparison levels if they exists
    if (length(names(limma_value$top_genes)) > 1) {
      for (i in 2:length(names(limma_value$top_genes))) {
        temp <- limma_value$top_genes[[i]]
        colnames(temp) <- c(
          (paste(limma_value$comparisons[[i]], "logFC", sep = "_")),
          (paste(limma_value$comparisons[[i]], "adjPval", sep = "_"))
        )
        temp$ensembl_ID <- rownames(temp)
        degs_data <- dplyr::inner_join(degs_data, temp, by = "ensembl_ID")
      }
    }

    # connect to gene symbols and original user id
    processed_data <- as.data.frame(processed_data)
    processed_data$ensembl_ID <- rownames(processed_data)
    degs_data$"Processed data:" <- "" # add a empty column
    degs_data <- dplyr::full_join(degs_data, gene_names, by = "ensembl_ID")
    degs_data <- dplyr::full_join(degs_data, processed_data, by = "ensembl_ID")


    degs_data <- degs_data |>
      dplyr::relocate(User_ID) |>
      dplyr::relocate(ensembl_ID) |>
      dplyr::relocate(symbol)
  }
  return(list(degs_data, limma_value$Results))
}

#' UI component to customize gene labels
#'
#' This component contains an action button to activate the pop-up modal to
#' customize how genes are labeled on the volcano or ma plot.
#'
#' @param id Namespace ID
#'
#' @return Shiny module
#' @export
mod_label_ui <- function(id) {
  ns <- shiny::NS(id)
  tagList(
    actionButton(
      inputId = ns("customize_labels"),
      label = "Customize gene labels"
    )
  )
}

#' Server component to customize gene labels
#'
#' This component contains the pop-up modal and data processing to customize
#'  how genes are labeled on the volcano or ma plot. Users have the option to
#'  not label genes, label specific genes, label top n genes by a certain value,
#'  or label genes above certain threshold.
#'
#' @param id Namespace ID
#' @param data_list List of gene data from \code{volcano_data()}
#' @param method String designating if the results are for the volcano plot or
#'   ma plot
#'
#' @return A shiny module.
#' @export
mod_label_server <- function(id, data_list, method = c("volcano", "ma")) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    if (method != "volcano" & method != "ma") {
      stop(
        "The method parameter is misspecified. It must either be 'volcano' or 'ma'."
      )
    }

    if (method == "volcano") {
      choice_list <- list(
        "Absolute LFC" = 1,
        "-log10( Adjusted p-Val )" = 2,
        "Distance from the origin" = 3
      )
      name <- "-log10 (Adj. p-Val )"
    } else {
      choice_list <- list(
        "Absolute LFC" = 1,
        "Average Expression" = 2,
        "Distance from the origin" = 3
      )
      name <- "Average expression"
    }

    # pop up modal for selections ----
    observeEvent(input$customize_labels, {
      shiny::showModal(
        shiny::modalDialog(
          size = "m",
          p("Customize which genes are labeled."),
          selectInput(
            inputId = ns("gene_label_type"),
            label = "Gene selection",
            choices = list(
              "Do not label genes" = 1,
              "Label specific gene(s) from list" = 2,
              "Label top n genes" = 3,
              "Label genes above a certain threshold" = 4
            ),
            selected = 1
          ),
          conditionalPanel(
            condition = "input.gene_label_type == 2",
            selectInput(
              inputId = ns("vol_genes"),
              label = "Label Genes",
              choices = data_list()$anotate_genes$Row.names,
              multiple = TRUE,
              selectize = TRUE
            ),
            ns = ns
          ),
          conditionalPanel(
            condition = "input.gene_label_type == 3",
            fluidRow(
              column(
                width = 6,
                numericInput(
                  inputId = ns("num_genes"),
                  label = "Label top n genes",
                  min = 1,
                  max = 25,
                  value = 5
                )
              ),
              column(
                width = 6,
                selectInput(
                  inputId = ns("sort_type"),
                  label = "By",
                  choices = choice_list,
                  selected = 1
                )
              )
            ),
            ns = ns
          ),
          conditionalPanel(
            condition = "input.gene_label_type == 4",
            fluidRow(
              p("Label genes with"),
              column(
                width = 6,
                numericInput(
                  inputId = ns("min_lfc"),
                  label = "Absolute LFC greater than",
                  min = 2,
                  max = 20,
                  value = 3
                )
              ),
              column(
                width = 6,
                numericInput(
                  inputId = ns("min_value"),
                  label = paste0(name, " greater than"),
                  min = 5,
                  max = 60,
                  value = 20
                )
              )
            ),
            ns = ns
          ),
          p("Only genes that were identified as differently expressed can be labeled.")
        )
      )
    })

    # filter genes based on selections ----
    gene_labels <- reactive({
      req(data_list())
      data <- data_list()$data |>
        dplyr::filter(upOrDown != "None")

      if (method == "volcano") {
        data <- data |>
          dplyr::mutate(
            var = -log10(FDR),
            Fold = abs(Fold),
            dist = sqrt(Fold^2 + var^2)
          )
      } else {
        data <- data |>
          dplyr::mutate(
            var = Average,
            Fold = abs(Fold),
            dist = sqrt(Fold^2 + var^2)
          )
      }


      if (is.null(input$gene_label_type)) {
        genes <- NULL
      } else if (input$gene_label_type == 1) {
        genes <- NULL
      } else if (input$gene_label_type == 2) {
        genes <- input$vol_genes
      } else if (input$gene_label_type == 3) {
        if (input$sort_type == 1) {
          sorted <- data |>
            dplyr::arrange(dplyr::desc(Fold)) |>
            dplyr::slice(1:input$num_genes)
          genes <- sorted |>
            dplyr::pull(Row.names)
        } else if (input$sort_type == 2) {
          sorted <- data |>
            dplyr::arrange(dplyr::desc(var)) |>
            dplyr::slice(1:input$num_genes)
          genes <- sorted |>
            dplyr::pull(Row.names)
        } else if (input$sort_type == 3) {
          sorted <- data |>
            dplyr::arrange(dplyr::desc(dist)) |>
            dplyr::slice(1:input$num_genes)
          genes <- sorted |>
            dplyr::pull(Row.names)
        }
      } else if (input$gene_label_type == 4) {
        sorted <- data |>
          dplyr::filter(Fold >= input$min_lfc & var > input$min_value)
        genes <- sorted |>
          dplyr::pull(Row.names)
      }

      return(genes)
    })

    return(reactive({
      gene_labels()
    }))
  })
}
