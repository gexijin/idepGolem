#' fct_02_pre_process.R This file holds all of the main data analysis functions
#' associated with second tab of the iDEP website.
#'
#'
#' @section fct_02_pre_process.R functions:
#' \code{plot_genes}
#'
#'
#' @name fct_02_pre_process.R
NULL

#' Get x-axis label size based on number of samples
#'
#' @param n_samples Number of samples in the data
#'
#' @return Numeric value for x-axis label size
#'
#' @family preprocess functions
get_x_axis_label_size <- function(n_samples) {
  if (n_samples < 31) {
    return(16)
  } else if (n_samples <= 60) {
    return(12)
  } else if (n_samples <= 100) {
        return(10)
  } else {
    return(10)
  }
}

#' @title Pre-Process the data
#'
#' @description This function takes in user-defined values to
#' process the data for the EDA. Processing steps depend on data format, but
#' generally include missing value imputation, data filtering, and data
#' transformations.
#'
#' @param data Matrix of data that has already gone through
#'   \code{\link{convert_data}()}
#' @param missing_value String indicating method to deal with missing data. This
#'   should be one of "geneMedian", "treatAsZero", or "geneMedianInGroup"
#' @param data_file_format Integer indicating the data format. This should be
#'   one of 1 for read counts data, 2 for normalized expression, 3 for fold
#'   changes with adjusted P-values, or 4 for fold changes without P-values
#' @param low_filter_fpkm Integer for low count filter if
#'   \code{data_file_format} is normalized expression, \code{NULL} otherwise
#' @param n_min_samples_fpkm Integer for minimum samples if
#'   \code{data_file_format} is normalized expression, \code{NULL} otherwise
#' @param log_transform_fpkm TRUE/FALSE if a log transformation should be
#'   applied to normalized expression data
#' @param log_start_fpkm Integer added to log transformation if
#'   \code{data_file_format} is normalized expression, \code{NULL} otherwise
#' @param min_counts Numeric value for minimum count if
#'   \code{data_file_format} is read counts
#' @param n_min_samples_count Integer for minimum libraries with
#'   \code{min_counts} if \code{data_file_format} is read counts
#' @param counts_transform Integer to indicate which transformation to make if
#'   \code{data_file_format} is read counts. This should be one of 1 for
#'   log2(CPM+c) (EdgeR), 2 for variance stabilizing transformation (VST), or 3
#'   for regularized log (rlog)
#' @param counts_log_start Integer added to log if \code{counts_transform} is
#'   log2(CPM + 2)
#' @param no_fdr TRUE/FALSE to indicate fold-changes-only data with no p-values
#'   if \code{data_file_format} is fold changes
#'
#' @export
#' @return A list containing the transformed data, the mean kurtosis,
#' the raw counts, a data type warning, the size of the original data,
#' and p-values.
#'
#' @family preprocess functions
#' @seealso
#' * \code{\link[edgeR]{cpm}()} for information on calculating counts per
#'   million
#' * \code{\link[DESeq2]{vst}()} for information on variance
#'   stabilizing transformation
#' * \code{\link[DESeq2]{rlog}()} for
#'   information on the regularized log transformation
#' @md
#'
pre_process <- function(data,
                        missing_value = c("geneMedian", "treatAsZero", "geneMedianInGroup"),
                        data_file_format = c(1, 2, 3, 4),
                        low_filter_fpkm,
                        n_min_samples_fpkm,
                        log_transform_fpkm,
                        log_start_fpkm,
                        min_counts,
                        n_min_samples_count,
                        counts_transform,
                        counts_log_start,
                        no_fdr) {
  data_size_original <- dim(data)
  kurtosis_log <- 50

  results <- list(
    data = as.matrix(data),
    mean_kurtosis = NULL,
    raw_counts = NULL,
    data_type_warning = 0,
    data_size = c(data_size_original),
    p_vals = NULL,
    descr = generate_descr(
      missing_value,
      data_file_format,
      low_filter_fpkm,
      n_min_samples_fpkm,
      log_transform_fpkm,
      log_start_fpkm,
      min_counts,
      n_min_samples_count,
      counts_transform,
      counts_log_start,
      no_fdr
    )
  )

  # Sort by standard deviation -----------
  data <- data[order(-apply(
    data[, 1:dim(data)[2]],
    1,
    function(x) sd(x, na.rm = TRUE)
  )), ]

  # Missing values in expression data ----------
  if (sum(is.na(data)) > 0) {
    if (missing_value == "geneMedian") {
      row_medians <- apply(data, 1, function(y) median(y, na.rm = T))
      for (i in 1:ncol(data)) {
        val_miss_row <- which(is.na(data[, i]))
        data[val_miss_row, i] <- row_medians[val_miss_row]
      }
    } else if (missing_value == "treatAsZero") {
      data[is.na(data)] <- 0
    } else if (missing_value == "geneMedianInGroup") {
      sample_groups <- detect_groups(colnames(data), preserve_original = TRUE)
      for (group in unique(sample_groups)) {
        samples <- which(sample_groups == group)
        row_medians <- apply(
          data[, samples, drop = F],
          1,
          function(y) median(y, na.rm = T)
        )
        for (i in samples) {
          missing <- which(is.na(data[, i]))
          if (length(missing) > 0) {
            data[missing, i] <- row_medians[missing]
          }
        }
      }
      if (sum(is.na(data)) > 0) {
        row_medians <- apply(
          data,
          1,
          function(y) median(y, na.rm = T)
        )
        for (i in 1:ncol(data)) {
          missing <- which(is.na(data[, i]))
          data[missing, i] <- row_medians[missing]
        }
      }
    }
  }

  # Compute kurtosis ---------
  results$mean_kurtosis <- mean(apply(data, 2, e1071::kurtosis), na.rm = TRUE)

  # Pre-processing for each file format ----------
  if (data_file_format == 2) {
    if (is.integer(data)) {
      results$data_type_warning <- 1
    }

    # Filters ----------
    # Not enough counts
    data <- data[which(apply(
      data,
      1,
      function(y) sum(y >= low_filter_fpkm)
    ) >= n_min_samples_fpkm), ]

    # Same levels in every entry
    data <- data[which(apply(
      data,
      1,
      function(y) max(y) - min(y)
    ) > 0), ]

    # Takes log if log is selected OR kurtosis is bigger than 50
    if (
      (log_transform_fpkm == TRUE) ||
        (results$mean_kurtosis > kurtosis_log)
    ) {
      data <- log(data + abs(log_start_fpkm), 2)
    }

    std_dev <- apply(data, 1, sd)
    data <- data[order(-std_dev), ]
  } else if (data_file_format == 1) {
    if (!is.integer(data) && results$mean_kurtosis < kurtosis_log) {
      results$data_type_warning <- -1
    }

    data <- round(data, 0)
    # Check if any columns have all zeros
    if (any(apply(data, 2, function(col) all(col == 0)))) {
      results$data_type_warning <- -2
      return(results)
    }


    data <- data[which(apply(
      edgeR::cpm(edgeR::DGEList(counts = data)),
      1,
      function(y) sum(y >= min_counts)
    ) >= n_min_samples_count), ]

    # R cannot handle integers larger than 3 billion, which can impact popular 
    # packages such as DESeq2. R still uses 32-bit integers, and the largest 
    # allowable integer is 2^32 −1. In an unusual case, a user's RNA-Seq counts 
    # matrix included a count of 4 billion for a single gene, which was converted 
    # to NA, leading to an error in DESeq2. This issue also caused the iDEP app to crash. 
    if(max(data) > 2e9) {
      scale_factor <- max(data) / (2^32 - 1)
      # Round up scale_factor to the nearest integer
      scale_factor <- ceiling(scale_factor / 10 + 1) * 10 #  just to be safe.
      # Divide by the scale factor and round the entire matrix to the nearest integer
      data <- round(data / scale_factor)       
      # Warning message via the showNotification function
      showNotification(paste("Data includes counts bigger than 2.15 billion (2^32). The data was divided by ", scale_factor,". Double-check the data file. Is this really counts data? Proceed with caution."), duration = 60, type = "warning")
    }

    results$raw_counts <- data

    # Construct DESeqExpression Object ----------
    tem <- rep("A", dim(data)[2])
    tem[1] <- "B"
    col_data <- cbind(colnames(data), tem)
    colnames(col_data) <- c("sample", "groups")

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = data,
      colData = col_data,
      design = ~groups
    )
    dds <- DESeq2::estimateSizeFactors(dds)

    # Counts Transformation ------------
    if (counts_transform == 3) {
      data <- DESeq2::rlog(dds, blind = TRUE)
      data <- SummarizedExperiment::assay(data)
    } else if (counts_transform == 2) {
      data <- DESeq2::vst(dds, blind = TRUE)
      data <- SummarizedExperiment::assay(data)
    } else {
      data <- log2(BiocGenerics::counts(
        dds,
        normalized = TRUE
      ) + counts_log_start)
    }
  } else if (data_file_format %in% c(3, 4)) { # LFC summary statistics
    n2 <- (ncol(data) %/% 2)
    results$raw_counts <- data
    if (!no_fdr) {
      results$p_vals <- data[, 2 * (1:n2), drop = FALSE]
      data <- data[, 2 * (1:n2) - 1, drop = FALSE]
      if (ncol(data) == 1) {
        placeholder <- rep(1, dim(data)[1])
        results$p_vals <- cbind(results$p_vals, placeholder)
        Ctrl <- rep(0, dim(data)[1])
        data <- cbind(data, Ctrl)
      }
    }

    # Detect ratio data (fold-changes as ratios instead of log2)
    # Ratio data characteristics:
    # 1. No negative values
    # 2. High right skew (values > 1 for upregulation, 0-1 for downregulation)
    # If detected, apply log2 transformation
    fc_cols <- if (!no_fdr) {
      2 * (1:n2) - 1  # Fold-change columns (odd columns)
    } else {
      1:ncol(data)    # All columns are fold-changes
    }

    # Extract fold-change columns from raw data
    fc_data_raw <- results$raw_counts[, fc_cols, drop = FALSE]

    # Check each fold-change column for ratio characteristics
    for (i in seq_along(fc_cols)) {
      col_data <- fc_data_raw[, i]
      col_data <- col_data[!is.na(col_data)]  # Remove NAs for analysis

      if (length(col_data) > 0) {
        # Detect ratio data:
        # 1. No negative values (or very few, < 1%)
        # 2. High skewness to the right (indicating ratio scale)
        has_no_negatives <- sum(col_data < 0) / length(col_data) < 0.01

        # Calculate skewness: ratio data typically has skewness > 1
        # Ratio data: many values 0-1 (downreg), some values > 1 (upreg)
        # This creates right skew
        if (has_no_negatives && length(col_data) >= 10) {
          skewness <- e1071::skewness(col_data)

          # If skewness > 1 and no negatives, likely ratio data
          if (skewness > 1) {
            # Transform this column: log2(ratio)
            # Add small constant to avoid log(0)
            data[, i] <- log2(results$raw_counts[, fc_cols[i]] + 0.01)

            # Show notification to the user
            showNotification(
              paste0("Column ", colnames(results$raw_counts)[fc_cols[i]],
                     " detected as ratio data (no negatives, right skew = ",
                     round(skewness, 2), "). Applied log2 transformation."),
              duration = 10,
              type = "warning"
            )
          }
        }
      }
    }
  }
  results$data_size <- c(results$data_size, dim(data))

  validate(
    need(
      nrow(data) > 5 && ncol(data) >= 1,
      "Data file not recognized. Please double-check."
    )
  )

  data <- data[order(-apply(
    data[, 1:dim(data)[2]],
    1,
    sd
  )), ]

  results$data <- as.matrix(data)


  return(results)
}

#' Creates a barplot of the count data
#'
#' This function takes in either raw count or processed data and creates a
#' formatted barplot as a \code{ggplot} object that shows the number
#' of genes mapped to each sample in millions. This function is only used for
#' read counts data.
#'
#' @param counts_data Matrix of raw counts from gene expression data
#' @param sample_info Matrix of experiment design information for grouping
#'  samples
#' @param type String designating the type of data to be used in the title.
#'  Commonly either "Raw" or "Transformed"
#' @param plots_color_select Vector of colors for plots
#'
#' @export
#' @return A list containing 'plot' (a ggplot object) and 'warning' (a
#'  character string with ANOVA warning message, or NULL if no significant
#'  differences detected)
#'
#' @family preprocess functions
#' @family plots
#'
#'
total_counts_ggplot <- function(counts_data,
                                sample_info,
                                type = "",
                                plots_color_select) {
  counts <- counts_data

  groups <- as.factor(
    detect_groups(colnames(counts), sample_info)
  )

  x_axis_labels <- get_x_axis_label_size(ncol(counts))

  if (length(unique(groups)) <= 1 || length(unique(groups)) > 20) {
    plot_data <- data.frame(
      sample = as.factor(colnames(counts)),
      counts = colSums(counts) / 1e6,
      group = groups
    )

    plot <- ggplot2::ggplot(
      data = plot_data,
      ggplot2::aes(x = sample, y = counts)
    )
  } else {
    grouping <- groups

    color_palette <- generate_colors(n = length(unique(grouping)), palette_name = plots_color_select)

    plot_data <- data.frame(
      sample = as.factor(colnames(counts)),
      counts = colSums(counts) / 1e6,
      group = groups,
      grouping = grouping
    )

    plot <- ggplot2::ggplot(
      data = plot_data,
      ggplot2::aes(x = sample, y = counts, fill = grouping)
    ) +
    ggplot2::scale_fill_manual(values = color_palette)
  }

  plot <- plot +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
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
      title = paste("Total", type, "Read Counts (Millions)"),
      y = paste(type, "Counts (Millions)"),
      fill = NULL
    )

  return(plot)
}



#' Creates a barplot of the count of genes by type
#'
#' This function works with either raw count or normalized expression matrices.
#' It links the supplied gene metadata to tally how many genes fall into each
#' gene type and returns the counts as a bar chart.
#'
#' @param counts_data Matrix of expression values used to align genes with
#'  metadata (values themselves are ignored)
#' @param sample_info Matrix of experiment design information for grouping
#'  samples
#' @param type String designating the type of data to be used in the title.
#'  Commonly either "Raw" or "Transformed"
#' @param all_gene_info Gene info, including chr., gene type etc.
#' @param plots_color_select Vector of colors for plots
#' @export
#' @return A barplot as a \code{ggplot} object
#'
#' @family preprocess functions
#' @family plots
#'
#'
gene_counts_ggplot <- function(counts_data,
                                sample_info,
                                type = "",
                                all_gene_info,
                                plots_color_select) {
  df <- merge(
    counts_data,
    all_gene_info,
    by.x = "row.names",
    by.y = "ensembl_gene_id"
  )
  df$gene_biotype <- gsub(".*pseudogene", "Pseudogene", df$gene_biotype)
  df$gene_biotype <- gsub("TEC", "Unknown", df$gene_biotype)
  df$gene_biotype <- gsub("artifact", "Artifact", df$gene_biotype)
  df$gene_biotype <- gsub("IG_.*", "IG", df$gene_biotype)
  df$gene_biotype <- gsub("TR_.*", "TR", df$gene_biotype)
  df$gene_biotype <- gsub("protein_coding", "Coding", df$gene_biotype)

  counts <- table(df$gene_biotype)
  # Convert the vector to a dataframe
  data <- data.frame(
    category = names(counts),
    value = as.numeric(counts)
  )
  # Order the categories by value
  data <- data[order(-data$value), ]

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = reorder(category, value), y = value, fill = category)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_y_log10(limits = c(1, 2 * max(data$value))) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Number of genes") +       
    ggplot2::geom_text(ggplot2::aes(label = value), hjust = -0.1, vjust = 0.5) 

  plot <- plot +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 16
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
    ) 


  return(plot)
}



#' Creates a barplot of the count data by gene type
#'
#' This function takes in either raw count or processed data and creates a
#' formatted barplot as a \code{ggplot} object that shows the number
#' of genes mapped to each sample in millions. This function is only used for
#' read counts data.
#'
#' @param counts_data Matrix of raw counts from gene expression data
#' @param sample_info Matrix of experiment design information for grouping
#'  samples
#' @param type String designating the type of data to be used in the title.
#'  Commonly either "Raw" or "Transformed"
#' @param all_gene_info Gene info, including chr., gene type etc.
#' @param plots_color_select Vector of colors for plots
#' @export
#' @return A barplot as a \code{ggplot} object
#'
#' @family preprocess functions
#' @family plots
#'
#'
rRNA_counts_ggplot <- function(counts_data,
                                sample_info,
                                type = "",
                                all_gene_info,
                                plots_color_select) {
  counts <- counts_data

  groups <- as.factor(
    detect_groups(colnames(counts), sample_info)
  )

  x_axis_labels <- get_x_axis_label_size(ncol(counts))

  df <- merge(
    counts_data,
    all_gene_info,
    by.x = "row.names",
    by.y = "ensembl_gene_id"
  )
  df$gene_biotype <- gsub(".*pseudogene", "Pseudogene", df$gene_biotype)
  df$gene_biotype <- gsub("TEC", "Unknown", df$gene_biotype)
  df$gene_biotype <- gsub("IG_.*", "IG", df$gene_biotype)
  df$gene_biotype <- gsub("TR_.*", "TR", df$gene_biotype)
  df$gene_biotype <- gsub("protein_coding", "Coding", df$gene_biotype)

  counts_by_type <- aggregate(
    df[, colnames(counts_data)],
    by = list(df$gene_biotype),
    FUN = sum
  )
  colnames(counts_by_type)[1] = "Gene_Type"

  df <- cbind(Gene_Type = counts_by_type[, 1], sweep(counts_by_type[-1], 2, 0.01 * colSums(counts_by_type[-1]), "/"))

  # remove categories less than 0.5%
  df <- df[which(apply(df[, -1], 1, max) > 0.5), ]

  plot_data <- reshape2::melt(df, id.vars = "Gene_Type")
  plot_data$groups <- groups[match(plot_data$variable, colnames(counts_data))]

  avg_pct <- rowMeans(df[, -1, drop = FALSE], na.rm = TRUE)
  gene_type_order <- df$Gene_Type[order(avg_pct, decreasing = TRUE)]
  plot_data$Gene_Type <- factor(plot_data$Gene_Type, levels = gene_type_order)

  color_palette <- generate_colors(n = length(gene_type_order), palette_name = plots_color_select)
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = variable, y = value, fill = Gene_Type)) +
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_stack(reverse = TRUE)
    ) +
    ggplot2::labs(
      x = NULL,
      y = "% Reads",
      title = "% Reads by gene type",
      fill = NULL
    ) +
    ggplot2::scale_fill_manual(values = color_palette, drop = FALSE) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))

  plot <- plot +
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_stack(reverse = TRUE)
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
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
    ) 

  # Run ANOVA analysis if groups are defined and less than half of samples
  warning_message <- NULL
  n_samples <- length(groups)
  n_groups <- length(unique(groups))
  groups_well_defined <- (n_groups > 1) && (n_groups <= 20) && (n_samples >= 2 * n_groups)

  if (groups_well_defined && (n_groups < n_samples / 2)) {
    significant_types <- c()

    for (gene_type in unique(plot_data$Gene_Type)) {
      type_data <- plot_data[plot_data$Gene_Type == gene_type, ]
      type_data <- type_data[!is.na(type_data$groups), ]

      if (nrow(type_data) > 0 && length(unique(type_data$groups)) > 1) {
        tryCatch({
          anova_result <- summary(aov(value ~ groups, data = type_data))
          p_value <- anova_result[[1]][["Pr(>F)"]][1]

          group_means <- aggregate(value ~ groups, data = type_data, mean)
          max_mean <- max(group_means$value)
          min_mean <- min(group_means$value)
          diff_percent <- abs(max_mean - min_mean)
          # at least 5% difference and p < 0.01
          if (!is.na(p_value) && p_value < 0.01 && diff_percent > 5) {
            significant_types <- c(
              significant_types,
              paste0(gene_type, " (P = ", format(p_value, scientific = TRUE, digits = 2), ")")
            )
          }
        }, error = function(e) {
          # Skip if ANOVA fails
        })
      }
    }

    if (length(significant_types) > 0) {
      warning_message <- paste0(
        "Significant difference detected for % reads mapped across sample groups for gene type(s): ",
        paste(significant_types, collapse = ", "),
        "."
      )
    }
  }

  return(list(plot = plot, warning = warning_message))
}



#' Creates a barplot or boxplot of the count data by chr
#'
#' This function takes in either raw count or processed data and creates a
#' formatted plot as a \code{ggplot} object that shows the percentage of
#' reads mapped to each chromosome. By default shows barplot; can use boxplot
#' with jitter when user selects and sample groups are well-defined.
#'
#' @param counts_data Matrix of raw counts from gene expression data
#' @param sample_info Matrix of experiment design information for grouping
#'  samples
#' @param type String designating the type of data to be used in the title.
#'  Commonly either "Raw" or "Transformed"
#' @param all_gene_info Gene info, including chr., gene type etc.
#' @param plots_color_select Vector of colors for plots (used for boxplot groups)
#' @param use_boxplot Logical indicating whether to use boxplot (TRUE) or
#'  barplot (FALSE, default)
#' @export
#' @return A list containing 'plot' (a ggplot object) and 'warning' (a
#'  character string with ANOVA warning message, or NULL if no significant
#'  differences detected)
#'
#' @family preprocess functions
#' @family plots
#'
#'
chr_counts_ggplot <- function(counts_data,
                                sample_info,
                                type = "",
                                all_gene_info,
                                plots_color_select = "Set1",
                                use_boxplot = FALSE) {
  counts <- counts_data

  groups <- as.factor(
    detect_groups(colnames(counts), sample_info)
  )

  x_axis_labels <- get_x_axis_label_size(ncol(counts))

  df <- merge(
    counts_data,
    all_gene_info,
    by.x = "row.names",
    by.y = "ensembl_gene_id"
  )
  colnames(df)[which(colnames(df) == "Row.names")] <- "ensembl_gene_id"
  # only coding genes?
  ignore_non_coding = FALSE
  if (ignore_non_coding) {
    df <- subset(df, gene_biotype == "protein_coding")
  }

  # If no chromosomes found. For example if user do not convert gene IDs.
  if (dim(df)[1] < 5) {
    return(NULL)
  }

  # gather info on chrosomes
  tem <- sort(table(df$chromosome_name), decreasing = TRUE)
  # ch with less than 100 genes are excluded
  hide_chr <- TRUE
  if (hide_chr) {
    ch <- names(tem[tem >= 5])
  } else {
    ch <- names(tem[tem >= 1])
  }

  if (length(ch) > 50) {
    # At most 50 ch
    ch <- ch[1:50]
  }

  ch <- ch[order(as.numeric(ch))]
  tem <- ch
  # The numbers are continous from 1 to length(ch)
  ch <- 1:(length(ch))
  # The names are real chr. names
  names(ch) <- tem

  nchar_cutoff <- 3 * quantile(nchar(tem), .25)

  # remove long chr. patch.., longer than 3 times of the length of the first quantile
  ch <- ch[nchar(tem) < nchar_cutoff]

  df <- df[which(df$chromosome_name %in% names(ch)), ]
  df <- droplevels(df)
  # Numeric encoding
  df$chNum <- 1
  df$chNum <- ch[df$chromosome_name]
  # change order of chr.
  df$chromosome_name <- factor(df$chromosome_name, levels = names(ch))

  counts_by_chr <- aggregate(
    df[, colnames(counts_data)],
    by = list(df$chromosome_name),
    FUN = sum
  )
  colnames(counts_by_chr)[1] = "Chr"

  df <- cbind(
    Chr = counts_by_chr[, 1],
    sweep(
      counts_by_chr[-1],
      2,
      0.01 * colSums(counts_by_chr[-1]), "/"
    )
  )

  # remove categories less than 0.5%
  df <- df[which(apply(df[, -1], 1, max) > 0.01), ]

  plot_data <- reshape2::melt(df, id.vars = "Chr")

  # Determine if groups are well-defined
  # Only use boxplot if user requested it AND groups are well-defined
  n_samples <- ncol(counts)
  n_groups <- length(unique(groups))
  groups_well_defined <- (n_groups > 1) && (n_groups <= 20) && (n_samples >= 2 * n_groups)
  show_boxplot <- use_boxplot && groups_well_defined

  if (show_boxplot) {
    # Add group information to plot_data
    plot_data$groups <- groups[match(plot_data$variable, colnames(counts))]

    # Generate color palette for groups
    color_palette <- generate_colors(n = n_groups, palette_name = plots_color_select)

    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = groups, y = value, fill = groups)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2) +
      ggplot2::labs(
        x = NULL,
        y = "% Reads",
        title = "% Reads by Chromosomes",
        fill = NULL
      ) +
      ggplot2::scale_fill_manual(values = color_palette)
  } else {
    # Use original barplot
    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = variable, y = value)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::labs(x = NULL, y = "% Reads", title = "% Reads by Chromosomes") +
      ggplot2::scale_fill_brewer(palette = "Set1")
  }

  if(ncol(counts_data) < 10) {
    plot <- plot + ggplot2::facet_wrap (. ~ Chr)
  } else if (ncol(counts_data) < 20){
    plot <- plot + ggplot2::facet_wrap (. ~ Chr, ncol = 3)
  } else if (ncol(counts_data) < 40){
    plot <- plot + ggplot2::facet_wrap (. ~ Chr, ncol = 2)
  } else {
    plot <- plot + ggplot2::facet_wrap (. ~ Chr, ncol = 1)
  }

  if (!show_boxplot) {
    plot <- plot + ggplot2::geom_bar(stat = "identity")
  }

  plot <- plot +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.y = ggplot2::element_text(
        color = "black",
        face = "bold",
        size = 20
      ),
      axis.text.y = ggplot2::element_text(
        size = x_axis_labels
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      ),
      strip.text = ggplot2::element_text(
        size = 16,
        color = "black",
        face = "bold"
        )
    )

  # Run ANOVA analysis if groups are defined and less than half of samples
  warning_message <- NULL
  if (groups_well_defined && (n_groups < n_samples / 2)) {
    # Perform ANOVA for each chromosome
    significant_chrs <- c()

    for (chr in unique(df$Chr)) {
      chr_data <- df[df$Chr == chr, -1]  # Exclude Chr column
      chr_values <- as.numeric(chr_data)
      chr_groups <- rep(groups, each = 1)

      # Only run ANOVA if we have enough data
      if (length(chr_values) > 0 && length(unique(chr_groups)) > 1) {
        tryCatch({
          anova_result <- summary(aov(chr_values ~ chr_groups))
          p_value <- anova_result[[1]][["Pr(>F)"]][1]

          # Calculate group means
          group_means <- aggregate(chr_values, by = list(chr_groups), mean)
          max_mean <- max(group_means$x)
          min_mean <- min(group_means$x)
          diff_expression <- abs(max_mean - min_mean)

          # Check if significant and difference > 5 units (% reads)
          if (!is.na(p_value) && p_value < 0.01 && diff_expression > 5) {
            significant_chrs <- c(
              significant_chrs,
                paste0(chr, " (P = ", format(p_value, scientific = TRUE, digits = 2), ")")
            )
          }
        }, error = function(e) {
          # Skip if ANOVA fails
        })
      }
    }

    # Build warning message if significant chromosomes found
    if (length(significant_chrs) > 0) {
      warning_message <- paste0(
        "Sgnificant difference detected for % reads mapped ",
        "to chromosome(s): ",
        paste(significant_chrs, collapse = ", "),
        "."
      )
    }
  }

  return(list(plot = plot, warning = warning_message))
}



#' Creates a barplot or boxplot of the normalized data by chr
#'
#' This function takes in either raw count or processed data and creates a
#' formatted plot as a \code{ggplot} object that shows the 75th percentile of
#' normalized expression by chromosome. By default shows barplot; can use boxplot
#' with jitter when user selects and sample groups are well-defined.
#'
#' @param counts_data Matrix of raw counts from gene expression data
#' @param sample_info Matrix of experiment design information for grouping
#'  samples
#' @param type String designating the type of data to be used in the title.
#'  Commonly either "Raw" or "Transformed"
#' @param all_gene_info Gene info, including chr., gene type etc.
#' @param plots_color_select Vector of colors for plots (used for boxplot groups)
#' @param use_boxplot Logical indicating whether to use boxplot (TRUE) or
#'  barplot (FALSE, default)
#' @export
#' @return A list containing 'plot' (a ggplot object) and 'warning' (a
#'  character string with ANOVA warning message, or NULL if no significant
#'  differences detected)
#'
#' @family preprocess functions
#' @family plots
#'
#'
chr_normalized_ggplot <- function(counts_data,
                                sample_info,
                                type = "",
                                all_gene_info,
                                plots_color_select = "Set1",
                                use_boxplot = FALSE) {
  counts <- counts_data

  groups <- as.factor(
    detect_groups(colnames(counts), sample_info)
  )

  x_axis_labels <- get_x_axis_label_size(ncol(counts))

  df <- merge(
    counts_data,
    all_gene_info,
    by.x = "row.names",
    by.y = "ensembl_gene_id"
  )
  colnames(df)[which(colnames(df) == "Row.names")] <- "ensembl_gene_id"
  # only coding genes?
  ignore_non_coding = FALSE
  if (ignore_non_coding) {
    df <- subset(df, gene_biotype == "protein_coding")
  }

  # If no chromosomes found. For example if user do not convert gene IDs.
  if (dim(df)[1] < 5) {
    return(NULL)
  }

  # gather info on chrosomes
  tem <- sort(table(df$chromosome_name), decreasing = TRUE)
  # ch with less than 100 genes are excluded
  hide_chr <- TRUE
  if (hide_chr) {
    ch <- names(tem[tem >= 5])
  } else {
    ch <- names(tem[tem >= 1])
  }

  if (length(ch) > 50) {
    # At most 50 ch
    ch <- ch[1:50]
  }

  ch <- ch[order(as.numeric(ch))]
  tem <- ch
  # The numbers are continous from 1 to length(ch)
  ch <- 1:(length(ch))
  # The names are real chr. names
  names(ch) <- tem

  nchar_cutoff <- 3 * quantile(nchar(tem), .25)

  # remove long chr. patch.., longer than 3 times of the length of the first quantile
  ch <- ch[nchar(tem) < nchar_cutoff]

  df <- df[which(df$chromosome_name %in% names(ch)), ]
  df <- droplevels(df)
  # Numeric encoding
  df$chNum <- 1
  df$chNum <- ch[df$chromosome_name]
  # change order of chr.
  df$chromosome_name <- factor(df$chromosome_name, levels = names(ch))

  counts_by_chr <- aggregate(
    df[, colnames(counts_data)],
    by = list(df$chromosome_name),
    FUN = function(x) {
      quantile(x, 0.75, na.rm = TRUE)
    }
  )
  colnames(counts_by_chr)[1] = "Chr"

  df <- cbind(
    Chr = counts_by_chr[, 1], 
    sweep(
      counts_by_chr[-1], 
      2, 
      0.01 * colSums(counts_by_chr[-1]), "/"
    )
  )

  # remove categories less than 0.5%
  df <- df[which(apply(df[, -1], 1, max) > 0.1), ]

  plot_data <- reshape2::melt(df, id.vars = "Chr")

  # Determine if groups are well-defined
  # Only use boxplot if user requested it AND groups are well-defined
  n_samples <- ncol(counts)
  n_groups <- length(unique(groups))
  groups_well_defined <- (n_groups > 1) && (n_groups <= 20) && (n_samples >= 2 * n_groups)
  show_boxplot <- use_boxplot && groups_well_defined

  if (show_boxplot) {
    # Add group information to plot_data
    plot_data$groups <- groups[match(plot_data$variable, colnames(counts))]

    # Generate color palette for groups
    color_palette <- generate_colors(n = n_groups, palette_name = plots_color_select)

    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = groups, y = value, fill = groups)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_jitter(width = 0.2, height = 0, alpha = 0.6, size = 2) +
      ggplot2::labs(
        x = NULL,
        y = "Normalized Expression",
        title = "75th percentile of normalized expression by chromosomes",
        fill = NULL
      ) +
      ggplot2::scale_fill_manual(values = color_palette)
  } else {
    # Use original barplot
    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = variable, y = value)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::labs(x = NULL, y = "Normalized Expression", title = "75th percentile of normalized expression by chromosomes") +
      ggplot2::scale_fill_brewer(palette = "Set1")
  }

  if(ncol(counts_data) < 10) {
    plot <- plot + ggplot2::facet_wrap (. ~ Chr)
  } else if (ncol(counts_data) < 20){
    plot <- plot + ggplot2::facet_wrap (. ~ Chr, ncol = 3)
  } else if (ncol(counts_data) < 40){
    plot <- plot + ggplot2::facet_wrap (. ~ Chr, ncol = 2)
  } else {
    plot <- plot + ggplot2::facet_wrap (. ~ Chr, ncol = 1)
  }

  if (!show_boxplot) {
    plot <- plot + ggplot2::geom_bar(stat = "identity")
  }

  plot <- plot +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.y = ggplot2::element_text(
        color = "black",
        face = "bold",
        size = 20
      ),
      axis.text.y = ggplot2::element_text(
        size = x_axis_labels
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
      ),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      ),
      strip.text = ggplot2::element_text(
        size = 16,
        color = "black",
        face = "bold"
        )
    )

  # Run ANOVA analysis if groups are defined and less than half of samples
  warning_message <- NULL
  if (groups_well_defined && (n_groups < n_samples / 2)) {
    # Perform ANOVA for each chromosome
    significant_chrs <- c()

    for (chr in unique(df$Chr)) {
      chr_data <- df[df$Chr == chr, -1]  # Exclude Chr column
      chr_values <- as.numeric(chr_data)
      chr_groups <- rep(groups, each = 1)

      # Only run ANOVA if we have enough data
      if (length(chr_values) > 0 && length(unique(chr_groups)) > 1) {
        tryCatch({
          anova_result <- summary(aov(chr_values ~ chr_groups))
          p_value <- anova_result[[1]][["Pr(>F)"]][1]

          # Calculate group means
          group_means <- aggregate(chr_values, by = list(chr_groups), mean)
          max_mean <- max(group_means$x)
          min_mean <- min(group_means$x)
          diff_expression <- abs(max_mean - min_mean)

          # Check if significant and difference > 0.5 units
          if (!is.na(p_value) && p_value < 0.01 && diff_expression > 0.5) {
            significant_chrs <- c(
              significant_chrs,
                paste0(chr, " (P = ", format(p_value, scientific = TRUE, digits = 2), ")")
            )
          }
        }, error = function(e) {
          # Skip if ANOVA fails
        })
      }
    }

    # Build warning message if significant chromosomes found
    if (length(significant_chrs) > 0) {
      warning_message <- paste0(
        "Significant difference detected expression levels among genes on chromosome(s): ",
        "on chromosome(s): ",
        paste(significant_chrs, collapse = ", "),
        "."
      )
    }
  }

  return(list(plot = plot, warning = warning_message))
}




#' Scatterplot for EDA on processed data
#'
#' This function takes the data after it has been pre-processed and
#' creates a scatterplot of the counts for two samples that are
#' indicated by the user.
#'
#' @param processed_data Matrix of data that has gone through
#'  \code{\link{prep_process}()}
#' @param plot_xaxis Character string indicating sample to plot on the x-axis
#' @param plot_yaxis Character string indicating sample to plot on the y axis
#'
#' @export
#' @return A scatterplot a \code{ggplot} object
#'
#' @family plots
#' @family preprocess functions
eda_scatter <- function(processed_data,
                        plot_xaxis,
                        plot_yaxis) {
  plot_data <- as.data.frame(processed_data)
  scatter <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = get(plot_xaxis),
      y = get(plot_yaxis)
    )
  ) +
    ggplot2::geom_point(size = 1) +
    ggplot2::labs(
      title = paste0(
        "Scatter for ", plot_xaxis, " and ", plot_yaxis,
        " transfromed expression"
      ),
      x = paste0(plot_xaxis),
      y = paste0(plot_yaxis)
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "none",
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 16),
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      )
    ) +
    ggpubr::stat_cor(label.x.npc = "left", label.y.npc = "top",
      p.accuracy = 0.01, r.accuracy = 0.01)

  return(scatter)
}


#' Boxplot for processed data
#'
#' This function takes the processed data and creates
#' a boxplot of number of sequences mapped to each
#' tissue sample.
#'
#' @param processed_data Matrix of data that has gone through
#'  \code{\link{pre_process}()}
#' @param sample_info Matrix of experiment design information
#' @param plots_color_select Vector of colors for plots
#'
#' @export
#' @return Boxplot of the distribution of counts for each sample as a
#'  \code{ggplot} object.
#'
#' @family plots
#' @family preprocess functions
#'
eda_boxplot <- function(processed_data,
                        sample_info,
                        plots_color_select) {
  counts <- as.data.frame(processed_data)

  groups <- as.factor(
    detect_groups(colnames(counts), sample_info)
  )
  grouping <- groups
  x_axis_labels <- get_x_axis_label_size(ncol(counts))

  longer_data <- tidyr::pivot_longer(
    data = counts,
    colnames(counts),
    names_to = "sample",
    values_to = "expression"
  )

  longer_data$groups <- rep(groups, nrow(counts))
  longer_data$grouping <- rep(grouping, nrow(counts))

  color_palette <- generate_colors(n = length(unique(grouping)), palette_name = plots_color_select)

  plot <- ggplot2::ggplot(
    data = longer_data,
    ggplot2::aes(x = sample, y = expression, fill = grouping)
  ) +
    ggplot2::scale_fill_manual(values = color_palette) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = "right",
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        size = x_axis_labels
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
      title = "Distribution of Transformed Data",
      y = "Transformed Expression",
      fill = NULL
    )

  return(plot)
}

#' Density plot for the processed data
#'
#' This function takes in the processed data and sample info
#' and creates a density plot for the distribution of sequences
#' that are mapped to each sample.
#'
#' @param processed_data Matrix of gene data that has gone through
#'   \code{\link{pre_process}()}
#' @param sample_info Sample_info from the experiment file
#' @param plots_color_select Vector of colors for plots
#'
#' @export
#' @return A density plot as a \code{ggplot} object
#'
#' @family plots
#' @family preprocess functions
eda_density <- function(processed_data,
                        sample_info,
                        plots_color_select) {
  counts <- as.data.frame(processed_data)

  groups <- as.factor(
    detect_groups(colnames(counts), sample_info)
  )
  group_fill <- groups
  legend <- if (nlevels(groups) <= 1) "none" else "right"
  x_axis_labels <- get_x_axis_label_size(ncol(counts))

  longer_data <- tidyr::pivot_longer(
    data = counts,
    colnames(counts),
    names_to = "sample",
    values_to = "expression"
  )
  longer_data$groups <- rep(groups, nrow(counts))
  longer_data$group_fill <- rep(group_fill, nrow(counts))

  color_palette <- generate_colors(n = length(unique(group_fill)), palette_name = plots_color_select)

  plot <- ggplot2::ggplot(
    data = longer_data,
    ggplot2::aes(x = expression, color = group_fill, group = sample)
  ) +
    ggplot2::geom_density(size = 1) +
    ggplot2::theme_light() +
    ggplot2::theme(
      legend.position = legend,
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.text.x = ggplot2::element_text(
        size = x_axis_labels
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
      title = "Density Plot of Transformed Data",
      x = "Transformed Expression",
      y = "Density",
      color = "Sample"
    ) +
    ggplot2::scale_color_manual(values = color_palette)

  return(plot)
}



#' Individual plotting function for genes
#'
#' Takes in the merged data and other data to provide
#' plots on individual gene names. Depening on the selections
#' this function will either return a barplot from the
#' grouped data or a lineplot from the individual sample.
#'
#' @param individual_data Data that has been merged with the gene info
#' @param sample_info Matrix of experiment design information for
#'   grouping samples
#' @param gene_plot_box TRUE/FALSE for individual sample plot or grouped
#'   data plot
#' @param use_sd TRUE/FALSE for standard error or standard deviation bars on
#'   bar plot
#' @param lab_rotate Numeric value indicating what angle to rotate
#'   the x-axis labels
#' @param plots_color_select Vector of colors for plots
#' @param selected_gene list of genes selected for plotting
#' @param plot_raw logical value specifying whether raw count data should be 
#' used
#' @param plot_tukey logical value that specifies whether or not to run 
#' TukeyHSD tests within genes between sample groups
#'
#' @export
#'
#' @return A \code{ggplot} object. For \code{gene_plot_box = TRUE} the return
#'  will be a lineplot for the expression of each individual sample for the
#'  selected gene. If \code{gene_plot_box = FALSE} the return will be a barplot
#'  for the groups provided in the sample information.
#'
#' @seealso \code{\link{merge_data}()}
#' @family plots
#' @family preprocess functions
#'
individual_plots <- function(individual_data,
                             sample_info,
                             selected_gene,
                             gene_plot_box,
                             use_sd,
                             lab_rotate,
                             plots_color_select,
                             plot_raw,
                             plot_tukey,
                             max_groups = 12,
                             max_length = 30) {
  individual_data <- as.data.frame(individual_data)
  individual_data$symbol <- rownames(individual_data)
  
  plot_data <- individual_data |>
    dplyr::filter(symbol %in% selected_gene) |>
    tidyr::pivot_longer(!symbol, names_to = "sample", values_to = "value")
  
  sample_count <- length(unique(plot_data$sample))
  gene_count <- length(unique(plot_data$symbol))
  x_axis_labels <- if (gene_plot_box == 2) {
    get_x_axis_label_size(sample_count)
  } else {
    get_x_axis_label_size(gene_count)
  }

  if (gene_plot_box == 2) {
    plot_data$symbol <- factor(plot_data$symbol, levels = unique(plot_data$symbol))
    color_palette <- generate_colors(n = nlevels(as.factor(plot_data$symbol)), palette_name = plots_color_select)
    ind_line <- ggplot2::ggplot(
      data = plot_data,
      ggplot2::aes(x = sample, y = value, group = symbol, color = symbol)
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 5, fill = "white") +
      ggplot2::labs(
        y = ifelse(plot_raw, "Raw counts", "Normalized Expression")
      ) +
      ggplot2::coord_cartesian(ylim = c(0, max(plot_data$value))) +
      ggplot2::theme_light() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          color = "black",
          size = 16,
          face = "bold",
          hjust = .5
        ),
        axis.text.x = ggplot2::element_text(
          angle = as.numeric(lab_rotate),
          size = x_axis_labels,
          vjust = .5
        ),
        axis.text.y = ggplot2::element_text(size = 16),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(
          color = "black",
          size = 14
        ),
        legend.text = ggplot2::element_text(size = 12)
      ) +
    ggplot2::scale_color_manual(values = color_palette)

    return(ind_line)
  } else if (gene_plot_box == 1) {
    # Get unique sample names first to avoid duplicates when multiple genes selected
    unique_samples <- unique(plot_data$sample)
    sample_groups <- detect_groups(unique_samples, sample_info, max_groups, max_length)
    # Map groups back to all rows in plot_data
    sample_to_group <- setNames(sample_groups, unique_samples)
    plot_data$groups <- sample_to_group[plot_data$sample]

    summarized <- plot_data |>
      dplyr::group_by(groups, symbol) |>
      dplyr::summarise(Mean = mean(value), SD = sd(value), N = dplyr::n())

    summarized$SE <- summarized$SD / sqrt(summarized$N)
    
    color_palette <- generate_colors(n = length(unique(summarized$groups)), palette_name = plots_color_select)
  
    gene_bar <- ggplot2::ggplot(
      summarized,
      ggplot2::aes(x = symbol, y = Mean, fill = groups)
    ) +
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_dodge()
      ) +
    ggplot2::labs(
      y = ifelse(plot_raw, "Raw counts", "Normalized Expression"),
      fill = NULL
    ) +
    ggplot2::geom_dotplot(
      data = plot_data,
      ggplot2::aes(
        y = value,
        groups = groups
      ),
      fill = "black",
      position = ggplot2::position_dodge(),
      binaxis='y',
      stackdir='center',
      dotsize=0.5
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      ),
      axis.text.x = ggplot2::element_text(
        angle = as.numeric(lab_rotate),
        size = x_axis_labels,
        vjust = .5
      ),
      axis.text.y = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      legend.text = ggplot2::element_text(size = 12)
    ) +
    ggplot2::scale_fill_manual(values = color_palette)

    if (use_sd == TRUE) {
      gene_bar <- gene_bar + ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = Mean - SD,
          ymax = Mean + SD
        ),
        width = 0.2,
        position = ggplot2::position_dodge(.9)
      )
    } else {
      gene_bar <- gene_bar + ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = Mean - SE,
          ymax = Mean + SE
        ),
        width = 0.2,
        position = ggplot2::position_dodge(.9)
      )
    }
    
    if (plot_tukey == TRUE) {
      dfAOV <- get_tukey_data(gene_data = plot_data,
                              sample_info = sample_info,
                              summarized = TRUE)
      
      dfAOV <- dfAOV |>
        dplyr::group_by(symbol) |>
        dplyr::arrange(pval) |> # Sort p-values, take lowest 10 in gene
        dplyr::slice(1:10) |>
        dplyr::mutate(n_bracket = dplyr::row_number())
      
      # Get group levels and offsets for dodging
      group_levels <- levels(as.factor(plot_data$groups))
      dodge_width <- 0.9
      n_groups <- length(group_levels)
      bar_width <- dodge_width / n_groups
      offsets <- setNames(
        seq(-(dodge_width / 2) + bar_width/2, 
            (dodge_width / 2) - bar_width/2, 
            length.out = n_groups), 
        group_levels)
      
      # Filter only significant p-values and apply offset for brackets
      dfAOV <- dfAOV |>
        dplyr::filter(pval <= 0.05) |>
        dplyr::ungroup() |>
        dplyr::mutate(
          pval = dplyr::case_when(pval <= 0.05 & pval > 0.01 ~ "*",
                                  pval <= 0.01 & pval > 0.001 ~ "**",
                                  pval <= 0.001 ~ "***"),
          g1 = factor(gsub("-.*", "", comp), levels = group_levels),
          g2 = factor(gsub(".*-", "", comp), levels = group_levels)
        ) |>
        dplyr::mutate(
          symbol_num = as.numeric(as.factor(symbol)),
          g1_offset = offsets[as.character(g1)],
          g2_offset = offsets[as.character(g2)],
          edge1 = symbol_num + g1_offset,
          edge2 = symbol_num + g2_offset,
          label_loc = (edge1 + edge2) / 2
        )
     
      #Plot brackets and p-values above bars 
      gene_bar <- gene_bar +
        ggplot2::geom_errorbarh(data = dfAOV, 
                                ggplot2::aes(
                                  xmin = edge1, 
                                  xmax = edge2,
                                  y = avg + (max(avg)* 0.1) * n_bracket,
                                  height = (max(avg)* 0.05)),
                                inherit.aes = FALSE)+
        ggplot2::geom_text(data = dfAOV,
                           mapping = ggplot2::aes(
                             x = label_loc, 
                             y = avg + ((max(avg)* 0.1) * n_bracket),
                             label = pval,
                             vjust = -0.45),
                           size = 4,
                           inherit.aes = FALSE)
    }
    
    return(gene_bar)
  }
}

#' Get TukeyHSD Result Data
#' 
#' Takes a data frame of gene expression data and runs TukeyHSD tests between
#' sample groups within each gene. Returns as a data frame.
#'
#' @param sample_info additional sample info from experimental design file
#' @param gene_data data frame of gene expression data
#' @param selected_gene list of genes selected for analysis; Only required if 
#' summarized = FALSE
#' @param summarized TRUE/FALSE whether data has been summarized
#'
#' @returns Data frame containing TukeyHSD results
#' @export
#'
get_tukey_data <- function(gene_data,
                           sample_info,
                           selected_gene = NULL,
                           summarized){
  
  if (summarized == FALSE){
    gene_data <- as.data.frame(gene_data)
    gene_data$symbol <- rownames(gene_data)
    
    gene_data <- gene_data |>
      dplyr::filter(symbol %in% selected_gene) |>
      tidyr::pivot_longer(!symbol, names_to = "sample", values_to = "value")
  }
  
  # Run TukeyHSD test within each gene
  dfTukey <- gene_data |>
    dplyr::mutate(groups = detect_groups(sample, sample_info)) |>
    dplyr::group_by(symbol, groups) |>
    dplyr::mutate(avg = mean(value, na.rm = TRUE)) |>
    dplyr::group_by(symbol) |>
    dplyr::reframe(comp = rownames(TukeyHSD(aov(value ~ groups))$group),
                   diff = TukeyHSD(aov(value ~ groups))$group[,1],
                   lwr = TukeyHSD(aov(value ~ groups))$group[,2],
                   upr = TukeyHSD(aov(value ~ groups))$group[,3],
                   pval = TukeyHSD(aov(value ~ groups))$group[,4],
                   avg = max(avg))
  
  return(dfTukey)
}

#' Data processing message
#'
#' Creates a message about the size of the counts
#' data and the amount of IDs that were converted.
#'
#' @param data_size Data size matrix from \code{\link{pre_process}()}
#' @param all_gene_names Data frame with all gene names
#' @param n_matched Count of matched IDs after processing
#'
#' @export
#' @return Message about processed data
conversion_counts_message <- function(data_size,
                                      all_gene_names,
                                      n_matched) {
  if (ncol(all_gene_names) == 1) {
    return(paste(
      "After mapping, there are",
      data_size[1], "unique gene in", data_size[4], "samples. Of the",
      data_size[3], " genes passed filter. Original gene IDs used."
    ))
  } else {
    return(paste(
      "After mapping, there are",
      data_size[1], "unique genes in", data_size[4], "samples. Of the",
      data_size[3], " genes passed filter, ", n_matched,
      " were converted to Ensembl/STRING gene IDs in our database.
      The remaining ", data_size[3] - n_matched, " genes were
      kept in the data using original IDs."
    ))
  }
}

#' Create message for sequencing depth bias
#'
#' This function creates a warning message for the
#' UI to present to the user regarding the sequencing
#' depth bias.
#'
#' @param raw_counts Matrix of faw counts data from \code{\link{pre_process}()}
#' @param sample_info Matrix of experiment information about each sample
#'
#' @export
#' @return Message for the UI
counts_bias_message <- function(raw_counts,
                                data_file_format,
                                sample_info) {
  total_counts <- colSums(raw_counts)
  groups <- as.factor(
    detect_groups(
      colnames(raw_counts),
      sample_info,
      preserve_original = TRUE
    )
  )
  message <- NULL
      # all samples in one gropu         no grouping
  if(length(unique(groups)) == 1 || length(unique(groups)) == ncol(raw_counts)) {
    return(message)
  }
  # ANOVA of total read counts vs sample groups parsed by sample name
  pval <- summary(aov(total_counts ~ groups))[[1]][["Pr(>F)"]][1]
  means <- aggregate(total_counts, by = list(groups), mean)
  max_min_ratio <- max(means[, 2]) / min(means[, 2])

  if (is.null(pval)) {
    message <- NULL
  } else if (pval < 0.05) {
    message <- paste(
      "Warning! Total read counts are
       significantly different among sample groups
       (p=", sprintf("%-3.2e", pval), ") based on ANOVA.
       Total read counts max/min =", round(max_min_ratio, 2),
       "Some groups were sequenced deeper than others. Proceed analysis with caution."
    )
  }
  # ANOVA of total read counts vs factors in experiment design
  if (!is.null(sample_info) && data_file_format == 1) {
    y <- sample_info
    for (j in 1:ncol(y)) {
      pval <- summary(
        aov(
          total_counts ~ as.factor(y[, j])
        )
      )[[1]][["Pr(>F)"]][1]

      if (pval < 0.01) {
        message <- paste(
          message, " Total read counts seem to be correlated with factor",
          colnames(y)[j], "(p=", sprintf("%-3.2e", pval), ").  "
        )
      }
    }
  }
  return(message)
}

#' Mean vs. Standard Deviation plot
#'
#' Create a plot that shows the standard deviation as the
#' Y-axis across the mean of the counts data on the X-axis.
#' Option to make the X-axis the rank of the mean which
#' does a better job showing the spread of the data.
#'
#' @param processed_data Matrix of data that has gone through
#'  \code{\link{pre_process}()}
#' @param rank TRUE/FALSE whether to use the rank of the mean or not
#' @param heat_cols Heat color to use with black in the plot
#'
#' @export
#' @return A formatted \code{ggplot} hexplot of the mean and standard
#'  deviation of the processed data
#'
#' @family plots
#' @family preprocess functions
mean_sd_plot <- function(processed_data,
                         rank,
                         heat_cols) {
  table_data <- data.frame(
    "x_axis" = apply(
      processed_data,
      1,
      mean
    ),
    "y_axis" = apply(
      processed_data,
      1,
      sd
    )
  )

  if (rank) {
    table_data$x_axis <- rank(table_data$x_axis)
  }
  low_col <- "black"
  high_col <- heat_cols[[1]]

  hex_plot <- ggplot2::ggplot(
    table_data,
    ggplot2::aes(x = x_axis, y = y_axis)
  ) +
    ggplot2::geom_hex(bins = 75) +
    ggplot2::geom_smooth(
      method = "gam",
      formula = y ~ s(x, bs = "cs")
    ) +
    ggplot2::scale_fill_gradient2(
      mid = low_col,
      high = high_col
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = .5
      ),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.title.x = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      axis.title.y = ggplot2::element_text(
        color = "black",
        size = 14
      ),
      legend.text = ggplot2::element_text(size = 12)
    ) +
    ggplot2::labs(
      title = "Mean vs. Standard Deviation",
      y = "Standard Deviation"
    )

  if (rank) {
    hex_plot <- hex_plot +
      ggplot2::labs(x = "Rank of Mean")
  } else {
    hex_plot <- hex_plot +
      ggplot2::labs(x = "Mean")
  }

  return(hex_plot)
}

#' Write paragraph containing process details
#'
#' @param missing_value Method to deal with missing data
#' @param data_file_format Type of data being examined
#' @param low_filter_fpkm Low count filter for the fpkm data
#' @param n_min_samples_fpkm Min samples for fpkm data
#' @param log_transform_fpkm Type of transformation for fpkm data
#' @param log_start_fpkm Value added to log transformation for fpkm
#' @param min_counts Low count filter for count data
#' @param n_min_samples_count Min sample for count data
#' @param counts_transform Type of transformation for counts data
#' @param counts_log_start Value added to log for counts data
#' @param no_fdr Fold changes only data with no p values
#'
#' @return string with process summary
#'
generate_descr <- function(missing_value,
                           data_file_format,
                           low_filter_fpkm,
                           n_min_samples_fpkm,
                           log_transform_fpkm,
                           log_start_fpkm,
                           min_counts,
                           n_min_samples_count,
                           counts_transform,
                           counts_log_start,
                           no_fdr) {
  # read counts case
  if (data_file_format == 1) {
    part_2 <- switch(counts_transform,
      "1" = paste0("EdgeR using a pseudocount of ", counts_log_start),
      "2" = "VST: Variance Stabilizing Transformation",
      "3" = "Regularized log"
    )
    descr <- paste0(
      "Read Counts data were analyzed using iDEP v", packageVersion("idepGolem"), ". ",
      "The data was first filtered to remove reads below ", min_counts,
      " CPM in at least ", n_min_samples_count, " sample", ifelse(n_min_samples_count > 1, "s", ""), ". ",
      "Then the data was transformed with ", part_2, ". ",
      "Missing values were imputed using ", missing_value, "."
    )
  }
  # normalized expression values
  if (data_file_format == 2) {
    part_2 <- switch(toString(log_transform_fpkm),
      "FALSE" = "not log transformed",
      "TRUE" = paste0("log transformed with a psuedocount of ", log_start_fpkm, "")
    )

    descr <- paste0(
      "Normalized expression values were analyzed using iDEP v", packageVersion("idepGolem"), ". ",
      "The data was first filtered to remove genes below ", low_filter_fpkm,
      " in at least ", n_min_samples_fpkm, " sample", ifelse(n_min_samples_fpkm > 1, "s", ""), ". ",
      "The data was ", part_2, ". ",
      "Missing values were imputed using ", missing_value, "."
    )
  }
  # LFC and FDR
  if (data_file_format == 3) {
    descr <- paste0(
      "Log Fold Change and corrected p-value data were analyzed using iDEP v",
      packageVersion("idepGolem"), ". ",
      "Missing values were imputed using ", missing_value, "."
    )
  }
  if (data_file_format == 4) {
    descr <- paste0(
      "Log Fold Change (no p-value) data were analyzed using iDEP v",
      packageVersion("idepGolem"), ". ",
      "Missing values were imputed using ", missing_value, "."
    )
  }

  return(descr)
}


#' Generate colors for plots, solving the problem of not having enough colors for some selections
#'
#' Creates a message about the size of the counts
#' data and the amount of IDs that were converted.
#'
#' @param n Number of desired colors
#' @param palette_name Selected color set, i.e. Set1, Set2, etc.
#'
#' @export
#' @return a list of n colors, even if there are fewer colors in the palette
#' 
# Custom function to generate n colors from a palette
generate_colors <- function(n, palette_name = "Set1") {
  # Ensure RColorBrewer is available
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required but not installed.")
  }
  
  # Validate palette name and get the maximum number of colors available
  if (!palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
    stop("Invalid palette name provided. Please choose a valid palette name from RColorBrewer.")
  }
  
  max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
  
  if (n <= max_colors) {
    # If n is less than or equal to max_colors, use brewer.pal() directly
    return(RColorBrewer::brewer.pal(n, name = palette_name))
  } else {
    # If n is greater than max_colors, interpolate new colors
    base_colors <- RColorBrewer::brewer.pal(max_colors, name = palette_name)
    additional_colors_needed <- n - max_colors
    # Interpolate additional colors using grDevices
    additional_colors <- grDevices::colorRampPalette(base_colors)(n)
    return(additional_colors)
  }
}

