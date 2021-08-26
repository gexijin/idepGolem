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

##### work in process, want to rewrite other code that is relays on first
# plot the expression of one or more genes in the preprocess tab
plot_genes <- function(converted_data, all_gene_info, readSampleInfo, geneSearch, genePlotBox, useSD, selectOrg) {
  symbols <- rownames(converted_data)
  if (selectOrg != "NEW" && ncol(all_gene_info) != 1) {
    ix <- match(symbols, all_gene_info[, 1])
    # symbol really exists?
    if (sum(is.na(all_gene_info$symbol)) != nrow(all_gene_info)) {
      Symbols <- as.character(allGeneInfo$symbol[ix])
      Symbols[which(nchar(Symbols) <= 2)] <- rownames(x)[which(nchar(Symbols) <= 2)]
    }
  }
  x <- as.data.frame(x)
  x$Genes <- Symbols

  # matching from the beginning of symbol
  searchWord <- gsub("^ ", "", geneSearch)
  ix <- which(regexpr(paste("^", toupper(searchWord), sep = ""), toupper(x$Genes)) > 0)
  if (grepl(" $", searchWord)) { # if there is space character at the end, do exact match
    ix <- match(gsub(" ", "", toupper(searchWord)), toupper(x$Genes))
  }

  if (grepl(",|;", searchWord)) { # if there is comma or semicolon, split into multiple words
    Words <- unlist(strsplit(searchWord, ",|;")) # split words
    Words <- gsub(" ", "", Words)
    ix <- match(toupper(Words), toupper(x$Genes))
  }
  ix <- ix[!is.na(ix)] # remove NAs
  # too few or too many genes found
  if (length(ix) == 0 | length(ix) > 50) {
    return(NULL)
  }
  # no genes found

  mdf <- melt(x[ix, ], id.vars = "Genes", value.name = "value", variable.name = "samples")
  # bar plot of individual samples
  p1 <- ggplot(data = mdf, aes(x = samples, y = value, group = Genes, shape = Genes, colour = Genes)) +
    geom_line() +
    geom_point(size = 5, fill = "white") + # shape=21  circle
    # theme(axis.text.x = element_text(size=16,angle = 45, hjust = 1)) +
    labs(y = "Transformed expression level") +
    coord_cartesian(ylim = c(0, max(mdf$value)))
  p1 <- p1 + theme(plot.title = element_text(size = 16, hjust = 0.5)) + # theme(aspect.ratio=1) +
    theme(
      axis.text.x = element_text(angle = 45, size = 16, hjust = 1),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 16)
    ) +
    theme(legend.text = element_text(size = 12))


  # ggplotly(p) %>% layout(margin = list(b = 250,l=100))  # prevent cutoff of sample names

  # Barplot with error bars
  mdf$count <- 1
  g <- detectGroups(mdf$samples, readSampleInfo)
  mdf$g <- g

  options(dplyr.summarise.inform = FALSE)
  # calculate mean, SD, N, per gene per condition
  summarized <- mdf %>%
    group_by(g, Genes) %>%
    summarise(Mean = mean(value), SD = sd(value), N = sum(count))
  colnames(summarized) <- c("Samples", "Genes", "Mean", "SD", "N")
  summarized$SE <- summarized$SD / sqrt(summarized$N)

  if (grepl(",|;", searchWord)) { # re-order according to user input, not alphabetically
    levels <- unique(summarized$Genes)
    iy <- match(toupper(Words), toupper(levels))
    levels <- levels[iy]
    summarized$Genes <- factor(summarized$Genes, levels = levels)
  }

  # http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization
  p2 <- ggplot(summarized, aes(x = Genes, y = Mean, fill = Samples)) + # data & aesthetic mapping
    geom_bar(stat = "identity", position = position_dodge()) + # bars represent average
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, position = position_dodge(.9)) +
    labs(y = "Expression Level")
  if (useSD == 1) {
    p2 <- ggplot(summarized, aes(x = Genes, y = Mean, fill = Samples)) + # data & aesthetic mapping
      geom_bar(stat = "identity", position = position_dodge()) + # bars represent average
      geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(.9)) +
      labs(y = "Expression Level")
  }

  p2 <- p2 + theme(plot.title = element_text(size = 16, hjust = 0.5)) + # theme(aspect.ratio=1) +
    theme(
      axis.text.x = element_text(angle = 45, size = 16, hjust = 1),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 16)
    ) +
    theme(legend.text = element_text(size = 16))

  if (genePlotBox == 1) {
    return(p1)
  } else {
    return(p2)
  }
}

#' @title Pre-Process the data
#' 
#' @description This function takes in user defined values to
#' process the data for the EDA that occurs on the second page.
#' All the filtering and transformation that the app does will
#' occur on this page.
#' 
#' @param data Data that has already gone through the load data
#' 
#' 
#' 
pre_process <- function(
  data,
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
) {

  withProgress(message = "Reading and pre-processing ", {
    data_type_warning <- 0
    data_size_original <- dim(data)
    kurtosis_log <- 50

    # Sort by standard deviation -----------
    data <- data[order(-apply(
      data[, 1:dim(data)[2]],
      1,
      sd
    )), ]

    # Missng values in expression data ----------
    if (sum(is.na(data)) > 0) {
      if (missing_value == "geneMedian") {
        row_medians <- apply(data, 1, function(y) median(y, na.rm = T))
        for (i in 1:dim(data)[2]) {
          val_miss_row <- which(is.na(data[, i]))
          data[val_miss_row, i] <- row_medians[val_miss_row]
        }
      } else if (missing_value == "treatAsZero") {
        data[is.na(data)] <- 0
      } else if (missing_value == "geneMedianInGroup") {
        sample_groups <- detect_groups(colnames(data))
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
              data[missing, i] <- row_medians[misssing]
            }
          }
        }
        if (sum(is.na(data)) > 0) {
          row_medians <- apply(
            data,
            1,
            function(y) median(y, na.rm = T)
          )
          for (i in 1:dim(data)[2]) {
            missing <- which(is.na(data[, i]))
            data[missing, i] <- row_medians[missing]
          }
        }
      }
    }
    # Compute kurtosis ---------
    mean_kurtosis <- mean(apply(data, 2, e1071::kurtosis), na.rm = T)
    raw_counts <- NULL
    pvals <- NULL

    # Pre-processing for each file format ----------
    if (data_file_format == 2) {
      incProgress(1 / 3, "Pre-processing data")
      if (is.integer(data)){
        data_type_warning <- 1
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
        (log_transform_fpkm == TRUE) |
        (mean_kurtosis > kurtosis_log)
      ) {
        data <- log(data + abs(log_start_fpkm), 2)
      }

      std_dev <- apply(data, 1, sd)
      data <- data[order(-std_dev), ]
    } else if (data_file_format == 1) {
      incProgress(1 / 3, "Pre-processing counts data")

      if (!is.integer(data) & mean_kurtosis < kurtosis_log) {
        data_type_warning <- -1
      }

      data <- round(data, 0)

      data <- data[which(apply(
        edgeR::cpm(edgeR::DGEList(counts = data)),
        1,
        function(y) sum(y >= min_counts)
      ) >= n_min_samples_count), ]

      raw_counts <- data 

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

      incProgress(1 / 2, "transforming raw counts")

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
    } else if (data_file_format == 3) {
      n2 <- (dim(data)[2] %/% 2)
      if (!input$no_fdr) {
        pvals <- data[, 2 * (1:n2), drop = FALSE]
        data <- data[, 2 * (1:n2) - 1, drop = FALSE]
        if (dim(data)[2] == 1) {
          placeholder <- rep(1, dim(data)[1])
          pvals <- cbind(pvals, placeholder)
          zero_placeholder <- rep(0, dim(data)[1])
          data <- cbind(data, zero_placeholder)
        }
      }
    }
    data_size <- dim(data)

    validate(
      need(
        dim(data)[1] > 5 & dim(data)[2] >= 1,
        "Data file not recognized. Please double check."
      )
    )

    incProgress(1, "Done.")
  })

  results <- list(
    data = as.matrix(data),
    mean_kurtosis = mean_kurtosis,
    raw_counts = raw_counts,
    data_type_warning = data_type_warning,
    data_size = c(data_size_original, data_size),
    pvals = pvals
  )

  return(results)
}
  

