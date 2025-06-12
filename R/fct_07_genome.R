#' fct_07_genome.R This file holds all of the main data analysis functions
#' associated with eighth tab of the iDEP website.
#'
#'
#' @section fct_07_genome.R functions:
#'
#'
#' @name fct_07_genome.R
NULL

#' Plotly of chromosome position
#'
#' Calculate a plotly that shows the chromosome position of significant genes
#' and significantly enriched regions.
#'
#' @param limma Return from \code{limma_value} function
#' @param select_contrast DEG contrast to examine
#' @param all_gene_info Gene information return from \code{get_gene_info}
#' @param label_gene_symbol Paste the gene symbol label on the plot
#' @param ma_window_size Moving average window size for a chromosome
#'   (1, 2, 4, 6, 8, 10, 15, 20)
#' @param ma_window_steps Number of moving average window steps (1, 2, 3, 4)
#' @param ch_region_p_val P-value to use for finding significant chromosome
#'  region enrichment
#' @param hide_patches Boolean to indicate to only keep within 2 MAD from the
#'  median (TRUE/FALSE)
#' @param hide_chr Boolean to indicate if chromosomes with less than 100 genes
#'  are excluded (TRUE/FALSE)
#' @param x data frame of significant genes and corresponding chromosomes
#' @param x0 additional plotting data -required-
#' @param moving_average data frame of enriched region, moving average data
#' @param ch_length_table table/df of chromosome lengths
#'
#' @export
#' @return Plotly visualization of chromosomes and significantly enriched genes
chromosome_plotly <- function(limma,
                              select_contrast,
                              all_gene_info,
                              label_gene_symbol,
                              ma_window_size,
                              ma_window_steps,
                              ch_region_p_val,
                              hide_patches,
                              hide_chr,
                              x,
                              x0,
                              moving_average,
                              ch_length_table) {
  # Default plot
  fake <- data.frame(a = 1:3, b = 1:3)
  p <- ggplot2::ggplot(fake, ggplot2::aes(x = a, y = b)) +
    ggplot2::geom_blank() +
    ggplot2::ggtitle("No genes with position info.") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )
  if (length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]
  } else {
    top <- limma$top_genes
    ix <- match(select_contrast, names(top))
    if (is.na(ix)) {
      return(plotly::ggplotly(p))
    }
    top_1 <- top[[ix]]
  }
  if (dim(top_1)[1] == 0) {
    return(plotly::ggplotly(p))
  }

  # Species in STRING-db do not have chr. location data
  if (sum(!is.na(all_gene_info$start_position)) < 5) {
    return(p)
  }
  
  # If no chromosomes found. For example if user do not convert gene IDs.
  if (dim(x)[1] > 5) {

    # gather info on chrosomes
    tem <- sort(table(x$chromosome_name), decreasing = TRUE)
    # ch with less than 100 genes are excluded
    if (hide_chr) {
      ch <- names(tem[tem >= 5])
    } else {
      ch <- names(tem[tem >= 1])
    }

    if (length(ch) > 50) {
      # At most 50 ch
      ch <- ch[1:50]
    }

    if (hide_patches) {
      tem <- nchar(ch)
      # only keep within 2 MAD from the median
      ix <- which(tem <= median(tem) + 2 * mad(tem) + 0.0001)
      ch <- ch[ix]
      # ch. name less than 10 characters
      #    ch <- ch[nchar(ch) <= 12]
    }
    ch <- ch[order(as.numeric(ch))]
    tem <- ch
    # The numbers are continous from 1 to length(ch)
    ch <- 1:(length(ch))
    # The names are real chr. names
    names(ch) <- tem

    # Distance between chs.
    chD <- 30
    # Max log2 fold
    fold_cutoff <- 4
    ch_total <- dim(ch_length_table)[1]
    
    # Don't define x and y, so that we could plot use two datasets
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = x,
        ggplot2::aes(
          x = x,
          y = y,
          color = R,
          text = paste0(
            "Symbol: ",
            symbol,
            "\nRegulation: ",
            R,
            "\nChr Pos: ",
            x
          )
        ),
        shape = 20,
        size = 0.2
      )

    if (label_gene_symbol) {
      p <- p + ggplot2::geom_text(
        data = x,
        ggplot2::aes(x = x, y = y, label = symbol),
        check_overlap = FALSE,
        angle = 45,
        size = 2,
        vjust = 0,
        nudge_y = 0.35
      )
    }
    # Label y with ch names
    p <- p + ggplot2::scale_y_continuous(
      labels = paste(
        "chr",
        names(ch[ch_length_table$chNum]),
        sep = ""
      ),
      breaks = chD * (1:ch_total),
      limits = c(0, chD * (ch_total + 1) + 5)
    )

    # Draw horizontal lines for each ch.
    for (i in 1:dim(ch_length_table)[1]) {
      p <- p + ggplot2::annotate(
        "segment",
        x = 0,
        xend = ch_length_table$start_position[i],
        y = ch_length_table$chNum[i] * chD,
        yend = ch_length_table$chNum[i] * chD
      )
    }
    # Change legend  http://ggplot2.tidyverse.org/reference/scale_manual.html
    # Customize legend text
    p <- p + ggplot2::scale_color_manual(
      name = "",
      values = c("red", "blue"),
      breaks = c("Up", "Down"),
      labels = c("Up", "Dn")
    )
    p <- p + ggplot2::xlab("Position on chrs. (Mbp)") +
      ggplot2::theme(axis.title.y = ggplot2::element_blank())
    p <- p + ggplot2::theme(legend.position = "none")

    window_size <- as.numeric(ma_window_size)
    # Step size is then windowSize / steps
    steps <- as.numeric(ma_window_steps)

    # Significant regions are marked as horizontal error bars
    if (dim(moving_average)[1] > 0) {
      p <- p + ggplot2::geom_errorbarh(
        data = moving_average,
        ggplot2::aes(
          x = x,
          y = y,
          xmin = xmin,
          xmax = xmax,
          colour = ma,
          text = paste0(
            paste0(
              "Window: ",
              xmin,
              " - ",
              xmax,
              "\nRegulation: ",
              ma,
              "\nChr Pos: ",
              x
            )
          )
        ),
        size = 2,
        height = 15
      )

      # Label significant regions
      sig_ch <- sort(table(moving_average$chNum), decreasing = TRUE)
      sig_ch <- names(ch)[as.numeric(names(sig_ch))]
      if (length(sig_ch) <= 5) {
        # More than 5 just show 5
        sig_ch <- paste0("chr", sig_ch, collapse = ", ")
      } else {
        sig_ch <- sig_ch[1:5]
        sig_ch <- paste0("chr", sig_ch, collapse = ", ")
        sig_ch <- paste0(sig_ch, ", ...")
      }

      sig_ch <- paste(
        dim(moving_average)[1],
        " enriched regions \n(",
        round(
          sum(ch_length_table$start_position) /
            window_size * steps * as.numeric(ch_region_p_val),
          2
        ),
        " expected)  detected on:\n ",
        sig_ch
      )

      p <- p + ggplot2::annotate(
        geom = "text",
        x = max(x$x) * 0.70,
        y = max(x$y) * 0.90,
        label = sig_ch
      )
    }
  }
  plotly::ggplotly(p, tooltip = "text")
}

#' Chromosome Plot Data
#' 
#' Creates data sets used for gene-chromosome segment plot. Returns four data
#' frames of significant genes, enriched region boundaries, enriched genes, and
#' other necessary plotting data.
#'
#' @param limma Return from \code{limma_value} function
#' @param select_contrast DEG contrast to examine
#' @param all_gene_info Gene information return from \code{get_gene_info}
#' @param ignore_non_coding When TRUE only use protein coding genes
#' @param limma_p_val_viz Adjusted p-value to use for significant genes
#' @param limma_fc_viz Minimum fold-change value to filter with
#' @param ma_window_size Moving average window size for a chromosome
#'   (1, 2, 4, 6, 8, 10, 15, 20)
#' @param ma_window_steps Number of moving average window steps (1, 2, 3, 4)
#' @param ch_region_p_val P-value to use for finding significant chromosome
#'  region enrichment
#' @param hide_patches Boolean to indicate to only keep within 2 MAD from the
#'  median (TRUE/FALSE)
#' @param hide_chr Boolean to indicate if chromosomes with less than 100 genes
#'  are excluded (TRUE/FALSE)
#'
#' @returns Chromosome data used in interactive plot
#' @export
chromosome_data <- function(limma,
                            select_contrast,
                            all_gene_info,
                            ignore_non_coding,
                            limma_p_val_viz,
                            limma_fc_viz,
                            ma_window_size,
                            ma_window_steps,
                            ch_region_p_val,
                            hide_patches,
                            hide_chr){
  
  if (length(limma$comparisons) == 1) {
    top_1 <- limma$top_genes[[1]]
  } else {
    top <- limma$top_genes
    ix <- match(select_contrast, names(top))
    if (is.na(ix)) {
      return(NULL)
    }
    top_1 <- top[[ix]]
  }
  if (dim(top_1)[1] == 0) {
    return(NULL)
  }

  # Species in STRING-db do not have chr. location data
  if (sum(!is.na(all_gene_info$start_position)) < 5) {
    return(NULL)
  }

  colnames(top_1) <- c("Fold", "FDR")

  x <- merge(
    top_1,
    all_gene_info,
    by.x = "row.names",
    by.y = "ensembl_gene_id"
  )

  colnames(x)[which(colnames(x) == "Row.names")] <- "ensembl_gene_id"

  # only coding genes?
  if (ignore_non_coding) {
    x <- subset(x, gene_biotype == "protein_coding")
  }
  
  # If no chromosomes found. For example if user do not convert gene IDs.
  if (dim(x)[1] > 5) {
    x <- x[order(x$chromosome_name, x$start_position), ]
    x$ensembl_gene_id <- as.character(x$ensembl_gene_id)
    # If symbol is missing use Ensembl id
    x$symbol <- as.character(x$symbol)
    ix <- which(is.na(x$symbol))
    ix2 <- which(nchar(as.character(x$symbol)) <= 2)
    ix3 <- which(duplicated(x$symbol))
    ix <- unique(c(ix, ix2, ix3))
    x$symbol[ix] <- x$ensembl_gene_id[ix]

    x <- x[!is.na(x$chromosome_name), ]
    x <- x[!is.na(x$start_position), ]

    # Only keep significant genes
    ix <- which(
      (x$FDR < as.numeric(limma_p_val_viz)) &
        (abs(x$Fold) > log2(as.numeric(limma_fc_viz)))
    )

    if (length(ix) < 5) {
      return(NULL)
    }
    # Remove non-significant / not selected genes
    # Keep a copy
    x0 <- x
    x <- x[ix, ]

    # gather info on chrosomes
    tem <- sort(table(x$chromosome_name), decreasing = TRUE)
    # ch with less than 100 genes are excluded
    if (hide_chr) {
      ch <- names(tem[tem >= 5])
    } else {
      ch <- names(tem[tem >= 1])
    }

    if (length(ch) > 50) {
      # At most 50 ch
      ch <- ch[1:50]
    }

    if (hide_patches) {
      tem <- nchar(ch)
      # only keep within 2 MAD from the median
      ix <- which(tem <= median(tem) + 2 * mad(tem) + 0.0001)
      ch <- ch[ix]
      # ch. name less than 10 characters
      #    ch <- ch[nchar(ch) <= 12]
    }
    ch <- ch[order(as.numeric(ch))]
    tem <- ch
    # The numbers are continous from 1 to length(ch)
    ch <- 1:(length(ch))
    # The names are real chr. names
    names(ch) <- tem
    x <- x[which(x$chromosome_name %in% names(ch)), ]
    x <- droplevels(x)
    # Numeric encoding
    x$chNum <- 1
    x$chNum <- ch[x$chromosome_name]

    # Use max position as chr. length before filtering
    ch_length_table <- aggregate(start_position ~ chromosome_name, data = x, max)
    # Add chr. numer
    ch_length_table$chNum <- ch[ch_length_table$chromosome_name]
    ch_length_table <- ch_length_table[!is.na(ch_length_table$chNum), ]
    ch_length_table <- ch_length_table[order(ch_length_table$chNum), c(3, 2)]
    ch_length_table <- ch_length_table[order(ch_length_table$chNum), ]
    ch_length_table$start_position <- ch_length_table$start_position / 1e6

    # Prepare coordinates
    # Mbp
    x$start_position <- x$start_position / 1000000
    # Distance between chs.
    chD <- 30
    # Max log2 fold
    fold_cutoff <- 4

    # log2 fold within -5 to 5
    x$Fold[which(x$Fold > fold_cutoff)] <- fold_cutoff
    x$Fold[which(x$Fold < -1 * fold_cutoff)] <- -1 * fold_cutoff
    x$Fold <- 4 * x$Fold

    x$y <- x$chNum * chD + x$Fold
    ch_total <- dim(ch_length_table)[1]
    x$R <- as.factor(sign(x$Fold))

    colnames(x)[which(colnames(x) == "start_position")] <- "x"
    x$R <- as.character(x$R)
    x$R[x$R == "-1"] <- "Down"
    x$R[x$R == "1"] <- "Up"
    x$R <- as.factor(x$R)
    
    # remove chromosomes with no genes left
    x <- droplevels(x)
    
    x0 <- x0[x0$chromosome_name %in% unique(x$chromosome_name), ]
    # Numeric encoding
    x0$chNum <- 1
    x0$chNum <- ch[x0$chromosome_name]
    # Mbp
    x0$start_position <- x0$start_position / 1e6

    window_size <- as.numeric(ma_window_size)
    # Step size is then windowSize / steps
    steps <- as.numeric(ma_window_steps)

    # Centering
    x0$Fold <- x0$Fold - mean(x0$Fold)

    for (i in 0:(steps - 1)) {
      # Step size is windowSize/steps
      # If windowSize=10 and steps = 2; then step size is 5Mb
      # 1.3 becomes 5, 11.2 -> 15 for step 1
      # 1.3 -> -5
      x0$x <- (
        floor((x0$start_position - i * window_size / steps) / window_size)
        + 0.5 + i / steps
      ) * window_size

      moving_average_start <- x0 |>
        dplyr::select(chNum, x, Fold) |>
        # Beginning bin can be negative for first bin in the 2nd step
        dplyr::filter(x >= 0) |>
        dplyr::group_by(chNum, x) |>
        dplyr::summarize(
          ma = mean(Fold),
          n = dplyr::n(),
          pval = ifelse(
            dplyr::n() >= 3 && sd(Fold) > 0, stats::t.test(Fold)$p.value, 0
          )
        ) |>
        # na when only 1 data point?
        dplyr::filter(!is.na(pval))

      if (i == 0) {
        moving_average <- moving_average_start
      } else {
        moving_average <- rbind(moving_average, moving_average_start)
      }
    }
    # Translate fold to y coordinates
    moving_average <- moving_average |>
      dplyr::filter(n >= 3) |>
      dplyr::mutate(pval = stats::p.adjust(pval, method = "fdr")) |>
      dplyr::filter(pval < as.numeric(ch_region_p_val)) |>
      dplyr::mutate(y = ifelse(ma > 0, 1, -1),
                    ma = ifelse(ma > 0, 1, -1)) |>
      dplyr::mutate(y = chNum * chD + 3 * y,
                    ma = as.factor(ma),
                    xmin = x - window_size / 2,
                    xmax = x + window_size / 2)

    moving_average$ma <- as.character(moving_average$ma)
    moving_average$ma[moving_average$ma == "-1"] <- "Down"
    moving_average$ma[moving_average$ma == "1"] <- "Up"
    moving_average$ma <- as.factor(moving_average$ma)
    
    # Enriched intervals
    gr_enriched <- GenomicRanges::GRanges(
      seqnames = moving_average$chNum,
      ranges = IRanges::IRanges(start = moving_average$xmin,
                                end = moving_average$xmax)
    )
    
    # Raw genes as 1-base pair intervals
    gr_raw <- GenomicRanges::GRanges(
      seqnames = x$chNum,
      ranges = IRanges::IRanges(start = x$x, end = x$x),
      gene_id = x$ensembl_gene_id
    )
    
    # Overlap: matches any gene position falling in any enriched region
    hits <- GenomicRanges::findOverlaps(gr_raw, gr_enriched)
    
    # Filter raw gene data
    filtered_genes <- x[S4Vectors::queryHits(hits), ]
    filtered_genes <- filtered_genes[!duplicated(filtered_genes), ]
    
    return(list("chr_data" = x,
                "enriched_regions" = moving_average,
                "other" = x0,
                "enriched_genes" = filtered_genes,
                "ch_length" = ch_length_table))
  }
}

#' Filter chromosome data
#'
#' Filters gene chromosome data by regulation and list of chromosomes selected.
#' "All" selected for either negates filtering. Returns data frame of filtered
#' chromosomal data.
#' 
#' @param chr_data chromosome data
#' @param regulation Selected gene regulation
#' @param chr_select list of selected chromosomes
#'
#' @returns Filtered chromosome data
#' @export
#'
chr_filter <- function(chr_data,
                       regulation,
                       chr_select){
  
  if (regulation != "All"){
    chr_data <- dplyr::filter(chr_data, R == regulation)
  }
  
  if (!("All" %in% chr_select)){
    chr_data <- dplyr::filter(chr_data, chromosome_name %in% chr_select)
  }
  
  return(chr_data)
}