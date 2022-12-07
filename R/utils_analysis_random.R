#' utils_analysis_random.R Miscellaneous helper functions here,
#' we find a place for them later.These are functions don't manipulate data,
#' but there were used in data analysis functions for other purposes.
#'
#'
#' @section utils_analysis_random.R functions:
#' add later....
#'
#'
#' @name utils_analysis_random.R
NULL

#' List of method functions for hierarchical clustering
#'
#'
#' @description Returns a list of functions for hierarchical clustering
#' including average linkage, ward.d, ward.d2, single linkage, mcquitty,
#' median, and centroid.
#'
#' @export
#' @return  A list of functions
#'
#'
hcluster_functions <- function() {
  # Functions for hierarchical clustering
  hclust_average <- function(x, method = "average", ...) { # average linkage
    hclust(x, method = method, ...)
  }
  hclust_ward_D <- function(x, method = "ward.D", ...) { # ward.D linkage # nolint
    hclust(x, method = method, ...)
  }
  hclust_ward_D2 <- function(x, method = "ward.D2", ...) { # ward.D2 linkage # nolint
    hclust(x, method = method, ...)
  }
  hclust_single <- function(x, method = "single", ...) { # single linkage
    hclust(x, method = method, ...)
  }
  hclust_mcquitty <- function(x, method = "mcquitty", ...) { # mcquitty linkage # nolint
    hclust(x, method = method, ...)
  }
  hclust_median <- function(x, method = "median", ...) { # median linkage
    hclust(x, method = method, ...)
  }
  hclust_centroid <- function(x, method = "centroid", ...) { # centroid linkage # nolint
    hclust(x, method = method, ...)
  }

  return(list(
    complete = hclust,
    average = hclust_average,
    ward_D = hclust_ward_D,
    ward_D2 = hclust_ward_D2,
    single = hclust_single,
    mcquitty = hclust_mcquitty,
    median = hclust_median,
    centroid = hclust_centroid
  ))
}


#' List of distance functions for hierarchical clustering
#'
#'
#' @description Distance functions for clustering
#'
#' @export
#' @return A list of distance functions for hierarchical clustering including
#'  pearson, euclidian, and absolute pearson
#'
#'
dist_functions <- function() {
  dist_pcc <- function(x, ...) {
    # distance function = 1-PCC (Pearson's correlation coefficient)
    as.dist(1 - cor(t(x), method = "pearson"))
  }

  dist_abs_pcc <- function(x, ...) {
    # distance function = 1-abs(PCC) (Pearson's correlation coefficient)
    as.dist(1 - abs(cor(t(x), method = "pearson")))
  }

  return(list(
    Pearson = dist_pcc,
    Euclidean = dist,
    Absolute_Pearson = dist_abs_pcc
  ))
}


#' Dynamic Range
#'
#' Find the difference between 2nd largest and 2nd smallest number in set
#'
#' @param num_set Vector of values
#'
#' @return A numeric value indicating the difference
#'
#' @examples
#' dynamic_range(1:10)
#'
dynamic_range <- function(num_set) {
  # Given a set of numbers,
  # find the difference between 2nd largest and 2nd smallest
  sorted_set <- sort(num_set)
  if (length(num_set) >= 4) {
    k <- 2
  } else {
    k <- 1
  }
  return(sorted_set[length(num_set) - k + 1] - sorted_set[k])
}


#' Detect groups by sample names
#'
#' Detects groups from column names in sample info file so that they can be used
#' for things such as coloring plots or building the model for DEG analysis.
#'
#' @param sample_names Vector of column headings from data file or design file
#' @param sample_info Matrix of the experiment design information
#'
#' @export
#' @return A character vector with the groups
#' @note This function is mainly called internally in other idepGolem functions.
#'
detect_groups <- function(sample_names, sample_info = NULL) {
  # sample_names are col names parsing samples by either the name
  # or using a data frame of sample infos.
  # Note that each row of the sample_info data frame represents a sample.
  sample_group <- NULL
  if (is.null(sample_info)) {
    # Remove all numbers from end
    # remove "_" from end
    # remove "_Rep" from end
    # remove "_rep" from end
    # remove "_REP" from end
    sample_group <- gsub(
      "[0-9]*$", "",
      sample_names
    )
    sample_group <- gsub("_$", "", sample_group)
    sample_group <- gsub("_Rep$", "", sample_group)
    sample_group <- gsub("_rep$", "", sample_group)
    sample_group <- gsub("_REP$", "", sample_group)
  } else {
    # the orders of samples might not be the same.
    # The total number of samples might also differ
    match_sample <- match(sample_names, row.names(sample_info))
    sample_info2 <- sample_info[match_sample, , drop = FALSE]
    if (ncol(sample_info2) == 1) {
      # if there's only one factor
      sample_group <- sample_info2[, 1]
    } else {
      # multiple columns/factors
      foo <- function(y) paste(y, collapse = "_")
      sample_group <- unlist(apply(
        X = sample_info2,
        MARGIN = 1,
        FUN = foo
      ))
      names(sample_group) <- row.names(sample_info2)
      if (min(table(sample_group)) == 1) { # no replicates?
        sample_group <- sample_info2[, 1]
      }
    }
  }
  return(as.character(sample_group))
}


#' Data cleaning on gene sets
#'
#' Clean up a set of genes by removing duplicates, spaces, and other control characters from
#' gene names
#'
#' @param gene_set Vector of genes to be cleans.
#'
#' @export
#' @return Vector of cleaned genes.
#'
#' @note This function is called internally in other idepGolem functions.
clean_gene_set <- function(gene_set) {
  # remove duplicate; upper case; remove special characters
  gene_set <- unique(toupper(gsub("\n| ", "", gene_set)))
  # genes should have at least two characters
  gene_set <- gene_set[which(nchar(gene_set) > 1)]
  return(gene_set)
}


#' Clean gene query
#'
#' Clean and prepare a gene query by adding tabs and new lines after each gene
#'
#' @param query_input Vector with list of genes
#'
#' @export
#' @return Vector of cleaned gene list
#'
#' @note This function is called internally in other idepGolem functions.
clean_query <- function(query_input) {
  return(clean_gene_set(unlist(strsplit(
    x = toupper(query_input),
    split = "\t| |\n|\\,"
  ))))
}


#' Read .gmt file
#'
#' This functions cleans and converts gene names to upper case in  a .gmt file
#'
#' @param file_path String with path to file
#'
#' @return Data frame with data from gmt file
#' @export
#'
read_gmt <- function(file_path) { # size restriction
  # Read in the first file
  gmt_data <- scan(file = file_path, what = "", sep = "\n")
  gmt_data <- gsub(
    pattern = " ",
    replacement = "",
    x = gmt_data
  )
  gmt_data <- toupper(gmt_data)

  #----Process the first file
  # Separate elements by one or more whitespace
  gmt_data <- strsplit(x = gmt_data, split = "\t")
  # Extract the first vector element and set it as the list element name
  extract <- function(x) x[[1]]
  names(gmt_data) <- sapply(X = gmt_data, FUN = extract)
  # Remove the first vector element from each list element
  extract2 <- function(x) x[-1]
  gmt_data <- lapply(X = gmt_data, FUN = extract2)
  # remove duplicated elements
  ## need to try to get this to lappy later
  for (i in 1:length(gmt_data)) {
    gmt_data[[i]] <- clean_gene_set(gmt_data[[i]])
  }
  # check the distribution of the size of gene lists sapply(y, length) hold a vector of sizes
  gene_size <- max(sapply(X = gmt_data, FUN = length))
  if (gene_size < 5) {
    cat("Warning! Gene sets have very small number of genes!\n
     Please double check format.")
  }
  # gene sets smaller than 1 is ignored!!!
  gmt_data <- gmt_data[which(sapply(X = gmt_data, FUN = length) > 1)]
  return(gmt_data)
}


#' Capitalize first letter
#'
#' This function converts gene set names to a more readable format
#'
#' @param x A string
#'
#' @return Cleaned string
#' @examples
#' proper("GOBP_mmu_mgi_GO:0000183_chromatin_silencing_at_rDNA")
proper <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))


#' Extract a word
#'
#' @param word_list word list
#'
#' @return extracted words
#' @examples
#' # ADD_EXAMPLES_HERE
extract_word <- function(word_list) {
  words <- unlist(strsplit(word_list, "_"))
  if (length(words) <= 4) {
    return(gsub("_", " ", word_list))
  } else {
    words <- words[-c(1:4)]
    return(proper(paste(words, collapse = " ")))
  }
}


#' Check object state
#'
#' check_object_state an Utility function to simplify checking object states.
#'
#'
#' A utility function to evaluates the state of objects, and sends messages from
#' server to UI of a shiny app.
#'
#'
#' @param check_exp An expression that should be evaluated
#' @param true_message Message displayed on console and
#'  returned if \code{check_exp} is true
#' @param false_message Optional message returned if \code{check_exp} is false,
#'  default to blank string
#'
#' @return A list is returned, with the following elements:
#'  \code{bool} is either true or false depending on \code{check_exp} evaluation
#'  \code{content} either message and depended on \code{check_exp} evaluation
#'
#' @export
#' @examples
#' check <- check_object_state(
#'   check_exp = (is.null(NULL)),
#'   true_message = as.data.frame("This is NULL")
#' )
#' # will eval to true and display message to console
#'
#' check <- check_object_state(
#'   check_exp = (length(c(1, 2)) == 0),
#'   true_message = "this has 0 elements",
#'   false_message = "this doesn't have 0 elements"
#' )
#' # This will not return check and check$content will be false_message
#'
check_object_state <- function(check_exp,
                               true_message,
                               false_message = "") {
  if (check_exp) {
    message(true_message)
    return(list(
      bool = TRUE,
      content = true_message
    ))
  } else {
    return(list(
      bool = FALSE,
      content = false_message
    ))
  }
}

#' \code{ggplot} colors
#'
#' Create a vector of hex colors to be used in \code{ggplot}s
#'
#' @param n Number of colors to return
#'
#' @return Vector of hex color codes for a plot.
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Swap rowname IDs for data matrix
#'
#' This function uses the all_gene_names dataframe
#' to swap the current rownames of a data matrix with
#' the desired gene ID. For instance, if the rownames
#' were currently ensembl, this function is able to
#' switch them back to the original form.
#'
#' @param data_matrix Data matrix with ensembl or user gene ID rownames
#' @param all_gene_names Data frame of gene names
#' @param select_gene_id String designating the desired ID type for rownames
#'   ("User_ID", "ensembl_ID", "symbol")
#'
#' @export
#' @return Data matrix with changed rownames
rowname_id_swap <- function(data_matrix,
                            all_gene_names,
                            select_gene_id) {
  if (select_gene_id == "User_ID" && ncol(all_gene_names) == 1) {
    return(data_matrix)
  } else if (select_gene_id == "User_ID") {
    data_matrix <- as.data.frame(data_matrix)
    data_matrix$order <- seq(1, nrow(data_matrix), 1)
    new_data <- merge(
      data_matrix,
      all_gene_names,
      by.x = "row.names",
      by.y = "ensembl_ID",
      all.x = T
    )
    rownames(new_data) <- new_data$User_ID
    nums <- unlist(lapply(new_data, is.numeric))
    new_data <- new_data[, nums]
    new_data <- new_data[order(new_data$order), ]
    new_data <- dplyr::select(new_data, -order)
    new_data <- as.matrix(new_data)
    return(new_data)
  } else if (select_gene_id == "ensembl_ID") {
    data_matrix <- as.data.frame(data_matrix)
    data_matrix$order <- seq(1, nrow(data_matrix), 1)
    new_data <- merge(
      data_matrix,
      all_gene_names,
      by.x = "row.names",
      by.y = "ensembl_ID",
      all.x = T
    )
    rownames(new_data) <- new_data$Row.names
    nums <- unlist(lapply(new_data, is.numeric))
    new_data <- new_data[, nums]
    new_data <- new_data[order(new_data$order), ]
    new_data <- dplyr::select(new_data, -order)
    new_data <- as.matrix(new_data)
    return(new_data)
  } else if (select_gene_id == "symbol") {
    data_matrix <- as.data.frame(data_matrix)
    data_matrix$order <- seq(1, nrow(data_matrix), 1)
    new_data <- merge(
      data_matrix,
      all_gene_names,
      by.x = "row.names",
      by.y = "ensembl_ID",
      all.x = TRUE
    )
    rownames(new_data) <- new_data$symbol
    nums <- unlist(lapply(new_data, is.numeric))
    new_data <- new_data[, nums]
    new_data <- new_data[order(new_data$order), ]
    new_data <- dplyr::select(new_data, -order)
    new_data <- as.matrix(new_data)
    return(new_data)
  }
}

#' Merge data from gene info and processed data
#'
#' This function takes in the gene info data and merges it
#' with the data that has gone through the processing function.
#' The returned data contains the gene names as well as the
#' ensembl id in the first two columns.
#'
#' @param all_gene_names All matched gene names from idep data
#' @param data Data matrix with rownames to merge with gene names
#' @param merge_ID String designating the type gene id such as "User_ID",
#'  "ensembl_ID", "symbol"
#'
#'
#' @export
#' @return Data frame with all gene name information.
merge_data <- function(all_gene_names,
                       data,
                       merge_ID) {
  isolate({
    if (dim(all_gene_names)[2] == 1) {
      new_data <- data
      new_data <- as.data.frame(new_data)
      new_data$User_id <- rownames(new_data)
      new_data <- dplyr::select(
        new_data,
        User_id,
        tidyselect::everything()
      )
      rownames(new_data) <- seq(1, nrow(new_data), 1)
      tmp <- apply(new_data[, 2:dim(new_data)[2]], 1, sd)
      new_data <- new_data[order(-tmp), ]

      return(new_data)
    } else if (dim(all_gene_names)[2] == 2) {
      new_data <- merge(
        all_gene_names,
        data,
        by.x = merge_ID,
        by.y = "row.names",
        all.y = T
      )
      new_data <- dplyr::select(
        new_data,
        User_ID,
        ensembl_ID,
        tidyselect::everything()
      )
      rownames(new_data) <- seq(1, nrow(new_data), 1)
      tmp <- apply(new_data[, 3:dim(new_data)[2]], 1, sd)
      new_data <- new_data[order(-tmp), ]

      return(new_data)
    } else {
      new_data <- merge(
        all_gene_names,
        data,
        by.x = merge_ID,
        by.y = "row.names",
        all.y = T
      )
      new_data <- dplyr::select(
        new_data,
        User_ID,
        ensembl_ID,
        symbol,
        tidyselect::everything()
      )
      rownames(new_data) <- seq(1, nrow(new_data), 1)
      tmp <- apply(new_data[, 4:dim(new_data)[2]], 1, sd)
      new_data <- new_data[order(-tmp), ]

      return(new_data)
    }
  })
}


#' Add a legend
#'
#' Add a sample legends to the main heatmap
#' @note https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
add_legend <- function(...) {
  opar <- par(
    fig = c(0, 1, 0, 1),
    oma = c(0, 0, 0, 0),
    mar = c(0, 0, 0, 6),
    new = TRUE
  )
  on.exit(par(opar))
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend(...)
}

#' Wrapping long text by adding a new line character
#'
#' https://stackoverflow.com/questions/7367138/text-wrap-for-plot-titles
wrap_strings <- function(vector_of_strings,
                         width = 30) {
  as.character(
    sapply(vector_of_strings, FUN = function(x) {
      paste(strwrap(x, width = width), collapse = "\n")
    })
  )
}


#' Replace underscore
#' Find any instances of "_" and replace them with a space
#' @param x String
#'
extract_under <- function(x) {
  words <- unlist(strsplit(x, "_"))
  if (length(words) <= 4) {
    return(gsub("_", " ", x))
  } else {
    words <- words[-c(1:4)]
    return(loose_rock_proper(paste(words, collapse = " ")))
  }
}

#' Capitalizes the first letter in all words of a string
#'
#' from: https://github.com/averissimo/loose.rock/blob/master/R/string.R
#'
#' @param x String
#'
#' @return String with capitalized letters
#'
#' #examples
#' loose_rock_proper("i saw a ParRot")
loose_rock_proper <- function(x) {
  return(gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(x), perl = TRUE))
}

#' Pathview function
#' https://rdrr.io/bioc/pathview/src/R/colorpanel2.R
#'
#' @param n Number of HEX color values to return
#' @param low String designating the low color
#' @param mid String designating the mid color
#' @param high String designating the high color
#'
#' @return A list the length of \code{n} with HEX color values.
colorpanel2 <- function(n,
                        low,
                        mid,
                        high) {
  if (missing(mid) || missing(high)) {
    low <- grDevices::col2rgb(low)
    if (missing(high)) {
      high <- grDevices::col2rgb(mid)
    } else {
      high <- grDevices::col2rgb(high)
    }
    red <- seq(low[1, 1], high[1, 1], length = n) / 255
    green <- seq(low[3, 1], high[3, 1], length = n) / 255
    blue <- seq(low[2, 1], high[2, 1], length = n) / 255
  } else {
    isodd <- n %% 2 == 1
    if (isodd) {
      n <- n + 1
    }
    low <- grDevices::col2rgb(low)
    mid <- grDevices::col2rgb(mid)
    high <- grDevices::col2rgb(high)
    lower <- floor(n / 2)
    upper <- n - lower
    red <- c(
      seq(low[1, 1], mid[1, 1], length = lower),
      seq(mid[1, 1], high[1, 1], length = upper)
    ) / 255
    green <- c(
      seq(low[3, 1], mid[3, 1], length = lower),
      seq(mid[3, 1], high[3, 1], length = upper)
    ) / 255
    blue <- c(
      seq(low[2, 1], mid[2, 1], length = lower),
      seq(mid[2, 1], high[2, 1], length = upper)
    ) / 255
    if (isodd) {
      red <- red[-(lower + 1)]
      green <- green[-(lower + 1)]
      blue <- blue[-(lower + 1)]
    }
  }
  grDevices::rgb(red, blue, green)
}

#' PATHVIEW SOURCE FUNCTION
#' https://rdrr.io/bioc/pathview/src/R/render.kegg.node.R
#'
#' @param plot.data Plot data
#' @param cols.ts Columns
#' @param img Image
#' @param same.layer TRUE/FALSE
#' @param type Type
#' @param text.col String designating text color
#' @param cex Size of text
#'
render.kegg.node <- function(plot.data,
                             cols.ts,
                             img,
                             same.layer = TRUE,
                             type = c("gene", "compound")[1],
                             text.col = "black",
                             cex = 0.25) {
  width <- ncol(img)
  height <- nrow(img)
  nn <- nrow(plot.data)
  pwids <- plot.data$width
  if (!all(pwids == max(pwids))) {
    message("Info: ", "some node width is different from others, and hence adjusted!")
    wc <- table(pwids)
    pwids <- plot.data$width <- as.numeric(names(wc)[which.max(wc)])
  }

  if (type == "gene") {
    if (same.layer != T) {
      rect.out <- sliced.shapes(plot.data$x + 0.5, height - plot.data$y, plot.data$width / 2 - 0.5, plot.data$height / 2 - 0.25, cols = cols.ts, draw.border = F, shape = "rectangle")
      text(plot.data$x + 0.5, height - plot.data$y,
        labels = as.character(plot.data$labels),
        cex = cex, col = text.col
      )
      return(invisible(1))
    } else {
      img2 <- img
      pidx <- cbind(
        ceiling(plot.data$x - plot.data$width / 2) + 1,
        floor(plot.data$x + plot.data$width / 2) + 1,
        ceiling(plot.data$y - plot.data$height / 2) + 1,
        floor(plot.data$y + plot.data$height / 2) + 1
      )
      cols.ts <- cbind(cols.ts)
      ns <- ncol(cols.ts)
      brk.x <- sapply(plot.data$width / 2, function(wi) seq(-wi, wi, length = ns + 1))
      for (k in 1:ns) {
        col.rgb <- col2rgb(cols.ts[, k]) / 255
        pxr <- t(apply(pidx[, 1:2], 1, function(x) x[1]:x[2])) - plot.data$x - 1
        sel <- pxr >= ceiling(brk.x[k, ]) & pxr <= floor(brk.x[k + 1, ])
        for (i in 1:nn) {
          sel.px <- (pidx[i, 1]:pidx[i, 2])[sel[i, ]]
          node.rgb <- img[pidx[i, 3]:pidx[i, 4], sel.px, 1:3]
          node.rgb.sum <- apply(node.rgb, c(1, 2), sum)
          blk.ind <- which(node.rgb.sum == 0 | node.rgb.sum == 1, arr.ind = T)
          node.rgb <- array(col.rgb[, i], dim(node.rgb)[3:1])
          node.rgb <- aperm(node.rgb, 3:1)
          for (j in 1:3) node.rgb[cbind(blk.ind, j)] <- 0
          img2[pidx[i, 3]:pidx[i, 4], sel.px, 1:3] <- node.rgb
        }
      }
      return(img2)
    }
  } else if (type == "compound") {
    if (same.layer != T) {
      nc.cols <- ncol(cbind(cols.ts))
      if (nc.cols > 2) { # block the background circle
        na.cols <- rep("#FFFFFF", nrow(plot.data))
        cir.out <- sliced.shapes(plot.data$x, height - plot.data$y, plot.data$width[1], plot.data$width[1], cols = na.cols, draw.border = F, shape = "ellipse", lwd = 0.2)
      }
      cir.out <- sliced.shapes(plot.data$x, height - plot.data$y, plot.data$width[1], plot.data$width[1], cols = cols.ts, shape = "ellipse", blwd = 0.2)
      return(invisible(1))
    } else {
      #    col.rgb=col2rgb(cols.ts)/255
      blk <- c(0, 0, 0)
      img2 <- img
      w <- ncol(img) # repeat
      h <- nrow(img) # repeat
      cidx <- rep(1:w, each = h)
      ridx <- rep(1:h, w)
      pidx <- lapply(1:nn, function(i) {
        ii <- which(
          (cidx - plot.data$x[i])^2 + (ridx - plot.data$y[i])^2 < (plot.data$width[i])^2
        )
        imat <- cbind(cbind(ridx, cidx)[rep(ii, each = 3), ], 1:3)
        imat[, 1:2] <- imat[, 1:2] + 1
        ib <- which(
          abs((cidx - plot.data$x[i])^2 + (ridx - plot.data$y[i])^2 - (plot.data$width[i])^2) <= 8
        )
        ibmat <- cbind(cbind(ridx, cidx)[rep(ib, each = 3), ], 1:3)
        ibmat[, 1:2] <- ibmat[, 1:2] + 1
        return(list(fill = imat, border = ibmat))
      })

      cols.ts <- cbind(cols.ts)
      ns <- ncol(cols.ts)
      brk.x <- sapply(plot.data$width, function(wi) seq(-wi, wi, length = ns + 1))
      for (i in 1:nn) {
        pxr <- pidx[[i]]$fill[, 2] - 1 - plot.data$x[i]
        col.rgb <- col2rgb(cols.ts[i, ]) / 255
        for (k in 1:ns) {
          sel <- pxr >= brk.x[k, i] & pxr <= brk.x[k + 1, i]
          img2[pidx[[i]]$fill[sel, ]] <- col.rgb[, k]
        }
        img2[pidx[[i]]$border] <- blk
      }
      return(img2)
    }
  } else {
    stop("unrecognized node type!")
  }
}

#' PATHVIEW SOURCE FUNCTION
#'
pathview.stamp <- function(x = NULL,
                           y = NULL,
                           position = "bottomright",
                           graph.sizes,
                           on.kegg = TRUE,
                           cex = 1) {
  if (on.kegg) {
    labels <- "Data on KEGG graph\nRendered by Pathview"
  } else {
    labels <- "-Data with KEGG pathway-\n-Rendered  by  Pathview-"
  }
  if (is.null(x) | is.null(y)) {
    x <- graph.sizes[1] * .80
    y <- graph.sizes[2] / 40
    if (length(grep("left", position)) == 1) x <- graph.sizes[1] / 40
    if (length(grep("top", position)) == 1) y <- graph.sizes[2] - y
  }
  text(x = x, y = y, labels = labels, adj = 0, cex = cex, font = 2)
}

#' Heatmap of the data
#'
#' Create a ComplexHeatmap from a data matrix.
#'
#' @param data Matrix of data
#' @param heatmap_color_select Vector of colors to use for the fill
#'  in the heatmap
#' @export
#' @return A drawn ComplexHeatmap.
#'
basic_heatmap <- function(data,
                          heatmap_color_select) {
  # Number of genes to show
  n_genes <- nrow(data)

  data <- as.matrix(data) - apply(data, 1, mean)
  cutoff <- median(unlist(data)) + 3 * sd(unlist(data))
  data[data > cutoff] <- cutoff
  cutoff <- median(unlist(data)) - 3 * sd(unlist(data))
  data[data < cutoff] <- cutoff

  data <- data[which(apply(data, 1, sd) > 0), ]

  # Color scale
  if (min(data) < 0) {
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
    clustering_distance_rows = function(x) {
      as.dist(
        1 - cor(t(x), method = "pearson")
      )
    },
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

#' Heatmap of User brush selection
#'
#' Create a heatmap from the brush selection of the main heatmap.
#' Used in iDEP to create an interactive heatmap and enable the
#' User to zoom in on areas they find interesting.
#'
#' @param ht_brush Brush information from the User on the main
#'  heatmap
#' @param ht Main heatmap to create the sub-heatmap from
#' @param ht_pos_main Position information from the main heatmap
#'  to use for the sub-heatmap
#' @param heatmap_data Data matrix that is being plotted in the
#'  main heatmap
#'
#' @export
#' @return A ComplexHeatmap object that will be inputted into the
#'  draw function in the server, the sub-heatmap data matrix, the
#'  group color mapping for the annotation, and the groups that
#'  the columns fall into.
basic_heat_sub <- function(ht_brush,
                           ht,
                           ht_pos_main,
                           heatmap_data) {
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


#' HTML code for sub-heatmap selected cell
#'
#' Create HTML code for a cell of information on the cell of the
#' sub-heatmap that the User clicks on. The cell contains the
#' expression value, the sample, the gene, and the group.
#'
#' @param click Information from the User clicking on a cell of
#'  the sub-heatmap
#' @param ht_sub The drawn sub-heatmap
#' @param ht_sub_obj The sub-heatmap ComplexHeatmap object
#' @param ht_pos_sub Position information for the sub-heatmap
#' @param sub_groups Vector of the groups that the samples
#'  belong to
#' @param group_colors color of the top annotation that
#'  is used for each group
#' @param data Sub data matrix that is plotted in the sub-heatmap
#'
#' @return HTML code that will be used in the shiny UI to tell
#'  the user the information of the cell they selected.
heat_click_info <- function(click,
                            ht_sub,
                            ht_sub_obj,
                            ht_pos_sub,
                            sub_groups,
                            group_colors,
                            data) {
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
</pre></div>")

  return(HTML(html))
}

#' Convert regular text to hypertext with links
#'
#'
#' @param click Information from the User clicking on a cell of
#'  the sub-heatmap
#' @param textVector A vector of characters
#' @param urlVector a set of URLs
#'
#' @export
#' @return HTML format hyper text
hyperText <- function(textVector, urlVector) {
  # for generating pathway lists that can be clicked.
  # Function that takes a vector of strings and a vector of URLs
  # and generate hyper text
  # add URL to Description
  # see https://stackoverflow.com/questions/30901027/convert-a-column-of-text-urls-into-active-hyperlinks-in-shiny
  # see https://stackoverflow.com/questions/21909826/r-shiny-open-the-urls-from-rendertable-in-a-new-tab
  if (sum(is.null(urlVector)) == length(urlVector)) {
    return(textVector)
  }

  if (length(textVector) != length(urlVector)) {
    return(textVector)
  }

  #------------------URL correction
  # URL changed from http://amigo.geneontology.org/cgi-bin/amigo/term_details?term=GO:0000077
  #                  http://amigo.geneontology.org/amigo/term/GO:0000077
  urlVector <- gsub(
    "cgi-bin/amigo/term_details\\?term=",
    "amigo/term/",
    urlVector
  )
  urlVector <- gsub(" ", "", urlVector)


  # first see if URL is contained in memo
  ix <- grepl("http:", urlVector, ignore.case = TRUE)
  if (sum(ix) > 0) { # at least one has http?
    tem <- paste0(
      "<a href='",
      urlVector, "' target='_blank'>",
      textVector,
      "</a>"
    )
    # only change the ones with URL
    textVector[ix] <- tem[ix]
  }
  return(textVector)
}


#' Change ggplot2 plots
#'
#'
#' @param p \code{ggplot} object
#' @param gridline TRUE/FALSE to indicate adding gridlines
#' @param ggplot2_theme String to indicate which \code{ggplot} theme to use.
#'  Should be one of "linedraw", "classic", "gray", "light", "dark", or "bw"
#'
#' @export
#' @return Formatted \code{ggplot} object
refine_ggplot2 <- function(p, gridline, ggplot2_theme = "light") {
  # apply theme based on selection
  p <- switch(ggplot2_theme,
    "linedraw" = p + ggplot2::theme_linedraw(),
    "classic" = p + ggplot2::theme_classic(),
    "gray" = p + ggplot2::theme_gray(),
    "light" = p + ggplot2::theme_light(),
    "dark" = p + ggplot2::theme_dark(),
    "bw" = p + ggplot2::theme_bw(),
    p # default, no change
  )

  if (!gridline) { # by default it has gridlines
    p <- p +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
  }

  return(p)
}

#' Remove list in a data frame
#'
#' Remove any list elements in a data frame
#'
#' @param data_object
#'
#' @return Data frame
data_frame_with_list <- function(data_object) {
  set_lists_to_chars <- function(x) {
    if (class(x) == "list") {
      y <- paste(unlist(x[1]), sep = "", collapse = ", ")
    } else {
      y <- x
    }
    return(y)
  }
  new_frame <- data.frame(
    lapply(data_object, set_lists_to_chars),
    stringsAsFactors = F
  )
  return(new_frame)
}
