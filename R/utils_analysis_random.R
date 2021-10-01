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

#' hcluster_functions
#'
#'
#' @description
#'
#'
#' @return
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


#' dist_functions
#'
#'
#' @description
#'
#'
#' @return
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
    euclidean = dist,
    pearson_correlation = dist_pcc,
    absolute_pcc = dist_abs_pcc
  ))
}


#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param num_set DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
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


#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param sample_names DESCRIPTION.
#' @param sample_info DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
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


# Clean up gene sets. Remove spaces and other control characters from gene names
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param gene_set DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
clean_gene_set <- function(gene_set) {
  # remove duplicate; upper case; remove special characters
  gene_set <- unique(toupper(gsub("\n| ", "", gene_set)))
  # genes should have at least two characters
  gene_set <- gene_set[which(nchar(gene_set) > 1)]
  return(gene_set)
}


#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param query_input DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
clean_query <- function(query_input) {
  return(clean_gene_set(unlist(strsplit(
    x = toupper(query_input),
    split = "\t| |\n|\\,"
  ))))
}


# Read gene sets GMT file
# This functions cleans and converts to upper case
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param file_path DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
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

### edit later
# This function convert gene set names
# x="GOBP_mmu_mgi_GO:0000183_chromatin_silencing_at_rDNA"
# chromatin silencing at rDNA
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param x DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
proper <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))


#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param word_list DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
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


#' check_object_state an Utility function to simplify checking object states.
#'
#'
#' This function is a booling utility function,
#' which helps evaluates the state of objects,
#' and helps with sending messages from server to UI logic of shiny app.
#'
#'
#' @param check_exp An expression that should be evaluated
#'
#' @param true_message Message displayed on console and
#'  returned if \code{check_exp} is true
#'
#' @param false_message Optional message returned if \code{check_exp} is false,
#'  default to blank string
#'
#'
#' @return A list is returned, with the following elements:
#'  \code{bool} is either true or false depending on \code{check_exp} evaluation
#'
#'  \code{content} either message and depended on \code{check_exp} evaluation
#' @examples
#' check <- check_object_state(
#'   check_exp = (is.null(NULL)),
#'   true_message = as.data.frame("This is NULL")
#' )
#' # will eval to true and display message to console
#' if (check$bool) {
#'   return(check)
#' }
#'
#' check <- check_object_state(
#'   check_exp = (length(c(1, 2)) == 0),
#'   true_message = "this has 0 elements",
#'   false_message = "this doesn't have 0 elements"
#' )
#' # This will not return check and check$content will be false_message
#' if (check$bool) {
#'   return(check)
#' }
check_object_state <- function(
  check_exp,
  true_message,
  false_message = ""
) {
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

#' ggplot colors function
#'
#' This function will return the colors that
#' the ggplot2 package uses for plots.
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
#' @param select_gene_id Desired ID type for rownames
#'   (User_ID, ensembl_ID, symbol)
#'
#' @return Data matrix with changed rownames
rowname_id_swap <- function(
  data_matrix,
  all_gene_names,
  select_gene_id
) {
  if (select_gene_id == "User_ID" && ncol(all_gene_names) == 1) {
    return(data_matrix)
  } else if (select_gene_id == "User_ID") {
    new_data <- merge(
      all_gene_names,
      data_matrix,
      by.x = "ensembl_ID",
      by.y = "row.names",
      all.y = T
    )
    rownames(new_data) <- new_data$User_ID
    nums <- unlist(lapply(new_data, is.numeric))
    new_data <- new_data[, nums]
    new_data <- as.matrix(new_data)

    new_data <- new_data[order(-apply(
      new_data[, 1:dim(new_data)[2]],
      1,
      sd
    )), ]

    return(new_data)
  } else if (select_gene_id == "ensembl_ID") {
    return(data_matrix)
  } else if (select_gene_id == "symbol") {
    new_data <- merge(
      all_gene_names,
      data_matrix,
      by.x = "ensembl_ID",
      by.y = "row.names",
      all.y = T
    )
    rownames(new_data) <- new_data$symbol
    nums <- unlist(lapply(new_data, is.numeric))
    new_data <- new_data[, nums]
    new_data <- as.matrix(new_data)

    new_data <- new_data[order(-apply(
      new_data[, 1:dim(new_data)[2]],
      1,
      sd
    )), ]

    return(new_data)
  }
}

#' Merge data from gene info and processed data
#'
#' This function takes in the gene info data and merges it
#' with the data that has gone through the processing fcn.
#' The returned data contains the gene names as well as the
#' ensembl id in the first two columns.
#'
#' @param all_gene_names All matched gene names from idep data
#' @param data Data matrix with rownames to merge with gene names
#'
#' @return Inputted data with all gene name information.
merge_data <- function(
  all_gene_names,
  data,
  merge_ID
) {
  isolate({
    if (dim(all_gene_names)[2] == 1) {
      new_data <- round(data, 2)
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
        round(data, 2),
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
        round(data, 2),
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
