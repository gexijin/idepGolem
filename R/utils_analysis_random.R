#' utils_analysis_random.R
#'
#'
#' @section utils_analysis_random.R functions:
#'
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
    averge = hclust_average,
    ward_D = hclust_ward_D,
    ward_D2 = hclust_ward_D2,
    single = hclust_single,
    mcquitty = hclust_mcquitty,
    median = hclust_median,
    centroid = hclust_centroid
  ))
}


#' hcluster_functions
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
        FUN = foo)
      )
      names(sample_group) <- row.names(sample_info2)
      if (min(table(sample_group)) == 1) { # no replicates?
        sample_group <- sample_info2[, 1]
      }
    }
  }
  return(as.character(sample_group))
}


# Clean up gene sets. Remove spaces and other control characters from gene names  
clean_gene_set <- function (gene_set) {
  # remove duplicate; upper case; remove special characters
  gene_set <- unique(toupper(gsub("\n| ", "", gene_set)))
  # genes should have at least two characters
  gene_set <- gene_set[which(nchar(gene_set) > 1)]
  return(gene_set)
}


# Read gene sets GMT file
# This functions cleans and converts to upper case
read_gmt <- function(file_path) { # size restriction
  # Read in the first file
  gmt_data <- scan(file = file_path, what = "", sep = "\n")
  gmt_data <- gsub(pattern = " ",
    replacement = "",
    x = gmt_data) 
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

