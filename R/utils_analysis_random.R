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
clean_gene_set <- function (gene_set){
  # remove duplicate; upper case; remove special characters
  gene_set <- unique(toupper(gsub("\n| ", "", gene_set)))
  # genes should have at least two characters
  gene_set <- gene_set[which(nchar(gene_set) > 1)]
  return(gene_set)
}


#HOW YOU WORK?????
# read gene set files in the GMT format, does NO cleaning. Assumes the GMT files are created with cleanGeneSet()
# See http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Gene_Set_Database_Formats
read_gmt <- function(file_name) {
  gmt_data <- scan(file = file_name, what = "", sep = "\n")
  gmt_data <- strsplit(x = gmt_data, split = "\t")
  # Extract the first vector element and set it as the list element name
  names(gmt_data) <- sapply(X = gmt_data, `[[`, 1)
  # 2nd element is comment, ignored
  gmt_data <- lapply(X = gmt_data, `[`, -c(1, 2)) 
  # gene sets smaller than 1 is ignored!!!
  gmt_data <- gmt_data[which(sapply(X = gmt_data, FUN = length) > 1)]
  return(gmt_data)
}