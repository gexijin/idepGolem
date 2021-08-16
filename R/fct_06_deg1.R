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
