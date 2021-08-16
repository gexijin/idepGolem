#' fct_07_deg2.R This file holds all of the main data analysis functions
#' associated with seventh tab of the iDEP website.
#'
#'
#' @section fct_07_deg2.R functions:
#' \code{short_species_names} makes short specie names for string api
#'
#'
#' @name fct_07_deg2.R
NULL



# Homo sapies --> hsapiens
# used for string api
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param specie_mame DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
short_species_names <- function(specie_mame) {
  specie_mame <- strsplit(x = as.character(specie_mame), split = " ") # nolint
  return(tolower(paste0( # nolint
    substr(x = specie_mame[[1]][1], start = 1, stop = 1),
    specie_mame[[1]][2]
  )))
}
