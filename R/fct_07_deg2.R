#' fct_07_deg2.R
#'
#'
#' @section fct_07_deg2.R functions:
#' short_species_names makes short specie names for string api
#'
#'
#' @name fct_07_deg2.R
NULL



# Homo sapies --> hsapiens
# used for string api
short_species_names <- function(specie_mame) {
  specie_mame <- strsplit(x = as.character(specie_mame), split = " ") # nolint
  return(tolower(paste0( # nolint
    substr(x = specie_mame[[1]][1], start = 1, stop = 1),
    specie_mame[[1]][2]
  )))
}
