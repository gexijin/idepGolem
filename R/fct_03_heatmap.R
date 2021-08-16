#' fct_03_heatmap.R This file holds all of the main data analysis functions
#' associated with third tab of the iDEP website.
#'
#'
#' @section fct_03_heatmap.R functions:
#' \code{add_legend}
#'
#'
#' @name fct_03_heatmap.R
NULL


# adding sample legends to heatmap; this is for the main heatmap
# https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param ... DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
add_legend <- function(...) {
  opar <- par(
    fig = c(0, 1, 0, 1),
    oma = c(0, 0, 0, 0),
    mar = c(0, 0, 0, 6),
    new = TRUE
  )
  on.exit(par(opar))
  plot(x = 0, y = 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend(...)
}
