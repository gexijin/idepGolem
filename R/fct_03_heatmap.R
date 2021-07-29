#' 03_heatmap 
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#'
#' @noRd

# adding sample legends to heatmap; this is for the main heatmap
# https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
add_legend <- function(...) {
  opar <- par(fig = c(0, 1, 0, 1),
      oma = c(0, 0, 0, 0),
      mar = c(0, 0, 0, 6),
      new = TRUE)
  on.exit(par(opar))
  plot(x = 0, y = 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend(...)
}