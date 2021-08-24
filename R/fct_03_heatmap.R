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

#' Get color choice list for the heatmap
#'
#' get_heatmap_colors
#'
#' Calling the function will provide a list of heatmap colors
#' for the select input bar to use and gives the colors that the
#' heatmap will use for the selected color set.
#'
#' @return Returns a list containing color choices and a data frame
#' of hex color values.
get_heatmap_colors <- function() {
  hmcols <- colorRampPalette(rev(c(
    "#D73027", "#FC8D59", "#FEE090", "#FFFFBF",
    "#E0F3F8", "#91BFDB", "#4575B4"
  )))(75)
  heat_colors <- rbind(
    gplots::greenred(75), gplots::bluered(75),
    gplots::colorpanel(75, "green", "black", "magenta"),
    gplots::colorpanel(75, "blue", "yellow", "red"), hmcols
  )
  rownames(heat_colors) <- c(
    "Green-Black-Red", "Blue-White-Red", "Green-Black-Magenta",
    "Blue-Yellow-Red", "Blue-white-brown"
  )
  color_choices <- setNames(1:dim(heat_colors)[1], rownames(heat_colors))

  return(
    list(
      color_choices = color_choices,
      heat_colors = heat_colors
    )
  )
}
