#' Utilities for dynamically sizing sample-heavy plots
#'
#' These helpers centralize the logic for creating horizontally scrollable
#' plot outputs whose width adapts to the number of samples being displayed.
#' They are deliberately lightweight so they can be reused across modules.

#' Compute a pixel width for sample-based plots
#'
#' @param n_samples Number of samples to display on the x-axis.
#' @param min_width Minimum width (in pixels) to use regardless of sample count.
#' @param per_sample Additional pixels to allocate per sample.
#' @param max_width Maximum width (in pixels) to avoid runaway sizes.
#' @return Numeric pixel width respecting the provided bounds.
#' @noRd
compute_scrollable_plot_width <- function(n_samples,
                                          min_width = 700,
                                          per_sample = 10,
                                          max_width = 60000) {
  if (is.null(n_samples) || !is.finite(n_samples) || n_samples <= 0) {
    return(NULL)
  }

  estimated_width <- min_width + (as.numeric(n_samples) * per_sample)
  max(min_width, min(estimated_width, max_width))
}

#' Wrap a plotOutput with a scrollable container
#'
#' @param outputId Plot output id (already namespaced when used inside modules).
#' @param ... Additional arguments forwarded to [shiny::plotOutput()].
#' @return A `div` containing the plot output with horizontal scrolling enabled.
#' @noRd
scrollable_plot_output <- function(outputId, ...) {
  shiny::div(
    class = "idep-scrollable-plot",
    shiny::plotOutput(outputId = outputId, ...)
  )
}

#' Send a sizing instruction for a scrollable plot
#'
#' @param session Shiny session.
#' @param output_id Plot output id (without namespace inside modules).
#' @param n_samples Number of samples to inform the width calculation.
#' @param min_width,per_sample,max_width See [compute_scrollable_plot_width()].
#' @noRd
update_scrollable_plot_width <- function(session,
                                         output_id,
                                         n_samples,
                                         min_width = 700,
                                         per_sample = 14,
                                         max_width = 6000) {
  if (is.null(session) || is.null(output_id)) {
    return(invisible(NULL))
  }

  width_px <- compute_scrollable_plot_width(
    n_samples = n_samples,
    min_width = min_width,
    per_sample = per_sample,
    max_width = max_width
  )

  session$onFlushed(function() {
    session$sendCustomMessage(
      type = "idep-set-plot-width",
      message = list(
        id = session$ns(output_id),
        width = width_px
      )
    )
  }, once = TRUE)

  invisible(width_px)
}
