#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(onStart = NULL,
                    options = list(),
                    enableBookmarking = NULL,
                    uiPattern = "/",
                    ...) {
  # === TIMING: App startup ===
  .app_start_time <- Sys.time()
  cat(sprintf("\n"))
  cat(sprintf("========================================\n"))
  cat(sprintf("iDEP APP STARTUP TIMING\n"))
  cat(sprintf("========================================\n"))
  cat(sprintf("[%s] APP: run_app() called\n", format(.app_start_time, "%H:%M:%OS3")))

  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}
