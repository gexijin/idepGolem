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

  # Initialize database configuration - shared across all user sessions
  # Using onStart callback ensures this runs once when app starts, not per user
  # See: https://shiny.posit.co/r/articles/improve/scoping/

  init_app <- function() {
    # Global database configuration
    db_ver <<- "data113"
    db_url <<- "http://bioinformatics.sdstate.edu/data/"

    # if environmental variable is not set, use relative path
    DATAPATH <<- Sys.getenv("IDEP_DATABASE")[1]
    # if not defined in the environment, use two levels above
    if (nchar(DATAPATH) == 0) {
      DATAPATH <<- paste0("../../data/")
    }
    # Add version
    DATAPATH <<- paste0(DATAPATH, "/", db_ver, "/")
    org_info_file <<- paste0(DATAPATH, "demo/orgInfo.db")
    if (!file.exists(org_info_file)) {
      DATAPATH <<- paste0("./", db_ver, "/")
      org_info_file <<- paste0(DATAPATH, "demo/orgInfo.db")
    }

    # Load static data files (species list, demo files, etc.) once at startup
    # This is shared across all user sessions for better performance
    # See https://github.com/ThinkR-open/golem/issues/6

    idep_data <<- get_idep_data()

    # Call user's onStart if provided
    if (!is.null(onStart)) {
      if (is.function(onStart)) {
        onStart()
      }
    }
  }

  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = init_app,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}
