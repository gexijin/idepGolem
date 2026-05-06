#' Version check utilities
#'
#' Compare the locally installed idepGolem version against the Version field
#' in the DESCRIPTION file on the master branch of gexijin/idepGolem and
#' surface a notification when a newer release is available.
#'
#' @noRd

# Package-internal cache so the network call only happens once per R session
.idep_version_cache <- new.env(parent = emptyenv())

#' Fetch and compare the remote DESCRIPTION version
#'
#' @param url Raw URL to the DESCRIPTION file on the default branch.
#' @param timeout_seconds Hard cap on the HTTP read; failures degrade silently.
#' @param force If TRUE, bypass the in-memory cache.
#' @return list(status, local, remote) where status is one of
#'   "current", "outdated", "ahead", "error".
#' @noRd
check_idep_version_remote <- function(
    url = "https://raw.githubusercontent.com/gexijin/idepGolem/master/DESCRIPTION",
    timeout_seconds = 5,
    force = FALSE) {

  if (!force && !is.null(.idep_version_cache$result)) {
    return(.idep_version_cache$result)
  }

  local <- tryCatch(
    as.character(utils::packageVersion("idepGolem")),
    error = function(e) NA_character_
  )

  remote <- tryCatch({
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = timeout_seconds)
    con <- url(url, "r")
    on.exit(close(con), add = TRUE)
    lines <- readLines(con, n = 50, warn = FALSE)
    version_line <- grep("^Version:\\s*", lines, value = TRUE)
    if (length(version_line) == 0) {
      stop("Version field not found in remote DESCRIPTION")
    }
    sub("^Version:\\s*", "", version_line[1])
  }, error = function(e) NA_character_)

  status <- if (is.na(local) || is.na(remote)) {
    "error"
  } else {
    cmp <- utils::compareVersion(local, remote)
    if (cmp < 0) "outdated" else if (cmp > 0) "ahead" else "current"
  }

  result <- list(status = status, local = local, remote = remote)
  .idep_version_cache$result <- result
  result
}

#' Show a Shiny notification when a newer iDEP version is available
#'
#' Reads the cached result from check_idep_version_remote() and emits a
#' dismissible warning notification. No-op unless status is "outdated".
#' @noRd
show_outdated_version_notification <- function() {
  result <- check_idep_version_remote()
  if (!identical(result$status, "outdated")) {
    return(invisible(NULL))
  }

  shiny::showNotification(
    shiny::tags$div(
      "A newer version of iDEP (",
      shiny::tags$strong(result$remote),
      ") is available. You are running ",
      shiny::tags$strong(result$local), ".",
      shiny::tags$br(),
      "Update from R with ",
      shiny::tags$code(
        'devtools::install_github("gexijin/idepGolem", upgrade = "never")'
      ),
      " or see the ",
      shiny::tags$a(
        href = "https://github.com/gexijin/idepGolem",
        target = "_blank",
        "GitHub README"
      ),
      "."
    ),
    type = "warning",
    duration = 20
  )
}
