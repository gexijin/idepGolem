# Set options here
{
  # === TIMING: Development startup ===
  .dev_start_time <- Sys.time()
  cat(sprintf("\n"))
  cat(sprintf("========================================\n"))
  cat(sprintf("DEVELOPMENT MODE STARTUP\n"))
  cat(sprintf("========================================\n"))
  cat(sprintf("[%s] DEV: Starting development mode\n", format(.dev_start_time, "%H:%M:%OS3")))

  options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode

  # Detach all loaded packages and clean your environment
  cat(sprintf("[%s] DEV: Detaching packages...\n", format(Sys.time(), "%H:%M:%OS3")))
  golem::detach_all_attached()

  .detach_time <- Sys.time()
  cat(sprintf("[%s] DEV: Packages detached (%.3fs)\n",
              format(.detach_time, "%H:%M:%OS3"),
              as.numeric(difftime(.detach_time, .dev_start_time, units = "secs"))))

  # Document and reload your package
  cat(sprintf("[%s] DEV: Documenting and reloading package...\n", format(Sys.time(), "%H:%M:%OS3")))
  golem::document_and_reload()

  .reload_time <- Sys.time()
  cat(sprintf("[%s] DEV: Package reloaded (%.3fs)\n",
              format(.reload_time, "%H:%M:%OS3"),
              as.numeric(difftime(.reload_time, .detach_time, units = "secs"))))

#  styler::style_pkg()

  # Run the application
  run_app()
  .reload_time2 <- Sys.time()
  cat(sprintf("[%s] DEV: App reloaded (%.3fs)\n",
              format(.reload_time2, "%H:%M:%OS3"),
              as.numeric(difftime(.reload_time, .reload_time, units = "secs"))))
}
