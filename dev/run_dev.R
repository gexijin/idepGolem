# Set options here
{
  options(golem.app.prod = FALSE) # TRUE = production mode, FALSE = development mode

  # Detach all loaded packages and clean your environment
  golem::detach_all_attached()

  # Document and reload your package
  golem::document_and_reload()

#  styler::style_pkg()

  # Run the application
  run_app()
}
