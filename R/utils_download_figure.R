#' Download Figure UI Function
#' Adapted from ottoPlots by Emma Spors.
#'
#' @description A shiny module that provides a "Download" button for plots.
#' Works with \code{\link{mod_download_figure_server}()}. On click, opens a
#' pop-up where users can set width/height and download as PDF, PNG, or SVG.
#' Supports both {ggplot2} and base R plots (base R plots must be captured
#' with \code{\link{recordPlot}()}).
#'
#' @param id String. Namespace id, internal parameter for {shiny}.
#'
#' @return UI elements for the download button.
#'
#' @noRd
#'
#' @importFrom shiny NS uiOutput
mod_download_figure_ui <- function(id) {
  ns <- NS(id)

  uiOutput(ns("download"))
}

#' Download Figure Server Function
#'
#' @description Server component of the download figure module. Works with
#' \code{\link{mod_download_figure_ui}()} to handle the download pop-up modal
#' and write the plot to PDF, PNG, or SVG. Supports both {ggplot2} and base R
#' plots (base R plots must be captured with \code{\link{recordPlot}()}).
#'
#' @param id String. Namespace id, internal parameter for {shiny}.
#' @param filename Character string for the output filename (without extension).
#' @param figure Reactive expression returning the plot object to download.
#' @param width Default plot width in inches. Must be between 2 and 30.
#'   Default is 8.
#' @param height Default plot height in inches. Must be between 2 and 30.
#'   Default is 6.
#' @param label Character string for the download button label.
#'   Default is "Download Plot".
#'
#' @return {shiny} server logic for the download pop-up modal.
#'
#' @noRd
#'
#' @importFrom grDevices dev.off pdf png svg
#'
mod_download_figure_server <- function(id,
                                       filename,
                                       figure,
                                       width = 8,
                                       height = 6,
                                       label = "Download Plot") {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$download <- renderUI({
      req(figure())
      tagList(
        actionButton(
          inputId = ns("download_popup"),
          icon = icon("download"),
          label = label
        ),
        tippy::tippy_this(
          ns("download_popup"),
          "Click to download plot in preferred format and size."
        )
      )
    })

    min_size <- 2
    max_size <- 30

    figure_width <- reactive({
      if (is.numeric(input$width)) {
        return(max(min_size, min(max_size, input$width, na.rm = TRUE)))
      } else {
        return(width)
      }
    })

    figure_height <- reactive({
      if (is.numeric(input$width)) {
        return(max(min_size, min(max_size, input$height, na.rm = TRUE)))
      } else {
        return(height)
      }
    })

    observeEvent(input$download_popup, {
      showModal(
        modalDialog(
          numericInput(
            inputId = ns("width"),
            label = "Width (in)",
            value = width,
            min = min_size,
            max = max_size
          ),
          numericInput(
            inputId = ns("height"),
            label = "Height (in)",
            value = height,
            min = min_size,
            max = max_size
          ),
          h5(
            "The plot will be rendered differently depending on size.",
            "When the dimensions are too small, an error or blank plot",
            "will be generated."
          ),
          downloadButton(outputId = ns("dl_pdf"), label = "PDF"),
          downloadButton(outputId = ns("dl_png"), label = "PNG"),
          downloadButton(outputId = ns("dl_svg"), label = "SVG"),
          size = "s"
        )
      )
    })

    output$dl_pdf <- downloadHandler(
      filename = paste0(filename, ".pdf"),
      content = function(file) {
        on.exit(removeModal())
        pdf(file, width = figure_width(), height = figure_height())
        print(figure())
        dev.off()
      }
    )

    output$dl_png <- downloadHandler(
      filename = paste0(filename, ".png"),
      content = function(file) {
        on.exit(removeModal())
        png(file, res = 360, width = figure_width(), height = figure_height(), units = "in")
        print(figure())
        dev.off()
      }
    )

    output$dl_svg <- downloadHandler(
      filename = paste0(filename, ".svg"),
      content = function(file) {
        on.exit(removeModal())
        svg(file, width = figure_width(), height = figure_height())
        print(figure())
        dev.off()
      }
    )
  })
}
