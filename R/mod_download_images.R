#' download_images UI Function
#'
#' @description A shiny Module to enable downloading figures. Main idea is to 
#' create the plot in a reactive function, which returns a plot object. 
#' Call this function to print to render for shiny app, 
#' or export as PDF, png, or SVG files.
#' This module can be reused throughout a Shiny app.
#' By Emma Spors, Ben Derenge 
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param label is the label for the download button
#' @noRd
#'
#' @importFrom shiny NS tagList
#' 
mod_download_images_ui <- function(id, label = "Download Plot") {
  ns <- NS(id)
  tagList(
    actionButton(
      inputId = ns("download_popup"),
      label = label
    )
  )
}


#' download_images Server Functions
#'
#' @noRd
#' @param id is module id
#' @param filename is a name of the output file without extensions
#' @param figure is a graphics object. Note that ggplot2 objects can 
#' be directly used, while base R graphics, we need to use 
#' the savePlot() function. 
#' 
#' @param width specifies a default width in inches
#' @param height specifies default height in inches
#' 
mod_download_images_server <- function(id, filename, figure, width = 8, height = 6) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    min_size <- 2 # min for width or height
    max_size <- 30 # max for width or height
    #Check width
    figure_width <- reactive({
      if (is.numeric(input$width)) {
        # cutoff the input value to a range as shiny does not enforce
        return(max(min_size, min(max_size, input$width, na.rm = TRUE)))
      } else {
        # use default value when entered text
        return(width)
      }
    })

    #Check height
    figure_height <- reactive({
      if (is.numeric(input$width)) {
        # cutoff the input value to a range as shiny does not enforce
        return(max(min_size, min(max_size, input$height, na.rm = TRUE)))
      } else {
        # use default value when entered text
        return(height)
      }
    })

    # Generate a popup UI when clicked. You can do this in server function?
    observeEvent(
      input$download_popup,
      {
        showModal(
          modalDialog(
            numericInput(               #Figure width
              inputId = ns("width"),
              label = "Width (in)",
              value = width,
              min = min_size,
              max = max_size
            ),
            numericInput(               #Figure Height
              inputId = ns("height"),
              label = "Height (in)",
              value = height,
              min = min_size,
              max = max_size
            ),
            h5("The plot will be rendered differently depending on size. When the dimensions are too small, 
            error or blank plot will be generated."),
            downloadButton(             #buttons
              outputId = ns("dl_pdf"),
              label = "PDF"
            ),
            downloadButton(
              outputId = ns("dl_png"),
              label = "PNG"
            ),
            downloadButton(
              outputId = ns("dl_svg"),
              label = "SVG"
            ),
            size = "s" # small dialog modal
          )
        )
      }
    )
    
    # Download PDF
    output$dl_pdf <- downloadHandler(
      filename = paste0(filename, ".pdf"),
      content = function(file) {
        # remove popup when downloaded
        #on.exit(removeModal())
        pdf(
          file,
          width = figure_width(),
          height = figure_height()
        )
        print(figure())
        dev.off()
      }
    )
    
    # Download PNG
    output$dl_png <- downloadHandler(
      filename = paste0(filename, ".png"),
      content = function(file) {
        # remove popup when downloaded
        #on.exit(removeModal())
        png(
          file,
          res = 360,
          width = figure_width(),
          height = figure_height(),
          units = "in"
        )
        print(figure())
        dev.off()
      }
    )

    #download SVG
    output$dl_svg <- downloadHandler(
      filename = paste0(filename, ".svg"),
      content = function(file) {
        # remove popup when downloaded
        #on.exit(removeModal())
        svg(
          file,
          width = figure_width(),
          height = figure_height()
        )
        print(figure())
        dev.off()
      }
    )
  })
}

## To be copied in the UI
# mod_download_images_ui("download_images_ui_1")

## To be copied in the server
# mod_download_images_server("download_images_ui_1")