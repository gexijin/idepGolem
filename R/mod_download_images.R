#' download_images UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_download_images_ui <- function(id){
  ns <- NS(id)
  tagList(
    actionButton(
      inputId = ns("download_popup"), 
      label = "Download Plot"
    )
  )
}


    
#' download_images Server Functions
#'
#' @noRd 
mod_download_images_server <- function(id, filename, figure){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    observeEvent(
      
      
      input$download_popup, 
      {
        showModal(modalDialog(
          numericInput(
            inputId = ns("width"), 
            label = "Width (in)", 
            value = 5, 
            min = 1, 
            max = 100
          ), 
          numericInput(
            inputId = ns("height"), 
            label = "Height (in)", 
            value = 4, 
            min = 1, 
            max = 100
          ),
          downloadButton(
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
          )
        ))
      },
      
    )
    
   
    output$dl_pdf <-downloadHandler(
      filename = paste0(filename, ".pdf"), 
      content = function(file){
        pdf(
          file,
          width = min(100, input$width, na.rm = TRUE), 
          height = min(100, input$height, na.rm = TRUE)
        )
        print(
          figure
        )
        dev.off()
      }
    )

    output$dl_png <- downloadHandler(
      filename = paste0(filename, ".png"), 
      content = function(file){
        png(
          file, 
          res = 360, 
          width = min(100, input$width, na.rm = TRUE), 
          height = min(100, input$height, na.rm = TRUE),
          units = "in"
        )
        print(
          figure
        )
        dev.off()
      }
    )
    
    output$dl_svg <- downloadHandler(
      filename = paste0(filename, ".svg"), 
      content = function(file){
        svg(
          file, 
          width = min(100, input$width, na.rm = TRUE), 
          height = min(100, input$height, na.rm = TRUE)
        )
        print(
          figure
        )
        dev.off()
      }
    )
 
  })
}
    
## To be copied in the UI
# mod_download_images_ui("download_images_ui_1")
    
## To be copied in the server
# mod_download_images_server("download_images_ui_1")
