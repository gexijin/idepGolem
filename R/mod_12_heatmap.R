#' 12_heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_12_heatmap_ui <- function(id){
  ns <- NS(id)
  tagList(

    fluidRow(
      column(
        width = 3,
        selectInput(
          inputId = ns("heatmap_color_select"),
          label = NULL,
          choices = "green-black-red",
          width = "100%"
        ),
        plotOutput(
          outputId = ns("main_heatmap"),
          height = "450px",
          width = "100%",
          brush = ns("ht_brush")
        ),
        br(),
        h5("Selected Cell (Submap):"),
        uiOutput(
          outputId = ns("ht_click_content")
        )
      ),
      column(
        width = 9,
        plotOutput(
          outputId = ns("sub_heatmap"),
          height = "650px",
          width = "100%",
          click = ns("ht_click")
        )
      )
    )
  )
}
    
#' 12_heatmap Server Functions
#'
#' @noRd 
mod_12_heatmap_server <- function(
  id,
  data,
  bar,
  all_gene_names,
  cluster_rows = FALSE
){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    # Interactive heatmap environment
    shiny_env <- new.env()

    # Heatmap Colors ----------
    heatmap_colors <- list(
      "Green-Black-Red" = c("green", "black", "red"),
      "Red-Black-Green" = c("red", "black", "red"), 
      "Blue-White-Red" = c("blue", "white", "red"),
      "Green-Black-Magenta" = c("green", "black", "magenta"),
      "Blue-Yellow-Red" = c("blue", "yellow", "red"),
      "Blue-White-Brown" = c("blue", "white", "brown"), 
      "Orange-White-Blue" = c("orange", "white", "blue")
    )

    heatmap_choices <- c(
      "Green-Black-Red",
      "Red-Black-Green", 
      "Blue-White-Red",
      "Green-Black-Magenta",
      "Blue-Yellow-Red",
      "Blue-White-Brown", 
      "Orange-White-Blue"
    )

    observe({
      updateSelectInput(
        session = session,
        inputId = "heatmap_color_select",
        choices = heatmap_choices
      )
    })

    output$main_heatmap <- renderPlot({
      req(!is.null(data()))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Creating Heatmap",
        color = "#000000"
      )

      # Assign heatmap to be used in multiple components
      shiny_env$ht <- deg_heatmap(
        data = data(),
        bar = bar,
        heatmap_color_select = heatmap_colors[[input$heatmap_color_select]],
        cluster_rows = cluster_rows
      )

      # Use heatmap position in multiple components
      shiny_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht)

      shinybusy::remove_modal_spinner()

      return(shiny_env$ht)
    })

    output$sub_heatmap <- renderPlot({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("Select a region on the heatmap to zoom in.", 0.5, 0.5)
      } else {
        heat_return <- deg_heat_sub(
          ht_brush = input$ht_brush,
          ht = shiny_env$ht,
          ht_pos_main = shiny_env$ht_pos_main,
          heatmap_data = data(),
          bar = bar,
          all_gene_names = all_gene_names()
        )

        shiny_env$ht_select <- heat_return$ht_select
        shiny_env$submap_data <- heat_return$submap_data
        shiny_env$group_colors <- heat_return$group_colors
        shiny_env$column_groups <- heat_return$column_groups
        shiny_env$bar <- heat_return$bar

        shiny_env$ht_sub <- ComplexHeatmap::draw(
          shiny_env$ht_select,
          annotation_legend_side = "top",
          heatmap_legend_side = "top"
        )

        shiny_env$ht_pos_sub <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht_sub)

        return(shiny_env$ht_sub)
      }
    })

    # Sub Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      if (is.null(input$ht_click)) { 
        "Click for Info."
      } else {
        deg_click_info(
          click = input$ht_click,
          ht_sub = shiny_env$ht_sub,
          ht_sub_obj = shiny_env$ht_select,
          ht_pos_sub = shiny_env$ht_pos_sub,
          sub_groups = shiny_env$column_groups,
          group_colors = shiny_env$group_colors,
          bar = shiny_env$bar,
          data = shiny_env$submap_data
        )
      }
    })



  })
}
    
## To be copied in the UI
# mod_12_heatmap_ui("12_heatmap_1")
    
## To be copied in the server
# mod_12_heatmap_server("12_heatmap_1")
