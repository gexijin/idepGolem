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
         plotOutput(
          outputId = ns("heatmap_main"),
          height = "500px",
          width = "100%",
          brush = ns("ht_brush")
        ),
        tableOutput(ns("heat_data"))
 
  )
}
    
#' 12_heatmap Server Functions
#'
#' @noRd 
mod_12_heatmap_server <- function(
  id,
  heatmap_data,
  cluster_meth,
  heatmap_cutoff,
  sample_info,
  select_factors_heatmap,
  dist_funs,
  dist_function,
  hclust_function,
  sample_clustering,
  heatmap_color_select,
  row_dend,
  k_clusters,
  re_run
){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
    # Interactive heatmap environment
    shiny_env <- new.env()

    output$heat_data <- renderTable({
      req(!is.null(heatmap_data()))
      head(heatmap_data())
    })

    # HEATMAP -----------
    # Information on interactivity
    # https://jokergoo.github.io/2020/05/15/interactive-complexheatmap/
    output$heatmap_main <- renderPlot({

      req(!is.null(heatmap_data()))
#      req(!is.null(select_factors_heatmap))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Creating Heatmap",
        color = "#000000"
      )

      # Assign heatmap to be used in multiple components
      shiny_env$ht <- heatmap_main(
        data = heatmap_data(),
        cluster_meth = cluster_meth,
        heatmap_cutoff = heatmap_cutoff,
        sample_info = sample_info(),
        select_factors_heatmap = select_factors_heatmap,
        dist_funs = dist_funs,
        dist_function = dist_function,
        hclust_function = hclust_function,
        sample_clustering = sample_clustering,
        heatmap_color_select = heatmap_color_select,
        row_dend = row_dend,
        k_clusters = k_clusters,
        re_run = re_run
      )

      # Use heatmap position in multiple components
      shiny_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht)

      shinybusy::remove_modal_spinner()

      return(shiny_env$ht)
    })

    # Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      
      if (is.null(input$ht_click) | is.null(shiny_env$ht_sub)) { 
        "Click on zoomed heatmap"
      } else {
        cluster_heat_click_info(
          click = input$ht_click,
          ht_sub = shiny_env$ht_sub,
          ht_sub_obj = shiny_env$ht_sub_obj,
          ht_pos_sub = shiny_env$ht_pos_sub,
          sub_groups = shiny_env$sub_groups,
          group_colors = shiny_env$group_colors,
          cluster_meth = input$cluster_meth,
          click_data = shiny_env$click_data
        )
      }
    })

    # Subheatmap creation ---------
    output$sub_heatmap <- renderPlot({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("Select a region on the heatmap to zoom in. 
        Gene IDs shows up when less than 60 genes are selected.", 0.5, 0.5)
      } else {
        submap_return <- heat_sub(
          ht_brush = input$ht_brush,
          ht = shiny_env$ht,
          ht_pos_main = shiny_env$ht_pos_main,
          heatmap_data = heatmap_data(),
          sample_info = pre_process$sample_info(),
          select_factors_heatmap = input$select_factors_heatmap,
          cluster_meth = input$cluster_meth
        )

        # Objects used in other components ----------
        shiny_env$ht_sub_obj <- submap_return$ht_select
        shiny_env$submap_data <- submap_return$submap_data
        shiny_env$sub_groups <- submap_return$sub_groups
        shiny_env$group_colors <- submap_return$group_colors
        shiny_env$click_data <- submap_return$click_data
        
        shiny_env$ht_sub <- ComplexHeatmap::draw(
          shiny_env$ht_sub_obj,
          annotation_legend_list = submap_return$lgd,
          annotation_legend_side = "top"
        )

        shiny_env$ht_pos_sub <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht_sub)

        return(shiny_env$ht_sub)
      }
    })


  })
}
    
## To be copied in the UI
# mod_12_heatmap_ui("12_heatmap_1")
    
## To be copied in the server
# mod_12_heatmap_server("12_heatmap_1")
