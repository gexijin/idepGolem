#' 08_bicluster UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_08_bicluster_ui <- function(id){
  ns <- NS(id)
  tabPanel(
    "Bicluster",
    sidebarLayout(
      sidebarPanel(
        h5(
          "Biclustering can discover genes correlated on subset of samples. 
           Only useful when  sample size is large(>10). 
           Uses methods implemented in the biclust R package."
        ),
        numericInput(
          inputId = ns("n_genes"), 
          label = h5("Most variable genes to include "), 
          min = 10, 
          max = 2000, 
          value = 1000
        ),
        selectInput(
          inputId = ns("biclust_method"),
          label = "Method:", 
          choices = list( 
            "BCCC" = "biclust::BCCC()",
            "QUBIC" = "QUBIC::BCQU()",
            "runibic" = "runibic::BCUnibic()",
            "BCXmotifs" = "biclust::BCXmotifs()",
            "BCPlaid" = "biclust::BCPlaid()",
            "BCSpectral" = "biclust::BCSpectral()",
            "BCBimax" = "biclust::BCBimax()",
            "BCQuest" = "biclust::BCQuest()"
          ),
          selected = "BCCC()"
        ),
        selectInput(
          inputId = ns("heatmap_color_select"),
          label = "Select Heatmap Color: ",
          choices = "green-black-red",
          width = "100%"
        ),
        htmlOutput(outputId = ns("list_biclusters")),
        h5("Enrichment database"),
        htmlOutput(outputId = ns("select_go_selector")),
        tags$style(
          type = "text/css",
          "#bicluster-select_go{ width:100%;   margin-top:-9px}"
        ),
        br(),
        br(),
        
        # Show biclust message
        actionButton(
          inputId = ns("show_messages"),
          label = "Show Biclust Summary"
        ),

        a(
          h5(
            "Questions?",
            align = "right"
          ),
          href = "https://idepsite.wordpress.com/biclustering/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Heatmap",
            fluidRow(
              column(
                width = 3,
                plotOutput(
                  outputId = ns("biclust_main_heatmap"),
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
                  outputId = ns("biclust_sub_heatmap"),
                  height = "650px",
                  width = "100%",
                  click = ns("ht_click")
                )
              )
            )
          ),
          tabPanel(
            "Enrichment",
            fluidRow(
              column(
                width = 4,
                checkboxInput(
                  inputId = ns("filtered_background"), 
                  label = "Use filtered data as background in enrichment (slow)", 
                  value = TRUE
                )
              ),
              column(
                width = 4,
                checkboxInput(
                  inputId = ns("remove_redudant"),
                  label = "Remove Redudant Gene Sets",
                  value = FALSE
                )
              )
            ),
            h3("Enriched gene sets in selected bicluster"),
            DT::dataTableOutput(
              outputId = ns("pathway_data_biclust")
            )
          ),
          tabPanel(
            "Cluster Gene Table",
            DT::dataTableOutput(
              outputId = ns("gene_list_bicluster")
            )
          )
        )
      )
    )
  )
}
    
#' 08_bicluster Server Functions
#'
#' @noRd 
mod_08_bicluster_server <- function(id, pre_process, idep_data, tab){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

    # Interactive heatmap environment
    biclust_env <- new.env()

    # GMT choices for enrichment ----------
    output$select_go_selector <- renderUI({
	    req(!is.null(pre_process$gmt_choices()))

	    selectInput(
        inputId = ns("select_go"),
        label = "Select Geneset:",
        choices = pre_process$gmt_choices(),
        selected = "GOBP"
      )
    })

    biclustering <- reactive({
      req(!is.null(pre_process$data()))

      get_biclustering(
        data = pre_process$data(),
        n_genes = input$n_genes,
        biclust_method = input$biclust_method
      )
    })

    output$list_biclusters <- renderUI({
      req(tab() == "Bicluster")
		  req(!is.null(biclustering()))
      req(biclustering()$res@Number != 0)
		
			selectInput(
        inputId = ns("select_bicluster"), 
				label = "Select a cluster",
			  choices = 1:biclustering()$res@Number  
			)		
	  })

    # Heatmap Colors ----------
    heatmap_colors <- list(
      "Green-Black-Red" = c("green", "black", "red"),
      "Blue-White-Red" = c("blue", "white", "red"),
      "Green-Black-Magenta" = c("green", "black", "magenta"),
      "Blue-Yellow-Red" = c("blue", "yellow", "red"),
      "Blue-White-Brown" = c("blue", "white", "brown")
    )
    heatmap_choices <- c(
      "Green-Black-Red",
      "Blue-White-Red",
      "Green-Black-Magenta",
      "Blue-Yellow-Red",
      "Blue-White-Brown"
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "heatmap_color_select",
        choices = heatmap_choices
      )
    })

    biclust_data <- reactive({
      req(!is.null(biclustering()) && !is.null(input$select_bicluster))
      req(biclustering()$res@Number != 0)

      return(
        biclust::bicluster(
          biclustering()$data,
          biclustering()$res,
          as.numeric(input$select_bicluster)
        )[[1]]
      )
    })

    output$biclust_main_heatmap <- renderPlot({
      req(!is.null(biclust_data()))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Creating Heatmap",
        color = "#000000"
      )

      # Assign heatmap to be used in multiple components
      biclust_env$ht <- basic_heatmap(
        data = biclust_data(),
        heatmap_color_select = heatmap_colors[[input$heatmap_color_select]]
      )

      # Use heatmap position in multiple components
      biclust_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(biclust_env$ht)

      shinybusy::remove_modal_spinner()

      return(biclust_env$ht)
    })

    output$biclust_sub_heatmap <- renderPlot({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("No region is selected.", 0.5, 0.5)
      } else {
        biclust_heat_return <- basic_heat_sub(
          ht_brush = input$ht_brush,
          ht = biclust_env$ht,
          ht_pos_main = biclust_env$ht_pos_main,
          heatmap_data = biclust_data()
        )

        biclust_env$ht_select <- biclust_heat_return$ht_select
        biclust_env$submap_data <- biclust_heat_return$submap_data
        biclust_env$group_colors <- biclust_heat_return$group_colors
        biclust_env$column_groups <- biclust_heat_return$column_groups
        
        biclust_env$ht_sub <- ComplexHeatmap::draw(
          biclust_env$ht_select,
          annotation_legend_side = "top",
          heatmap_legend_side = "top"
        )

        biclust_env$ht_pos_sub <- InteractiveComplexHeatmap::htPositionsOnDevice(
          biclust_env$ht_sub
        )

        return(biclust_env$ht_sub)
      }
    })

    # Sub Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      if (is.null(input$ht_click)) { 
        "Click for Info."
      } else {
        heat_click_info(
          click = input$ht_click,
          ht_sub = biclust_env$ht_sub,
          ht_sub_obj = biclust_env$ht_select,
          ht_pos_sub = biclust_env$ht_pos_sub,
          sub_groups = biclust_env$column_groups,
          group_colors = biclust_env$group_colors,
          data = biclust_env$submap_data
        )
      }
    })

    # Biclustering summary message -----------
    bicluster_info <- reactive({		
      req(!is.null(biclustering()) && !is.null(input$select_bicluster))

      bicluster_summary_message(
        biclustering = biclustering(),
        select_bicluster = input$select_bicluster
      )	
	  }) 

    # Show messages when on the Bicluster tab or button is clicked
    observe({
      req(input$show_messages || tab() == "Bicluster")
      req(!is.null(bicluster_info()))

      showNotification(
        ui = bicluster_info(),
        id = "biclust_summary",
        duration = NULL,
        type = "default"
      )
    })

    pathway_table_biclust <- reactive({
      req(!is.null(biclust_data()))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Running Analysis",
        color = "#000000"
      )

      gene_names <- merge_data(
        all_gene_names = pre_process$all_gene_names(),
        data = biclust_data(),
        merge_ID = "ensembl_ID"
      )
      # Only keep the gene names and scrap the data
      gene_names_query <- dplyr::select_if(gene_names, is.character)

      req(!is.null(input$select_go))

      gene_sets <- read_pathway_sets(
        all_gene_names_query = gene_names_query,
        converted = pre_process$converted(),
        go = input$select_go,
        select_org = pre_process$select_org(),
        gmt_file = pre_process$gmt_file(),
        idep_data = idep_data,
        gene_info = pre_process$all_gene_info()
      )

      pathway_info <- find_overlap(
        pathway_table = gene_sets$pathway_table,
        query_set = gene_sets$query_set,
        total_genes = gene_sets$total_genes,
        processed_data = pre_process$data(),
        gene_info = pre_process$all_gene_info(),
        go = input$select_go,
        idep_data = idep_data,
        select_org = pre_process$select_org(),
        sub_pathway_files = gene_sets$pathway_files,
        use_filtered_background = input$filtered_background,
        reduced = input$remove_redudant
      )

      shinybusy::remove_modal_spinner()

      return(pathway_info)
    })

    # Enrichment Data Table ----------
    output$pathway_data_biclust <- DT::renderDataTable({
      req(!is.null(pathway_table_biclust()))

      if(ncol(pathway_table_biclust()) > 1) {
        pathway_table <- pathway_table_biclust()[, 1:4]
      } else {
        pathway_table <- pathway_table_biclust()
      }

      DT::datatable(
        pathway_table,
        options = list(
          pageLength = 20,
          scrollX = "400px"
        ),
        rownames = FALSE
      )
    })

    output$gene_list_bicluster <- DT::renderDataTable({
      req(!is.null(biclust_data()))

      biclust_table <- get_biclust_table_data(
        res = biclustering()$res,
        biclust_data = biclust_data(),
        select_go = input$select_go,
        select_org = pre_process$select_org(),
        all_gene_info = pre_process$all_gene_info()
      )

      DT::datatable(
        biclust_table,
        options = list(
          pageLength = 20,
          scrollX = "400px"
        ),
        rownames = FALSE
      )
    })
  })
}
    
## To be copied in the UI
# mod_08_bicluster_ui("08_bicluster_ui_1")
    
## To be copied in the server
# mod_08_bicluster_server("08_bicluster_ui_1")