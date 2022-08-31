#' Shows results from enrichment analysis of one or more lists of genes
#'
#' @description The input is a list of genes
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_11_enrichment_ui <- function(id){
  ns <- NS(id)
  library(shinyBS)
  tagList(

    fluidRow(
      column(
        width = 4,
        htmlOutput(outputId = ns("select_go_selector"))
      ),
      column(
        width = 4,
        actionButton(ns("customize_button"), "Customize")
      ),
      column(
        width = 4,
        shinyBS::bsModal(
          id = ns("modalExample"),
          title = "Options for enrichment analysis",
          trigger = ns("customize_button"),
          size = "small",
          fluidRow(
            column(
              width = 6,
              selectInput(
                inputId = ns("sort_by"),
                label = NULL,
                choices = list(
                  "Sort by FDR" = "FDR",
                  "Sort by fold enriched" = "Fold"
                ),
                selected = "FDR"
              )
            ),
            column(
              width = 6,
              downloadButton(
                outputId = ns("download_enrichment"),
                label = "Enrichment"
              )
            )
          ),
          fluidRow(
            column(
              width = 6,
              checkboxInput(
                inputId = ns("filtered_background"),
                label = "Use filtered genes as background.",
                value = TRUE
              )
            ),
            column(
              width = 6,
              checkboxInput(
                inputId = ns("remove_redudant"),
                label = "Remove Redudant Gene Sets",
                value = FALSE
              )
            )
          )
        )
      )
    ),
    tabsetPanel(
      tabPanel(
        title = "Results",
        tableOutput(ns("show_enrichment"))
      ),
      tabPanel(
        title = "Tree",
        htmlOutput(outputId = ns("select_cluster_tree")), 
        plotOutput(ns("enrichment_tree"))
      ),
      tabPanel(
        title = "Network",
        br(),
        fluidRow(
          column(
            width = 3,
            htmlOutput(outputId = ns("select_cluster_network"))
          ),
          column(
            width = 1,
            h5("Cutoff:"),
            align = "right"
          ),
          column(
            width = 2,
            numericInput(
              inputId = ns("edge_cutoff_deg"),
              label = NULL,
              value = 0.30,
              min = 0,
              max = 1,
              step = .1
            ),
            align = "left"
          ),
          column(
            width = 2,
            actionButton(
              inputId = ns("layout_vis_deg"),
              label = "Change layout"
            )
          ),
          column(
            width = 2,
            checkboxInput(
              inputId = ns("wrap_text_network_deg"),
              label = "Wrap text",
              value = TRUE
            )
          )
        ),
        visNetwork::visNetworkOutput(
          outputId = ns("vis_network_deg"),
          height = "800px",
          width = "100%"
        ),
        h6(
          "Two pathways (nodes) are connected if they share 30% (default, adjustable) or more genes.
          Green and red represents down- and up-regulated pathways. You can move the nodes by 
          dragging them, zoom in and out by scrolling, and shift the entire network by click on an 
          empty point and drag. Darker nodes are more significantly enriched gene sets. Bigger nodes
          represent larger gene sets. Thicker edges represent more overlapped genes."
        )
      )
    )

  )
}
    
#' 11_enrichment Server Functions
#' @param id module Id
#' @param results a list containing results from pathway analysis.
#'  Each element is a data frame. But the last column is a list,
#' hosting genes.
#' @noRd 
mod_11_enrichment_server <- function(
  id,
  results, # pathway table with fdr, fold, etc
  gmt_choices, # list of pathway categories "GOBP"
  gene_lists, # list of genes, each element is a list
  processed_data,
  gene_info,
  idep_data,
  select_org,
  converted,
  gmt_file
  ){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

        # GMT choices for enrichment ----------
    output$select_go_selector <- renderUI({

	    req(!is.null(gmt_choices()))

      selected = gmt_choices()[1] # default, overwrite by below
      if("GOBP" %in% gmt_choices()) {
        selected <- "GOBP"
      }
      if("KEGG" %in% gmt_choices()) {
        selected <- "KEGG"
      }
	    selectInput(
        inputId = ns("select_go"),
        label = NULL,
        choices = gmt_choices(),
        selected = selected
      )
    })


    output$select_cluster_network <- renderUI({
	    req(!is.null(enrichment_dataframe()))
      choices <- sort(unique(enrichment_dataframe()$group))
      selected <- choices[1]
      if(length(choices) > 1) {
        choices <- c("All Groups", choices)
        # if two groups, defaults to both
        if(length(choices) == 3) {
          selected <- choices[1]
        }
      }
	    selectInput(
        inputId = ns("select_cluster_network"),
        label = NULL,
        choices = choices,
        selected = selected
      )
    })

    output$select_cluster_tree <- renderUI({
	    req(!is.null(enrichment_dataframe()))
      choices <- sort(unique(enrichment_dataframe()$group))
      selected <- choices[1]
      if(length(choices) > 1) {
        choices <- c("All Groups", choices)
        # if two groups, defaults to both
        if(length(choices) == 3) {
          selected <- choices[1]
        }
      }
	    selectInput(
        inputId = ns("select_cluster_tree"),
        label = NULL,
        choices = choices,
        selected = selected
      )
    })

    output$download_enrichment <- downloadHandler(
      filename = function() {
        "enrichment.csv"
      },
      content = function(file) {
        write.csv(enrichment_dataframe(), file)
      }
    )


    # Conduct Enrichment Analysis ----------
    # returns a list object
    pathway_table <- reactive({
      req(!is.null(gene_lists()))
      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Running Analysis",
        color = "#000000"
      )

      pathway_info <- list()

        # disregard user selection use clusters for enrichment
      for (i in 1:length(gene_lists())) {
        gene_names_query <- gene_lists()[[i]]
        req(!is.null(input$select_go))
        gene_sets <- read_pathway_sets(
          all_gene_names_query = gene_names_query,
          converted = converted(), #n
          go = input$select_go,
          select_org = select_org(),
          gmt_file = gmt_file(), #n
          idep_data = idep_data,
          gene_info = gene_info()
        )

        pathway_info[[names(gene_lists())[i]]] <- find_overlap(
          pathway_table = gene_sets$pathway_table,
          query_set = gene_sets$query_set,
          total_genes = gene_sets$total_genes,
          processed_data = processed_data(),
          gene_info = gene_info(),
          go = input$select_go,
          idep_data = idep_data,
          select_org = select_org(),
          sub_pathway_files = gene_sets$pathway_files,
          use_filtered_background = input$filtered_background,
          reduced = input$remove_redudant
        )

      }
      

      shinybusy::remove_modal_spinner()

      return(pathway_info)
    })

    # returns a data frame
    enrichment_dataframe <- reactive({
      req(!is.null(pathway_table()))

      results_all <- do.call(rbind,
        #combine multiple data frames that are elements of a list
        lapply(
          names(pathway_table()),
          function(x) {
            if(ncol(pathway_table()[[x]]) == 1) {
              return(NULL)
            }
            df1 <- data_frame_with_list(pathway_table()[[x]])
            df1$group <- x
            return(df1)
          }
        )
      )

      if(!is.null(results_all)){
        if(ncol(results_all) > 1) {
          results_all <- results_all[,
            c("group",
              colnames(results_all)[1:(ncol(results_all) - 1)]
            )
          ]
        }
      }
      return(results_all)
    })

    # returns a data frame, but last column stores genes as lists
    # this is for tree and network plots
    enrichment_dataframe_for_tree <- reactive({
      req(!is.null(pathway_table()))

      results_all <- do.call(rbind,
        #combine multiple data frames that are elements of a list
        lapply(
          names(pathway_table()),
          function(x) {
            if(ncol(pathway_table()[[x]]) == 1) {
              return(NULL)
            }
            df1 <- pathway_table()[[x]]
            df1$group <- x
            return(df1)
          }
        )
      )

      if(!is.null(results_all)){
        if(ncol(results_all) > 1) {
          results_all <- results_all[,
            c("group",
              colnames(results_all)[1:(ncol(results_all) - 1)]
            )
          ]
        }
      }

      results_all <- subset(
        results_all,
        select = c(group, FDR, nGenes, Pathway, Genes)
      )
      results_all$FDR <- as.numeric(results_all$FDR)
      colnames(results_all) <- c(
        "Direction", "adj_p_val", "Pathway.size", "Pathways",  "Genes"
      )

      return(results_all)
    })

    # Enrichment Tree -----------
    output$enrichment_tree <- renderPlot({
      req(!is.null(enrichment_dataframe_for_tree()))
       req(!is.null(input$select_cluster_tree))
      enrichment_tree_plot(
        go_table = enrichment_dataframe_for_tree(),
        group = input$select_cluster_tree,
        right_margin = 45
      )
    })

    # Define a Network
    network_data_deg <- reactive({
      req(!is.null(enrichment_dataframe_for_tree()))
      req(!is.null(input$select_cluster_network))

      network_data(
        network = enrichment_dataframe_for_tree(),
        up_down_reg_deg = input$select_cluster_network,
        wrap_text_network_deg = input$wrap_text_network_deg,
        layout_vis_deg = input$layout_vis_deg,
        edge_cutoff_deg = input$edge_cutoff_deg
      )
    })

    # Interactive vis network plot
    output$vis_network_deg <- visNetwork::renderVisNetwork({
      req(!is.null(network_data_deg()))
      
      vis_network_plot(
        network_data = network_data_deg()
      )
    })

    output$show_enrichment <- renderTable({
      if(is.null(enrichment_dataframe())) {
        return(as.data.frame("No significant enrichment found."))
      } 

      res <- enrichment_dataframe()
      colnames(res) <- gsub("\\.", " ", colnames(res))
      if(input$sort_by == "Fold") {
        res <- res[order(res$group, -res$'Fold enriched'), ]
      }
      # if only one group remove group column
      if(length(unique(res$group)) == 1) {
        res <- res[, -1]
      } else {
        # if multiple groups clean up 
        res$group[duplicated(res$group)] <- ""
      }

      res$'nGenes' <- as.character(res$'nGenes')
      res$'Fold enriched' <- as.character(round(res$'Fold enriched', 1))
      res$'Pathway size' <- as.character(
        res$'Pathway size'
      )

      res$'Pathway' <- hyperText(
        res$'Pathway',
        res$URL
      )
      res <- subset(res, select = -Genes)
      res <- subset(res, select = -URL)
      colnames(res)[ncol(res)] <- "Pathway (Click for more info)"
      return(res)

    },
    digits = -1,
    spacing = "s",
    striped = TRUE,
    bordered = TRUE,
    width = "auto",
    hover = TRUE,
    sanitize.text.function = function(x) x
    )

  })

}
    
## To be copied in the UI
# mod_11_enrichment_ui("11_enrichment_1")
    
## To be copied in the server
# mod_11_enrichment_server("11_enrichment_1")
