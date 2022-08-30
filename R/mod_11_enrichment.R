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
  tagList(

    fluidRow(
      column(
        width = 4,
        htmlOutput(outputId = ns("select_go_selector"))
      ),
      column(
        width = 4,
        checkboxInput(
          inputId = ns("filtered_background"),
          label = "Use filtered genes as background.",
          value = FALSE
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
    fluidRow(
      column(
        width = 3, 
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
        width = 2,
        downloadButton(
          outputId = ns("download_enrichment"),
          label = "Enrichment"
        )
      )
    ),
    tableOutput(ns("enrichment_table"))
  )
}
    
#' 11_enrichment Server Functions
#' @param id module Id
#' @param results a list containing results from pathway analysis.
#'  Each element is a data frame. But the last column is a list,
#' hosting genes.
#' @noRd 
mod_11_enrichment_server <- function(id, results, gmt_choices){
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

    full_table <- reactive({
      req(!is.null(results()))

      results_all <- do.call(rbind,
        #combine multiple data frames that are elements of a list
        lapply(
          names(results()),
          function(x) {
            if(ncol(results()[[x]]) == 1) {
              return(NULL)
            }
            df1 <- data_frame_with_list(results()[[x]])
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
    })

    output$enrichment_table <- renderTable({
      if(is.null(full_table())) {
        return(as.data.frame("No significant enrichment found."))
      } 

      res <- full_table()
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

    output$download_enrichment <- downloadHandler(
      filename = function() {
        "enrichment.csv"
      },
      content = function(file) {
        write.csv(full_table(), file)
      }
    )


  })

}
    
## To be copied in the UI
# mod_11_enrichment_ui("11_enrichment_1")
    
## To be copied in the server
# mod_11_enrichment_server("11_enrichment_1")
