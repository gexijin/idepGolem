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
    tableOutput(ns("enrichment_table")),
    downloadButton(
      outputId = ns("download_enrichment"),
      label = "Enrichment"
    )
  )
}
    
#' 11_enrichment Server Functions
#' @param id module Id
#' @param results a list containing results from pathway analysis.
#'  Each element is a data frame. But the last column is a list,
#' hosting genes.
#' @noRd 
mod_11_enrichment_server <- function(id, results){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

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

      # if only one group remove group column
      if(length(unique(res$group)) == 1) {
        res <- res[, -1]
      } else {
        # if multiple groups clean up 
        res$group[duplicated(res$group)] <- ""
      }
      colnames(res) <- gsub("\\.", " ", colnames(res))
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
