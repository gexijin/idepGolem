# Global variable: list of column names in the pathway enrichment
# results data frame, and their display name on the plot.
column_selection <- list(
  "-log10(FDR)" = "EnrichmentFDR",
  "Fold Enrichment" = "FoldEnrichment",
  "Genes" = "nGenes",
  "Category Name" = "Pathway"
)

#' Shows results from enrichment analysis of one or more lists of genes
#'
#' @description The input is a list of genes
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_11_enrichment_ui <- function(id) {
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
        htmlOutput(outputId = ns("select_cluster"))
      ),
      column(
        width = 4,
        checkboxInput(
          inputId = ns("customize_button"),
          label = "More options",
          value = FALSE
        )
      )
    ),
    fluidRow(
      column(
        width = 4,
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
        width = 4,
        checkboxInput(
          inputId = ns("filtered_background"),
          label = "Use filtered genes as background.",
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
    tabsetPanel(
      id = ns("subtab"),
      tabPanel(
        title = "Pathways",
        value = "Pathways",
        tableOutput(ns("show_enrichment")),
        downloadButton(
          outputId = ns("download_enrichment"),
          label = "CSV file"
        ),
        tippy::tippy_this(
          ns("download_enrichment"),
          "Download enrichment analysis",
          theme = "light-border"
        ),
      ),
      tabPanel(
        title = "Tree",
        plotOutput(ns("enrichment_tree")),
        br(),
        ottoPlots::mod_download_figure_ui(
          id = ns("dl_treeplot")
        )
      ),
      tabPanel(
        title = "Network",
        value = "Network",
        br(),
        fluidRow(
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
          "Two pathways (nodes) are connected if they
          share 30% (default, adjustable) or more genes.
          Green and red represents down- and up-regulated
           pathways. You can move the nodes by
          dragging them, zoom in and out by scrolling, and
          shift the entire network by click on an
          empty point and drag. Darker nodes are more
          significantly enriched gene sets. Bigger nodes
          represent larger gene sets. Thicker edges
          represent more overlapped genes."
        )
      ),
      tabPanel(
        title = "Plot",
        value = "Plot",
        plotOutput(ns("enrich_barchart"), width = "100%", height = "100%"),
        fluidRow(
          column(
            width = 3,
            selectInput(
              inputId = ns("pathway_order"),
              label = h5("Sort Pathway by"),
              choices = column_selection,
              selected = column_selection[1]
            )
          ),
          column(
            width = 3,
            selectInput(
              inputId = ns("order_x"),
              label = h5("x-axis"),
              choices = column_selection[1:3],
              selected = column_selection[1]
            )
          ),
          column(
            width = 3,
            selectInput(
              inputId = ns("plot_color"),
              label = h5("Color"),
              choices = column_selection[1:3],
              selected = column_selection[2]
            )
          ),
          column(
            width = 3,
            selectInput(
              inputId = ns("plot_size"),
              label = h5("Size"),
              choices = column_selection[1:3],
              selected = column_selection[3]
            )
          )
        ),
        fluidRow(
          column(
            width = 3,
            numericInput(
              inputId = ns("font_size"),
              label = h5("Font Size"),
              value = 12,
              min = 3,
              max = 18,
              step = 1
            )
          ),
          column(
            width = 3,
            numericInput(
              inputId = ns("marker_size"),
              label = h5("Circle Size"),
              value = 4,
              min = 0,
              max = 10,
              step = 1
            )
          ),
          column(
            width = 3,
            selectInput(
              inputId = ns("high_color"),
              label = h5("Color:High"),
              choices = c("red", "orange", "yellow", "green", "blue", "purple"),
              selected = "red"
            )
          ),
          column(
            width = 3,
            selectInput(
              inputId = ns("log_color"),
              label = h5("Color:Low"),
              choices = c("red", "orange", "yellow", "green", "blue", "purple"),
              selected = "blue"
            )
          )
        ),
        fluidRow(
          column(
            width = 3,
            selectInput(
              inputId = ns("chart_type"),
              label = h5("Chart type"),
              choices = c("lollipop", "dotplot", "barplot"),
              selected = "lollipop"
            )
          ),
          column(
            width = 3,
            selectInput(
              inputId = ns("aspect_ratio"),
              label = h5("Aspect Ratio"),
              choices = .1 * (5:30),
              selected = 2
            )
          ),
          column(
            width = 3,
            style = "margin-top: 30px;",
            ottoPlots::mod_download_figure_ui(
              id = ns("dl_barplot")
            )
          )
        ) # 3rd row
      ),
      tabPanel(
        title = "Genes",
        value = "Genes",
        fluidRow(
          column(
            width = 7,
            textOutput(ns("gene_counts"))
          ),
          column(
            width = 2,
            downloadButton(ns("download_gene_info"), "CSV file")
          ),
          tippy::tippy_this(
            ns("download_gene_info"),
            "Download enrichment analysis results",
            theme = "light-border"
          ),
          column(
            width = 3,
            checkboxInput(
              ns("show_detail"),
              "Detailed Desc.",
              value = FALSE
            )
          )
        ),
        tableOutput(ns("gene_info_table")),
        p("Note: In the gene type column, \"C\" indicates
        protein-coding genes, and \"P\" means pseduogenes.")
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
mod_11_enrichment_server <- function(id,
                                     results, # pathway table with fdr, fold, etc
                                     gmt_choices, # list of pathway categories "GOBP"
                                     gene_lists, # list of genes, each element is a list
                                     processed_data,
                                     gene_info,
                                     idep_data,
                                     select_org,
                                     converted,
                                     gmt_file) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    observe({
      shinyjs::toggle(id = "sort_by", condition = input$customize_button)
      shinyjs::toggle(id = "filtered_background", condition = input$customize_button)
      shinyjs::toggle(id = "remove_redudant", condition = input$customize_button)
    })
    # GMT choices for enrichment ----------
    output$select_go_selector <- renderUI({
      req(!is.null(gmt_choices()))

      selected <- gmt_choices()[1] # default, overwrite by below
      if ("GOBP" %in% gmt_choices()) {
        selected <- "GOBP"
      }
      if ("KEGG" %in% gmt_choices()) {
        selected <- "KEGG"
      }
      selectInput(
        inputId = ns("select_go"),
        label = NULL,
        choices = gmt_choices(),
        selected = selected
      )
    })

    observe({
      req(input$subtab == "Plot")
      choices <- sort(unique(enrichment_dataframe()$group))
      selected <- choices[1]
      updateSelectInput(
        session = session,
        inputId = "select_cluster",
        selected = selected
      )
    })

    output$select_cluster <- renderUI({
      req(!is.null(enrichment_dataframe()))
      choices <- sort(unique(enrichment_dataframe()$group))
      selected <- choices[1]
      if (length(choices) > 1) {
        choices <- c("All Groups", choices)
        # if two groups, defaults to both
        if (length(choices) == 3) {
          selected <- choices[1]
        }
      }
      selectInput(
        inputId = ns("select_cluster"),
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
      withProgress(message = "Enrichment Analysis", {
        pathway_info <- list()
        # disregard user selection use clusters for enrichment
        for (i in 1:length(gene_lists())) {
          incProgress(1 / length(gene_lists()))
          gene_names_query <- gene_lists()[[i]]
          req(!is.null(input$select_go))
          gene_sets <- read_pathway_sets(
            all_gene_names_query = gene_names_query,
            converted = converted(), # n
            go = input$select_go,
            select_org = select_org(),
            gmt_file = gmt_file(), # n
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
      })
      return(pathway_info)
    })

    # a dataframe with information on the genes for enrichment analysis
    cluster_gene_info <- reactive({
      req(!is.null(gene_lists()))
      req(gene_info())
      req(!is.null(input$select_cluster))

      df <- do.call(
        rbind,
        # combine multiple data frames that are elements of a list
        lapply(
          names(gene_lists()),
          function(x) {
            if (input$select_cluster != "All Groups" &&
              input$select_cluster != x
            ) {
              return(NULL)
            }
            ix <- which(gene_info()$ensembl_gene_id %in%
              gene_lists()[[x]]$ensembl_ID)
            df1 <- gene_info()[ix, ]
            df1$group <- x
            return(df1)
          }
        )
      )
      if (!is.null(df)) {
        df <- subset(
          df,
          select = c(
            group,
            ensembl_gene_id,
            symbol,
            entrezgene_id,
            chromosome_name,
            # start_position,
            gene_biotype,
            description
          )
        )
        # show detailed gene info for string species
        if (!input$show_detail) {
          df$description <- gsub(";.*|\\[.*", "", df$description)
        }


        # df$start_position <- round(df$start_position / 1e6, 2)
        # protein_coding --> coding; processed_pseduogene --> pseduogene

        df$gene_biotype <- gsub(".*pseudogene", "P", df$gene_biotype)
        # coding is not shown
        df$gene_biotype <- gsub("coding|protein_coding", "C", df$gene_biotype)
        # TR_J_gene  --> TR_J
        df$gene_biotype <- gsub("_gene", "", df$gene_biotype)

        # GL456211.1 ---> ""
        df$chromosome_name[nchar(df$chromosome_name) > 5] <- ""

        colnames(df) <- c(
          "Group",
          "Ensembl id",
          "Symbol",
          "Entrezgene id",
          "Chr",
          # "Position (Mbp)",
          "Type",
          "Description"
        )
        # Remove columns with all missing values;
        df <- df[, which(!apply(is.na(df), 2, sum) == nrow(df))]
      }

      return(df)
    })

    output$gene_info_table <- renderTable(
      {
        req(!is.null(cluster_gene_info()))
        df <- cluster_gene_info()
        if (length(unique(cluster_gene_info()$Group)) == 1) {
          df <- subset(df, select = -Group)
        } else {
          df$Group[duplicated(df$Group)] <- ""
        }
        # add link to Ensembl
        ix <- grepl("ENS", df$"Ensembl id")
        if (sum(ix) > 0) { # at least one has http?
          tem <- paste0(
            "<a href='http://www.ensembl.org/id/",
            df$"Ensembl id",
            "' target='_blank'>",
            df$"Ensembl id",
            "</a>"
          )
          # only change the ones with URL
          df$"Ensembl id"[ix] <- tem[ix]
        }

        # first see if it is Entrez gene ID-----------------------
        ix <- !is.na(as.numeric(df$"Entrezgene id"))
        if (sum(ix) > 0) { # at least one has http?
          tem <- paste0(
            "<a href='https://www.ncbi.nlm.nih.gov/gene/",
            df$"Entrezgene id",
            "' target='_blank'>",
            df$"Entrezgene id",
            "</a>"
          )
          # only change the ones with URL
          df$"Entrezgene id"[ix] <- tem[ix]
        }
        return(df)
      },
      digits = 2,
      spacing = "s",
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = TRUE,
      sanitize.text.function = function(x) x
    )

    output$download_gene_info <- downloadHandler(
      filename = function() {
        paste(input$select_cluster, "_geneInfo.csv")
      },
      content = function(file) {
        write.csv(
          cluster_gene_info(),
          file,
          row.names = FALSE
        )
      }
    )

    output$gene_counts <- renderText({
      req(!is.null(cluster_gene_info()))
      counts <- table(cluster_gene_info()$Group)
      txt <- paste0(names(counts), ":", counts, collapse = " genes; ")
    })


    # returns a data frame
    enrichment_dataframe <- reactive({
      req(!is.null(pathway_table()))

      results_all <- do.call(
        rbind,
        # combine multiple data frames that are elements of a list
        lapply(
          names(pathway_table()),
          function(x) {
            if (ncol(pathway_table()[[x]]) == 1) {
              return(NULL)
            }
            df1 <- data_frame_with_list(pathway_table()[[x]])
            df1$group <- x
            return(df1)
          }
        )
      )

      if (!is.null(results_all)) {
        if (ncol(results_all) > 1) {
          results_all <- results_all[
            ,
            c(
              "group",
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

      results_all <- do.call(
        rbind,
        # combine multiple data frames that are elements of a list
        lapply(
          names(pathway_table()),
          function(x) {
            if (ncol(pathway_table()[[x]]) == 1) {
              return(NULL)
            }
            df1 <- pathway_table()[[x]]
            df1$group <- x
            return(df1)
          }
        )
      )

      if (!is.null(results_all)) {
        if (ncol(results_all) > 1) {
          results_all <- results_all[
            ,
            c(
              "group",
              colnames(results_all)[1:(ncol(results_all) - 1)]
            )
          ]
        }

        results_all <- subset(
          results_all,
          select = c(group, FDR, nGenes, Pathway, Genes)
        )
        results_all$FDR <- as.numeric(results_all$FDR)
        colnames(results_all) <- c(
          "Direction", "adj_p_val", "Pathway.size", "Pathways", "Genes"
        )
      }
      return(results_all)
    })

    # Enrichment Tree -----------
    enrichment_tree_object <- reactive({
      req(!is.null(enrichment_dataframe_for_tree()))
      req(!is.null(input$select_cluster))
      enrichment_tree_plot(
        go_table = enrichment_dataframe_for_tree(),
        group = input$select_cluster,
        right_margin = 45
      )
    })

    # Enrichment Tree -----------
    output$enrichment_tree <- renderPlot({
      req(!is.null(enrichment_tree_object()))
      enrichment_tree_object()
    })

    dl_treeplot <- ottoPlots::mod_download_figure_server(
      id = "dl_treeplot",
      filename = "enrichment_tree",
      figure = reactive({
        enrichment_tree_object()
      }),
      width = 12,
      height = 6,
      label = ""
    )

    # Define a Network
    network_data_deg <- reactive({
      req(!is.null(enrichment_dataframe_for_tree()))
      req(!is.null(input$select_cluster))

      network_data(
        network = enrichment_dataframe_for_tree(),
        up_down_reg_deg = input$select_cluster,
        wrap_text_network_deg = input$wrap_text_network_deg,
        layout_vis_deg = input$layout_vis_deg,
        edge_cutoff_deg = input$edge_cutoff_deg
      )
    })

    # Interactive vis network plot
    output$vis_network_deg <- visNetwork::renderVisNetwork({
      req(!is.null(network_data_deg()))
      req(nrow(network_data_deg()$edges) > 2)
      req(nrow(network_data_deg()$nodes) > 2)
      vis_network_plot(
        network_data = network_data_deg()
      )
    })

    output$show_enrichment <- renderTable(
      {
        if (is.null(enrichment_dataframe())) {
          return(as.data.frame("No significant enrichment found."))
        }

        res <- enrichment_dataframe()
        colnames(res) <- gsub("\\.", " ", colnames(res))
        if (input$sort_by == "Fold") {
          res <- res[order(res$group, -res$"Fold enriched"), ]
        }
        # if only one group remove group column
        if (length(unique(res$group)) == 1) {
          res <- res[, -1]
        } else {
          # if multiple groups clean up
          res$group[duplicated(res$group)] <- ""
        }

        res$"nGenes" <- as.character(res$"nGenes")
        res$"Fold enriched" <- as.character(round(res$"Fold enriched", 1))
        res$"Pathway size" <- as.character(
          res$"Pathway size"
        )

        res$"Pathway" <- hyperText(
          res$"Pathway",
          res$URL
        )
        res <- subset(res, select = -Genes)
        res <- subset(res, select = -URL)

        # remove pathway size
        colnames(res) <- gsub("Pathway size", "PathwaySize", colnames(res))
        res <- subset(res, select = -PathwaySize)

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

    # ggplot2 object for the enrichment chart;
    # used both for display and download
    enrich_barplot_object <- reactive({
      if (is.null(enrichment_dataframe())) {
        grid::grid.newpage()
        return(grid::grid.text("Not available.", 0.5, 0.5))
      }
      req(input$pathway_order)
      req(input$order_x)
      req(input$plot_size)
      req(input$plot_color)
      req(input$font_size)
      req(input$marker_size)
      req(input$high_color)
      req(input$log_color)
      req(input$chart_type)
      req(input$aspect_ratio)
      req(input$select_cluster)

      enrich_barplot(
        enrichment_dataframe = enrichment_dataframe(),
        pathway_order = input$pathway_order,
        order_x = input$order_x,
        plot_size = input$plot_size,
        plot_color = input$plot_color,
        plot_font_size = input$font_size,
        plot_marker_size = input$marker_size,
        plot_high_color = input$high_color,
        plot_low_color = input$log_color,
        chart_type = input$chart_type,
        aspect_ratio = input$aspect_ratio,
        select_cluster = input$select_cluster
      )
    })

    # Enrichment plot for display on the screen
    # https://stackoverflow.com/questions/34792998/shiny-variable-height-of-renderplot
    output$enrich_barchart <- renderPlot(
      {
        enrich_barplot_object()
      },
      # height increases as the number of terms increase. max at 1200, min 350
      height = function() {
        round(max(350, min(2500, round(18 * as.numeric(20))))) # 20 is maxTerms
      },
      width = function() {
        round(max(350, min(2500, round(18 * as.numeric(20)))) * as.numeric(input$aspect_ratio))
      }
    )

    dl_barplot <- ottoPlots::mod_download_figure_server(
      id = "dl_barplot",
      filename = "enrichment_barplot",
      figure = reactive({
        enrich_barplot_object()
      }),
      width = 8,
      height = 6,
      label = "Download"
    )
  })
}

## To be copied in the UI
# mod_11_enrichment_ui("11_enrichment_1")

## To be copied in the server
# mod_11_enrichment_server("11_enrichment_1")
