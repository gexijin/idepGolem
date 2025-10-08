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
        ),
        tippy::tippy_this(
          ns("customize_button"),
          "Reveal additional settings for sorting and display.",
          theme = "light"
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
          selected = "FDR",
          selectize = FALSE
        ),
        tippy::tippy_this(
          ns("sort_by"),
          "Order pathways by FDR or fold enrichment.",
          theme = "light"
        )
      ),
      column(
        width = 4,
        checkboxInput(
          inputId = ns("filtered_background"),
          label = "Use filtered genes as background.",
          value = TRUE
        ),
        tippy::tippy_this(
          ns("filtered_background"),
          "Use only filtered genes as the background universe for enrichment.",
          theme = "light"
        )
      ),
      column(
        width = 4,
        checkboxInput(
          inputId = ns("remove_redudant"),
          label = "Remove Redudant Gene Sets",
          value = FALSE
        ),
        tippy::tippy_this(
          ns("remove_redudant"),
          "Collapse overlapping pathways to reduce redundancy.",
          theme = "light"
        )
      )
    ),
    fluidRow(
      column(
        width = 4,
        align = "left",
        # not that pathways are still ranked by FDR
        numericInput(
          inputId = ns("top_pathways"),
          label = "Top pathways",
          min = 1,
          max = 30,
          value = 10
        ),
        tippy::tippy_this(
          ns("top_pathways"),
          "Choose how many pathways to display across the outputs.",
          theme = "light"
        )
      ),
      column(
        width = 4,
        align = "left",
        checkboxInput(
          inputId = ns("show_pathway_id"),
          label = "Show pathway IDs",
          value = FALSE
        ),
        tippy::tippy_this(
          ns("show_pathway_id"),
          "Append pathway IDs (e.g., Path:mmu04115 or GO:0042770) to each pathway name.",
          theme = "light"
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
          "Download the enrichment results table.",
          theme = "light"
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
            tippy::tippy_this(
              ns("edge_cutoff_deg"),
              "Filter network edges by overlap strength.",
              theme = "light"
            ),
            align = "left"
          ),
          column(
            width = 2,
            actionButton(
              inputId = ns("layout_vis_deg"),
              label = "Change layout"
            ),
            tippy::tippy_this(
              ns("layout_vis_deg"),
              "Rearrange the enrichment network layout.",
              theme = "light"
            )
          ),
          column(
            width = 2,
            checkboxInput(
              inputId = ns("wrap_text_network_deg"),
              label = "Wrap text",
              value = TRUE
            ),
            tippy::tippy_this(
              ns("wrap_text_network_deg"),
              "Wrap long pathway names on the network nodes.",
              theme = "light"
            )
          )
        ),
        visNetwork::visNetworkOutput(
          outputId = ns("vis_network_deg"),
          height = "800px",
          width = "100%"
        ),
        conditionalPanel(
          condition = "input.select_cluster == 'All Groups'",
          h6(
            "Two pathways (nodes) are connected if they
            share 30% (default, adjustable) or more genes.
            Green and red represents up- and down-regulated
            pathways. You can move the nodes by
            dragging them, zoom in and out by scrolling, and
            shift the entire network by click on an
            empty point and drag. Darker nodes are more
            significantly enriched gene sets. Bigger nodes
            represent larger gene sets. Thicker edges
            represent more overlapped genes."
          ),
          ns = ns
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
              selected = column_selection[1],
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("pathway_order"),
              "Choose how pathways are ordered in the bar chart.",
              theme = "light"
            )
          ),
          column(
            width = 3,
            selectInput(
              inputId = ns("order_x"),
              label = h5("x-axis"),
              choices = column_selection[1:3],
              selected = column_selection[1],
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("order_x"),
              "Select the metric displayed on the x-axis.",
              theme = "light"
            )
          ),
          column(
            width = 3,
            selectInput(
              inputId = ns("plot_color"),
              label = h5("Color"),
              choices = column_selection[1:3],
              selected = column_selection[2],
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("plot_color"),
              "Choose which metric controls point color.",
              theme = "light"
            )
          ),
          column(
            width = 3,
            selectInput(
              inputId = ns("plot_size"),
              label = h5("Size"),
              choices = column_selection[1:3],
              selected = column_selection[3],
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("plot_size"),
              "Choose which metric controls point size.",
              theme = "light"
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
            ),
            tippy::tippy_this(
              ns("font_size"),
              "Adjust label text size on the chart.",
              theme = "light"
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
            ),
            tippy::tippy_this(
              ns("marker_size"),
              "Control the marker size for each pathway.",
              theme = "light"
            )
          ),
          column(
            width = 3,
            selectInput(
              inputId = ns("high_color"),
              label = h5("Color:High"),
              choices = c("red", "orange", "yellow", "green", "blue", "purple"),
              selected = "red",
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("high_color"),
              "Set the color for high values.",
              theme = "light"
            )
          ),
          column(
            width = 3,
            selectInput(
              inputId = ns("log_color"),
              label = h5("Color:Low"),
              choices = c("red", "orange", "yellow", "green", "blue", "purple"),
              selected = "blue",
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("log_color"),
              "Set the color for low values.",
              theme = "light"
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
              selected = "lollipop",
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("chart_type"),
              "Switch between lollipop, dot, or bar chart styles.",
              theme = "light"
            )
          ),
          column(
            width = 3,
            selectInput(
              inputId = ns("aspect_ratio"),
              label = h5("Aspect Ratio"),
              choices = .1 * (5:30),
              selected = 2,
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("aspect_ratio"),
              "Adjust the height-to-width ratio of the chart.",
              theme = "light"
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
            "Download the enrichment results table.",
            theme = "light"
          ),
          column(
            width = 3,
            checkboxInput(
              ns("show_detail"),
              "Detailed Desc.",
              value = FALSE
            ),
            tippy::tippy_this(
              ns("show_detail"),
              "Toggle detailed descriptions for each pathway.",
              theme = "light"
            )
          )
        ),
        tableOutput(ns("gene_info_table")),
        p("Note: In the gene type column, \"C\" indicates
        protein-coding genes, and \"P\" means pseduogenes.")
      ),
      tabPanel(
        title = "Details",
        value = "Details",
        h4("Details of Enrichment Analysis"),
        htmlOutput(ns("enrichment_details"))
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
                                     filter_size,
                                     gene_info,
                                     idep_data,
                                     select_org,
                                     converted,
                                     gmt_file,
                                     plot_grid_lines,
                                     ggplot2_theme,
                                     heat_colors) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    observe({
      shinyjs::toggle(id = "sort_by", condition = input$customize_button)
      shinyjs::toggle(id = "filtered_background", condition = input$customize_button)
      shinyjs::toggle(id = "remove_redudant", condition = input$customize_button)
      shinyjs::toggle(id = "top_pathways", condition = input$customize_button)
      shinyjs::toggle(id = "show_pathway_id", condition = input$customize_button)
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
        selected = selected,
        selectize = FALSE
      )
    })

    observe({
      req(input$subtab == "Plot")
      choices <- sort(unique(enrichment_dataframe()$group))
      selected <- choices[1]
      updateSelectInput(
        session = session,
        inputId = ns("select_cluster"),
        selected = selected
      )
    })
    
    observe({
      req(!is.null(filter_size()))
      
      if(filter_size() < 1000) {
        
        updateCheckboxInput(
          session = session,
          inputId = "filtered_background",
          value = FALSE
        )
      } else {
        updateCheckboxInput(
          session = session,
          inputId = "filtered_background",
          value = TRUE
        )
      }
      
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
      choices_name <- choices
      ix <- which(nchar(choices) == 1) # "1" <- "Cluster 1"
      choices_name[ix] <- paste("Cluster ", choices[ix])
      if (min(nchar(choices)) == 1) {
        choices <- setNames(choices, choices_name)
      }

      selectInput(
        inputId = ns("select_cluster"),
        label = NULL,
        choices = choices,
        selected = selected,
        selectize = FALSE
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
        incProgress(0.2)
        pathway_info <- list()
        # disregard user selection use clusters for enrichment
        for (i in 1:length(gene_lists())) {
          incProgress(0.2 + length(gene_lists()) / 20)
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
            # when GMT file is uploaded, pretend to be a new species
            select_org = ifelse(is.null(gmt_file()), select_org(), "NEW"),
            use_filtered_background = input$filtered_background,
            reduced = input$remove_redudant,
            max_terms = input$top_pathways,
            sort_by_fold = (input$sort_by == "Fold")
          )
        }
      })
      # remove pathway ID for GO and KEGG
      # Path:hsa00270 Cysteine and methionine metabolism 
      #           --> Cysteine and methionine metabolism
                                    # result is not NULL
      # pathway_info is a list of data frames: Selection, Upregulated, Downregulated, etc
      if (!input$show_pathway_id && select_org() > 0) {
        for (i in 1:length(pathway_info)) {
          pathway_info[[i]]$Pathway <- remove_pathway_id(
            strings = pathway_info[[i]]$Pathway,
            select_go = input$select_go
          )
        }
      }
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
            if (nrow(df1) == 0) {
              return(NULL)
            }
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
          results_all <- results_all[,
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
          results_all <- results_all[,
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
        right_margin = 30,
        leaf_color_choices = heat_colors()
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
        edge_cutoff_deg = input$edge_cutoff_deg,
        group_color = heat_colors()
      )
    })

    # Interactive vis network plot
    output$vis_network_deg <- visNetwork::renderVisNetwork({
      req(!is.null(network_data_deg()))
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
        #        if (input$sort_by == "Fold") {
        #          res <- res[order(res$group, -res$"Fold enriched"), ]
        #        }
        # if only one group remove group column
        if (length(unique(res$group)) == 1) {
          res <- res[, -1]
        } else {
          # if multiple groups clean up
          res$group[duplicated(res$group)] <- ""
        }
        # 2.1e-03  --> 2.1e-3
        res$FDR <- gsub("e-0", "e-", res$FDR)
        res$FDR <- gsub("e", "E", res$FDR)
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
        colnames(res) <- gsub("Fold enriched", "Fold", colnames(res))
        colnames(res) <- gsub("FDR", "Adj.Pval", colnames(res))
        colnames(res) <- gsub("group", "Grp.", colnames(res))
        res <- subset(res, select = -PathwaySize)

        if (input$select_go != "All"){
          colnames(res)[ncol(res)] <- "Pathway (Click for more info)"
        } else {
          colnames(res)[ncol(res) - 1]
        }
        res <- subset(res, select = -nGenes)
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

      p <- enrich_barplot(
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
      refine_ggplot2(
        p = p,
        gridline = plot_grid_lines(),
        ggplot2_theme = ggplot2_theme()
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

    output$enrichment_details <- renderUI ({

      info <- pathway_source_info(
        pathway_file = gmt_file(),
        go = input$select_go,
        select_org = select_org(),
        idep_data = idep_data
      )

      tagList(
        p("Enrichment P-value are calculated based on one-sided the hypergeometric test, which is then 
                adjusted for multiple testing using the Benjamini-Hochberg procedure and converted to 
                FDR(false discovery rate). FDR tells us how likely to observe the enrichment by chance. Due to increased statistical power,
                large pathways tend to have smaller FDRs.
                Fold Enrichment is defined as the percentage
                of genes in the list belonging to a pathway, divided by the corresponding percentage in the
                background. As a measure of effect size, Fold Enrichment indicates 
                how drastically genes of a certain pathway is overrepresented.
                "),
        if(input$filtered_background) {
          p("The background genes are filtered genes from the original gene list. 
          The filtered genes are those that are passed a low filter in RNA-seq.")
        } else {
          p("The background genes are all protein-coding genes.")
        },

        if(input$remove_redudant) {
          p("Similar pathways sharing 90% of genes are represented by the most significant pathway if
                they also share 50% of the words in their names.")
        },


        if(!is.null(info)) {
          p("A total of ", info$GeneSets, " gene sets are obtained from ", 
            a(
              info$Subtype.Database.name,
              href = info$Link
            ),
            "(", info$Database.Description, ").",
            "Dabase version/date is ", info$Version, ".",
            "For details on  this ", info$Type, "type of database, please refer to: ",
            info$Author, ", ", 
            info$PaperTitle, ", ",
            info$Citation.Reference, 
            "(",
            a(
              paste("PubMed ID: ", info$PMID),
              href = paste0("https://www.ncbi.nlm.nih.gov/pubmed/", info$PMID)
            ),
            ")."
          )#p
        } else {
          p("Gene sets are derived from ", paste0(input$select_go, "."))
        },
        p("After the analysis is done, pathways are first filtered based on a FDR cutoff (0.05).
          Then the siginificant pathways are sorted by.", input$sort_by, ".",
          " Only the top ", input$top_pathways, " pathways are shown.")
      )
    })    
  })


}

## To be copied in the UI
# mod_11_enrichment_ui("11_enrichment_1")

## To be copied in the server
# mod_11_enrichment_server("11_enrichment_1")
