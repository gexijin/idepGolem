#' 08_pathway UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_06_pathway_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Pathway",
    sidebarLayout(
      sidebarPanel(
        htmlOutput(
          outputId = ns("list_comparisons_pathway")
        ),
        tags$style(
          type = "text/css",
          "#pathway-list_comparisons_pathway{width:100%;   margin-top:-12px}"
        ),
        selectInput(
          inputId = ns("pathway_method"), 
          label = "Select method:", 
          choices = list(
            "GAGE" = 1, 
            "GSEA (preranked fgsea)" = 3,
            "PGSEA" = 2, 
            "PGSEA w/ all samples" = 4, 
            "ReactomePA" = 5
          ),
          selected = 1
        ),
        tags$style(
          type = "text/css",
          "#pathway-pathway_method{width:100%;   margin-top:-12px}"
        ),
        htmlOutput(
          outputId = ns("select_go")
        ),
        tags$style(
          type = "text/css",
          "#pathway-select_go{width:100%;   margin-top:-12px}"
        ),
        fluidRow( 
          column(
            width = 6,
            numericInput(
              inputId = ns("min_set_size"), 
              label = h5("Geneset size: Min."), 
              min = 5, 
              max = 30, 
              value = 15,
              step = 1
            )
          ),
          column(
            width = 6,
            numericInput(
              inputId = ns("max_set_size"), 
              label = h5("Max."), 
              min = 1000, 
              max = 2000, 
              value = 2000,
              step = 100
            )
          )
        ),
        tags$style(
          type = "text/css",
          "#pathway-min_set_size{width:100%;   margin-top:-12px}"
        ),
        tags$style(
          type = "text/css",
          "#pathway-max_set_size{width:100%;   margin-top:-12px}"
        ),
        numericInput(
          inputId = ns("pathway_p_val_cutoff"), 
          label = h5("Pathway signifiance cutoff (FDR)"),
          value = 0.2,
          min = 1e-20,
          max = 1,
          step = .05
        ),
        tags$style(
          type = "text/css",
          "#pathway-pathway_p_val_cutoff{width:100%;   margin-top:-12px}"
        ),
        numericInput(
          inputId = ns("n_pathway_show"), 
          label = h5("Number of top pathways to show"),
          value = 30, 
          min = 5,
          max = 100,
          step = 5
        ),
        tags$style(
          type = "text/css",
          "#pathway-n_pathway_show{width:100%;   margin-top:-12px}"
        ),
        checkboxInput(
          inputId = ns("absolute_fold"),
          label = "Use absolute values of fold changes for GSEA and GAGE",
          value = FALSE
        ),
        numericInput(
          inputId = ns("gene_p_val_cutoff"), 
          label = h5("Remove genes with big FDR before pathway analysis:"),
          value = 1,
          min = 1e-20,
          max = 1,
          step = .05
        ),
        tags$style(
          type = "text/css",
          "#pathway-gene_p_val_cutoff{width:100%;   margin-top:-12px}"
        ),

        # if pathway analysis using methods other than ReactomePA
        # GET RID AND PUT DIFFERENT CHECKS IN
        conditionalPanel(
          condition = "input.pathway_method == 1 | input.pathway_method == 2 |
                       input.pathway_method == 3 | input.pathway_method == 4",
          actionButton("ModalEnrichmentPlotPathway", "Pathway tree"),
          actionButton("ModalVisNetworkPA", "Network(New!)"),
          tags$head(tags$style("#ModalVisNetworkPA{color: red}")),
          downloadButton('downloadPathwayListData', "Pathway list w/ genes"),
          ns = ns          
        ),
        conditionalPanel(
          condition = "input.pathway_method == 4", 
          downloadButton("PGSEAplotAllSamples.Download","High-resolution figure"),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.pathway_method == 2", 
          downloadButton("PGSEAplot.Download", "High-resolution figure"),
          ns = ns
        ),
        h5("* Warning! The many combinations can lead to false positives in pathway analyses."),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/pathways/",
          target = "_blank"
        )    
      ),
      mainPanel(
        conditionalPanel(
          condition = "input.pathway_method == 1",
          tableOutput(outputId = ns("gage_pathway")),
          ns = ns
        ),  
        conditionalPanel(
          condition = "input.pathway_method == 2",
          h5("Red and blue indicates activated and suppressed pathways, respectively."),
          plotOutput(
            outputId = ns("pgsea_plot"),
            inline = TRUE
          ),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.pathway_method == 3",
          tableOutput(outputId = ns("fgsea_pathway")),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.pathwayMethod == 4",
          h5("Red and blue indicates activated and suppressed pathways, respectively."),
          plotOutput(
            outputId = ns("pgsea_plot_all_samples"),
            inline = TRUE
          ),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.pathway_method == 5",
          tableOutput(outputId = ns("reactome_pa_pathway")),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.pathway_method == 1 | input.pathway_method == 2 |
                       input.pathway_method == 3 | input.pathway_method == 4",
          htmlOutput(outputId = ns("list_sig_pathways")),
          downloadButton("downloadSelectedPathwayData", "Expression data for genes in selected pathway"),
          ns = ns
        ),
        conditionalPanel(
          condition = "(input.pathway_method == 1 | input.pathway_method == 2 | 
                        input.pathway_method == 3 | input.pathway_method == 4) &
                        input.select_go == 'KEGG'",
          h5("Red and green represent up- and down-regulated genes, respectively."),
          imageOutput(
            outputId = "kegg_image",
            width = "100%",
            height = "100%"
          ),
          ns = ns
        ),
        conditionalPanel(
          condition = "(input.pathway_method == 1 | input.pathway_method == 2 |
                        input.pathway_method == 3 | input.pathway_method == 4) &
                        input.select_go != 'KEGG'",
          plotOutput(outputId = ns("selected_pathway_heatmap")),
          ns = ns
        ),
        #bsModal(
         # "ModalEnrichmentPlotPathway1",
          #"Significant pathways",
          #"ModalEnrichmentPlotPathway",
          #size="large",
          #h5("Gene sets closer on the tree share more genes. Sizes of dot correspond to adjuested Pvalues"),
          #downloadButton('enrichmentPlotPathway4Download',"Figure"),
          #plotOutput('enrichmentPlotPathway')
        #),
        #bsModal(
          #"ModalEnrichmentPlotPahtway2",
          #"Significant pathways",
          #"ModalEnrichmentNetworkPathway",
          #size="large",
          #h5("Connected gene sets share more genes. Color of node correspond to adjuested Pvalues."),
          #actionButton("layoutButton3", "Change layout"),
          #plotOutput('enrichmentNetworkPlotPathway')
        #),
        # visNetwork ------------------------------
        #bsModal(
          #"ModalVisNetworkPA1",
          #"Related pathways",
          #"ModalVisNetworkPA",
          #size="large",
          #h5("Connected gene sets share more genes. Color of node correspond to adjuested Pvalues."),
          #fluidRow(
            #column(2, actionButton("layoutVisPA", "Change layout") ),
            #column(1, h5("Cutoff:"), align="right" ) ,
            #column(2, numericInput("edgeCutoffPA", label = NULL, value = 0.30, min = 0, max = 1, step = .1), align="left"  ), 
            #column(2, checkboxInput("wrapTextNetworkPA", "Wrap text", value = TRUE)), 
            #column(1, downloadButton("visNetworkPADownload","Network") ),
            #column(1, downloadButton("downloadEdgesPA", "Edges")) ,
            #column(1, downloadButton("downloadNodesPA", "Nodes"))
          #),
          #selectInput(
            #"upORdownRegPA",
            #NULL,
						#c("Both Up & Down" = "Both",
						  #"Up regulated" = "Up",
						  #"Down regulated" = "Down")
          #),
          #h6(
            #"This interactive plot also shows the relationship between enriched pathways. 
			      #Two pathways (nodes) are connected if they share 30% (default, adjustable) or
            #more genes. Green and red represents down- and up-regulated pathways. You can
            #move the nodes by dragging them, zoom in and out by scrolling, and shift the 
            #entire network by click on an empty point and drag. Darker nodes are more 
            #significantly enriched gene sets. Bigger nodes represent larger gene sets.  
			      #Thicker edges represent more overlapped genes."
          #),
          #visNetwork::visNetworkOutput("visNetworkPA",height = "800px", width = "800px")
        #)
      )
    )    
  ) 
}

#' 08_pathway Server Functions
#'
#' @noRd
mod_06_pathway_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    gene_sets <- reactive({
      if(pre_process$select_org() == "NEW" && !is.null(pre_process$gmt_file())) {
        in_file <- pre_process$gmt_file()
        in_file <- in_file$datapath
				return(read_gmt_robust(inFile))
      } else {
        read_gene_sets(
          converted = pre_process$converted(),
          all_gene_names = pre_process$all_gene_names(),
          go = input$select_go,
          select_org = pre_process$select_org(),
          idep_data = idep_data,
          my_range = c(input$min_set_size, input$max_set_size)
        )
      }
    })

    gage_pathway_data <- reactive({

    })
  })
}

## To be copied in the UI
# mod_08_pathway_ui("08_pathway_ui_1")

## To be copied in the server
# mod_08_pathway_server("08_pathway_ui_1")
