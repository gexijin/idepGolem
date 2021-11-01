#' 06_pathway UI Function
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
          "#pathway-list_comparisons_pathway { width:100%;   margin-top:-12px}"
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
          "#pathway-pathway_method { width:100%;   margin-top:-12px}"
        ),
        htmlOutput(
          outputId = ns("select_go_selector")
        ),
        tags$style(
          type = "text/css",
          "#pathway-select_go { width:100%;   margin-top:-12px}"
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
          "#pathway-min_set_size { width:100%;   margin-top:-12px}"
        ),
        tags$style(
          type = "text/css",
          "#pathway-max_set_size { width:100%;   margin-top:-12px}"
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
          "#pathway-pathway_p_val_cutoff { width:100%;   margin-top:-12px}"
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
          "#pathway-n_pathway_show { width:100%;   margin-top:-12px}"
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
          "#pathway-gene_p_val_cutoff { width:100%;   margin-top:-12px}"
        ),
        h5("* Warning! The many combinations can lead to false positives in pathway analyses."),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/pathways/",
          target = "_blank"
        )    
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Results",
            br(),
            conditionalPanel(
              condition = "input.pathway_method == 1",
              DT::dataTableOutput(outputId = ns("gage_pathway_table")),
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
              DT::dataTableOutput(outputId = ns("fgsea_pathway")),
              ns = ns
            ),
            conditionalPanel(
              condition = "input.pathway_method == 4",
              h5("Red and blue indicates activated and suppressed pathways, respectively."),
              plotOutput(
                outputId = ns("pgsea_plot_all_samples"),
                inline = TRUE
              ),
              ns = ns
            ),
            conditionalPanel(
              condition = "input.pathway_method == 5",
              DT::dataTableOutput(outputId = ns("reactome_pa_pathway")),
              ns = ns
            )
          ),
          tabPanel(
            "Heatmap",
            conditionalPanel(
              condition = "input.pathway_method == 1 | input.pathway_method == 2 |
                           input.pathway_method == 3 | input.pathway_method == 4",
              htmlOutput(outputId = ns("list_sig_pathways")),
              ns = ns
            ),
            conditionalPanel(
              condition = "(input.pathway_method == 1 | input.pathway_method == 2 |
                            input.pathway_method == 3 | input.pathway_method == 4) &
                            input.select_go != 'KEGG'",
              plotOutput(outputId = ns("selected_pathway_heatmap")),
              ns = ns
            )
          ),
          tabPanel(
            "KEGG",
            conditionalPanel(
              condition = "(input.pathway_method == 1 | input.pathway_method == 2 | 
                            input.pathway_method == 3 | input.pathway_method == 4) &
                            input.select_go == 'KEGG'",
              h5("Red and green represent up- and down-regulated genes, respectively."),
              imageOutput(
                outputId = ns("kegg_image"),
                width = "100%",
                height = "100%"
              ),
              ns = ns
            )
          )
        )
      )
    )
  )
}

#' mod_06_pathway Server Functions
#'
#' @noRd
mod_06_pathway_server <- function(id, pre_process, deg, idep_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

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

    output$list_comparisons_pathway <- renderUI({
      if(is.null(deg$limma()$comparisons)) {
        selectInput(
          inputId = ns("select_contrast"),
          label = NULL,
          choices = list("All" = "All"),
          selected = "All"
        )  
			}	else {
        selectInput(
          inputId = ns("select_contrast"),
          label = 
            "Select a comparison to examine. \"A-B\" means A vs. B (See heatmap).
            Interaction terms start with \"I:\"",
          choices = deg$limma()$comparisons
	     )
      } 
	  })

    output$list_sig_pathways <- renderUI({
	    if(tab() != "Pathway") {
        selectInput(
          inputId = ns("sig_pathways"),
          label = NULL, 
          choices = list("All" = "All"),
          selected = "All"
        )
      }	else {
        req(!is.null(input$pathway_method))
        # Default, sometimes these methods returns "No significant pathway found"
		    choices <- "All"  
		    if(input$pathway_method == 1) { 
			    if(!is.null(gage_pathway_data())) {
            if(dim(gage_pathway_data())[2] > 1) {
              choices <- gage_pathway_data()[, 2]
            } 
          }
		    } else if(input$pathway_method == 2) {
          if(!is.null(pgsea_plot_data())) {
            if(dim(pgsea_plot_data())[2] > 1) {
					 	  pathways <- as.data.frame(pgsea_plot_data())
						  choices <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
					  }
          }
				} else if(input$pathway_method == 3) {
          if(!is.null(fgsea_pathway_data())) {
            if(dim(fgsea_pathway_data())[2] > 1) {
              choices <- fgsea_pathway_data()[, 2]
            }
          }
        } else if(input$pathway_method == 4) {
          if(!is.null(pgsea_plot_all_samples_data())) {
            if(dim(pgsea_plot_all_samples_data())[2] > 1) {
					    pathways <- as.data.frame(pgsea_plot_all_samples_data())
					    choices <- substr(rownames(pathways), 10, nchar(rownames(pathways)))
				    }
          }
        } else if(input$pathway_method == 5) {
			    if(!is.null(reactome_pa_pathway_data())) {
            if(dim(reactome_pa_pathway_data())[2] > 1) {
              choices <- reactome_pa_pathway_data()[, 2]
            }  
          }
        }
        
        selectInput(
          inputId = ns("sig_pathways"),
          label = 
            "Select a pathway to show expression pattern of related genes on a heatmap or a
             KEGG pathway diagram:",
          choices = choices
        )
	    } 
	  })

    gene_sets <- reactive({
      req(tab() == "Pathway")
      req(!is.null(input$select_go))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Finding Gene Sets",
        color = "#000000"
      )

      if(pre_process$select_org() == "NEW" && !is.null(pre_process$gmt_file())) {
        in_file <- pre_process$gmt_file()
        in_file <- in_file$datapath
				gene_sets <- read_gmt_robust(in_file)
      } else {
        gene_sets <- read_gene_sets(
          converted = pre_process$converted(),
          all_gene_names = pre_process$all_gene_names(),
          go = input$select_go,
          select_org = pre_process$select_org(),
          idep_data = idep_data,
          my_range = c(input$min_set_size, input$max_set_size)
        )
      }

      shinybusy::remove_modal_spinner()

      return(gene_sets)
    })

    gage_pathway_data <- reactive({  
      req(input$pathway_method == 1)
      req(!is.null(deg$limma()))
      req(!is.null(gene_sets()))

      gage_data(
        select_go = input$select_go,
        select_contrast = input$select_contrast,
        min_set_size = input$min_set_size,
        max_set_size = input$max_set_size,
        limma = deg$limma(),
        gene_p_val_cutoff = input$gene_p_val_cutoff,
        gene_sets = gene_sets(),
        absolute_fold = input$absolute_fold,
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show
      )
    })

    output$gage_pathway_table <- DT::renderDataTable({
      req(!is.null(gage_pathway_data()))
      
      DT::datatable(
        gage_pathway_data(),
        options = list(
          pageLength = 15,
          scrollX = "500px"
        ),
        rownames = FALSE
      )
    })

    contrast_samples <- reactive({
      req(!is.null(input$select_contrast))
      req(!is.null(pre_process$data()))

      find_contrast_samples(
        select_contrast = input$select_contrast, 
		    all_sample_names = colnames(pre_process$data()),
		    sample_info = pre_process$sample_info(),
		    select_factors_model = deg$select_factors_model(),
		    select_model_comprions = deg$select_model_comprions(), 
		    reference_levels = deg$factor_reference_levels(),
		    counts_deg_method = deg$counts_deg_method(),
		    data_file_format = pre_process$data_file_format()
	    )
    })

    output$pgsea_plot <- renderPlot({
      req(input$pathway_method == 2)

      plot_pgsea(
        my_range = c(input$min_set_size, input$max_set_size),
        processed_data = pre_process$data(),
        contrast_samples = contrast_samples(),
        gene_sets = gene_sets(),
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show
      )
    }, 
      height = 800,
      width = 800
    )

    pgsea_plot_data <- reactive({
      req(!is.null(gene_sets()))

      get_pgsea_plot_data(
        my_range = c(input$min_set_size, input$max_set_size),
        data = pre_process$data(),
        select_contrast = input$select_contrast,
        gene_sets = gene_sets(),
        sample_info = pre_process$sample_info(),
        select_factors_model = deg$select_factors_model(),
        select_model_comprions = deg$select_model_comprions(),
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show
      )
    })

    fgsea_pathway_data <- reactive({
      req(input$pathway_method == 3)
      req(!is.null(deg$limma()))
      req(!is.null(gene_sets()))

      fgsea_data(
        select_contrast = input$select_contrast,
        my_range = c(input$min_set_size, input$max_set_size),
        limma = deg$limma(),
        gene_p_val_cutoff = input$gene_p_val_cutoff,
        gene_sets = gene_sets(),
        absolute_fold = input$absolute_fold,
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show
      )
    })

    output$fgsea_pathway <- renderTable({
      req(!is.null(fgsea_pathway_data()))
      
      DT::datatable(
        fgsea_pathway_data(),
        options = list(
          pageLength = 15,
          scrollX = "500px"
        ),
        rownames = FALSE
      )
	  })

    reactome_pa_pathway_data <- reactive({
      req(!is.null(deg$limma()))

      reactome_data(
        select_contrast = input$select_contrast,
        my_range = c(input$min_set_size, input$max_set_size),
        limma = deg$limma(),
        gene_p_val_cutoff = input$gene_p_val_cutoff,
        converted = pre_process$converted(),
        idep_data = idep_data,
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show,
        absolute_fold = input$absolute_fold
      )
    })

    output$pgsea_plot_all_samples <- renderPlot({
      req(input$pathway_method == 4)

      pgsea_plot_all(
        go = input$select_go,
        my_range = c(input$min_set_size, input$max_set_size),
        data = pre_process$data(),
        select_contrast = input$select_contrast,
        gene_sets = gene_sets(),
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show
      )
    }, height = 800, width = 800)

    pgsea_plot_all_samples_data <- reactive({
      req(!is.null(gene_sets()))

      get_pgsea_plot_all_samples_data(
        data = pre_process$data(),
        select_contrast = input$select_contrast,
        gene_sets = gene_sets(),
        my_range = c(input$min_set_size, input$max_set_size),
        pathway_p_val_cutoff = input$pathway_p_val_cutoff,
        n_pathway_show = input$n_pathway_show
      )
    })




    output$reactome_pa_pathway <- renderTable({
      req(!is.null(reactome_pa_pathway_data()))
	    
      DT::datatable(
        reactome_pa_pathway_data(),
        options = list(
          pageLength = 15,
          scrollX = "500px"
        ),
        rownames = FALSE
      )
	  })
  })
}

## To be copied in the UI
# mod_08_pathway_ui("08_pathway_ui_1")

## To be copied in the server
# mod_08_pathway_server("08_pathway_ui_1")
