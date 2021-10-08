#' 06_deg1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_06_deg1_ui <- function(id) {
  ns <- NS(id)
  navbarMenu(
    "DEG",
    tabPanel(
      title = "Design (Step 1)",
      sidebarLayout(
        sidebarPanel(
          conditionalPanel(
            condition = "output.data_file_format == 1",
            selectInput(
              inputId = ns("counts_deg_method"),
              label = "Method:", 
              choices = list(
                "DESeq2" = 3,
                "limma-voom" = 2,
                "limma-trend" = 1
              ),
              selected = 3
            ),
            tags$style(
              type = 'text/css',
              "#deg-counts_deg_method {width:100%;   margin-top:-12px}"
            )
          ),
          conditionalPanel(
            condition = "output.data_file_format == 2",
            h5("Using the limma package")
          ),
          fluidRow(
            column(
              width = 5,
              numericInput(
                inputId = "limma_pval",
                label = h5("FDR cutoff"), 
                value = 0.1,
                min = 1e-5,
                max = 1,
                step = .05
              )
            ),
            column(
              width = 7,
              numericInput(
                inputId = "limma_fc",
                label = h5("Min fold change"),
                value = 2,
                min = 1,
                max = 100,
                step = 0.5
              )
            ),
            tags$style(
              type = "text/css",
              "#deg-limma_pval { width:100%;   margin-top:-12px}"
            ),
            tags$style(
              type = "text/css",
              "#deg-limma_fc { width:100%;   margin-top:-12px}"
            )
          ),
          a(
            h5("Questions?", align = "right"),
            href = "https://idepsite.wordpress.com/degs/",
            target = "_blank"
          )
        ),
        mainPanel(
          tabsetPanel(
            tabPanel(
              title = "Overview",
              plotOutput(
                outputId = ns("sig_gene_stats")
              ),
              br(),
              br(),
              h4(
                "Numbers of differentially expressed genes for all comparisons.
                \"B-A\" means B vs. A. Interaction terms start with \"I:\" "
              ),
              tableOutput(
                outputId = ns("sig_gene_stats_table")
              )
            ),
            tabPanel(
              title = "Experiment Design",
              fluidRow(
                column(
                  width = 6,
                  htmlOutput(outputId = ns("list_factors_de"))
                ),
                column(
                  width = 6,
                  htmlOutput(outputId = ns("list_block_factors_de"))
                ) 
              ),
              fluidRow(
                column(
                  width = 6,
                  htmlOutput(outputId = ns("select_reference_levels"))
                ),
                column(
                  width = 6,
                  htmlOutput(outputId = ns("select_reference_levels_2"))
                ) 
              ),
              fluidRow(
                column(
                  width = 6,
                  htmlOutput(outputId = ns("select_reference_levels_3"))
                ),
                column(
                  width = 6,
                  htmlOutput(outputId = ns("select_reference_levels_4"))
                ) 
              ),
              fluidRow(
                column(
                  width = 6,
                  htmlOutput(outputId = ns("select_reference_levels_5"))
                ),
                column(
                  width = 6,
                  htmlOutput(outputId = ns("select_reference_levels_6"))
                )  
              ),
              htmlOutput(outputId = ns("list_interaction_terms")),
              textOutput(outputId = ns("experiment_design")),
              tags$head(tags$style(
                "#deg-experiment_design{color: red;font-size: 16px;}"
              )),
              htmlOutput(outputId = ns("list_model_comparisons")),
              actionButton(
                inputId = ns("submit_model_button"),
                label = "Submit & re-calculate",
                style = "float:center"
              ),
              tags$head(tags$style(
                "#deg-submit_model_button{font-size: 20px;}"
              )),
              a(
                h5("More info on DESeq2 experiment design", align = "right"),
                href = "http://rpubs.com/ge600/deseq2",
                target = "_blank"
              )
            ),
            tabPanel(
              title = "Venn Diagram",
              checkboxInput(
                inputId = ns("up_down_regulated"),
                label = "Split gene lists by up- or down-regulation",
                value = FALSE
              ),
              htmlOutput(outputId = ns("listComparisonsVenn")),
              plotOutput(outputId = ns("vennPlot"))
            )
          )
        )
      )
    ),
    tabPanel(
      title = "Analysis (Step 2)",
      sidebarLayout(
        sidebarPanel(
          h5("Examine the results of DEGs for each comparison"),
          htmlOutput(outputId = ns("list_comparisons")),
          br(),
          HTML(
            "<hr style='height:1px;border:none;color:
            #333;background-color:#333;' />"
          ),
          h5("Enrichment analysis for DEGs:"),
          htmlOutput(outputId = ns("select_go")),
          tags$style(
            type = "text/css",
            "#deg-select_go { width:100%;   margin-top:-9px}"
          )
        ),
        mainPanel(
          tabsetPanel(
            tabPanel(
              "Heatmap",
              plotOutput(
                outputId = ns("selected_heatmap"),
              ),
              br(),
              h4("Enriched pathways in DEGs for the selected comparison:"),
              tableOutput(
                outputId = ns("gene_list_go")
              ),
              h4("Top Genes for selected comparison:"),
              tableOutput(
                outputId = ns("gene_list")
              )
            ),
            tabPanel(
              "Volcano Plot",
              plotOutput(
                outputId = ns("volcano_plot")
              )  
            ),
            tabPanel(
              "MA Plot",
              plotOutput(
                outputId = ns("ma_plot")
              )
            ),
            tabPanel(
              "Scatter Plot",
              plotOutput(
                outputId = ns("scatter_plot")
              )
            ),
            tabPanel(
              "Enriched TF Motifs",
              NULL
            ),
            navbarMenu(
              "Enrichment",
              tabPanel(
                "Tree",
                NULL
              ),
              tabPanel(
                "Network",
                NULL
              )
            )
          )
        )
      )
    )
  )
}

#' 06_deg1 Server Functions
#'
#' @noRd
mod_06_deg1_server <- function(id, pre_process) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$data_file_format <- reactive({
      pre_process$data_file_format()
    })
    outputOptions(
      x = output,
      name = "data_file_format",
      suspendWhenHidden = FALSE
    )

    # Experiment Design UI Elements ------------
    output$list_factors_de <- renderUI({
      list_factors_ui(
        sample_info = pre_process$sample_info(),
        data_file_format = pre_process$data_file_format(),
        counts_deg_method = input$counts_deg_method,
        id = id
      )
	  })

    output$list_block_factors_de <- renderUI({ 
      list_block_factors_ui(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model,
        data_file_format = pre_process$data_file_format(),
        deg_method = input$counts_deg_method,
        id = id
      )
	  })

    output$list_model_comparisons <- renderUI({
      req(pre_process$data())

		  list_model_comparisons_ui(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model,
        processed_data = pre_process$data(),
        id = id
      )
	  })

    output$list_interaction_terms <- renderUI({
      list_interaction_terms_ui(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model,
        id = input$id
      )
	  }) 

	
	# set limits for selections of factors. 
#observe({
#		if(length(input$selectFactorsModel) > maxFactors) # less than 4 factors
#			updateCheckboxGroupInput(session, "selectFactorsModel", selected= tail(input$selectFactorsModel,maxFactors))
#
#		if( input$CountsDEGMethod !=3 ) { # if using the limma package
#			if(length(input$selectFactorsModel) > 2) # less than 2 factors
#
#			if(length(input$selectBlockFactorsModel) > 1) # less than 1 factors
#				updateCheckboxGroupInput(session, "selectBlockFactorsModel", selected= tail(input$selectBlockFactorsModel,2))
				
#		}
#		if(length(input$selectComparisonsVenn) >5 )
#					updateCheckboxGroupInput(session, "selectComparisonsVenn", selected= tail(input$selectComparisonsVenn,5))
		
#	})

  })
}

## To be copied in the UI
# mod_06_deg1_ui("06_deg1_ui_1")

## To be copied in the server
# mod_06_deg1_server("06_deg1_ui_1")
