#' 06_deg1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_05_deg_ui <- function(id) {
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
                inputId = ns("limma_p_val"),
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
                inputId = ns("limma_fc"),
                label = h5("Min fold change"),
                value = 2,
                min = 1,
                max = 100,
                step = 0.5
              )
            ),
            tags$style(
              type = "text/css",
              "#deg-limma_p_val { width:100%;   margin-top:-12px}"
            ),
            tags$style(
              type = "text/css",
              "#deg-limma_fc { width:100%;   margin-top:-12px}"
            )
          ),
          actionButton(
            inputId = ns("submit_model_button"),
            label = "Submit & re-calculate",
            style = "float:center"
          ),
          tags$head(tags$style(
            "#deg-submit_model_button{font-size: 20px;}"
          )),
          a(
            h5("Questions?", align = "right"),
            href = "https://idepsite.wordpress.com/degs/",
            target = "_blank"
          )
        ),
        mainPanel(
          tabsetPanel(
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
                htmlOutput(outputId = ns("select_reference_levels"))
              ),
              htmlOutput(outputId = ns("list_interaction_terms")),
              textOutput(outputId = ns("experiment_design")),
              tags$head(tags$style(
                "#deg-experiment_design{color: red;font-size: 16px;}"
              )),
              htmlOutput(outputId = ns("list_model_comparisons")),
              h3("Use the submit button in the sidebar once the desired design is selected!"),
              a(
                h5("More info on DESeq2 experiment design", align = "right"),
                href = "http://rpubs.com/ge600/deseq2",
                target = "_blank"
              )
            ),
            tabPanel(
              title = "Results",
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
              title = "Venn Diagram",
              checkboxInput(
                inputId = ns("up_down_regulated"),
                label = "Split gene lists by up- or down-regulation",
                value = FALSE
              ),
              htmlOutput(outputId = ns("list_comparisons_venn")),
              plotOutput(outputId = ns("venn_plot"))
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
            "#deg-select_go { width:100%; margin-top:-9px}"
          ),
          checkboxInput(
            ns("filtered_background"), 
            label = "Filtered genes as background for enrichment", 
            value = TRUE
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
mod_05_deg_server <- function(id, pre_process) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # DEG STEP 1 ----------
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

	
	  # Set limits for selections of factors
    observe({
    	if(length(input$select_factors_model) > 6) {
        updateCheckboxGroupInput(
          session,
          ns("select_factors_model"),
          selected = tail(input$select_factors_model, 6)
        )
      }
      if(input$counts_deg_method !=3 ) {
        if(length(input$select_factors_model) > 2) {
          updateCheckboxGroupInput(
            session,
            ns("select_factors_model"),
            selected = tail(input$select_factors_model, 2)
          )
        }
        if(length(input$select_block_factors_model) > 1) {
          updateCheckboxGroupInput(
            session,
            ns("select_block_factors_model"),
            selected = tail(input$select_block_factors_model, 1)
          )
        }
      }
      
      if(length(input$select_comparisons_venn) > 5) {
        updateCheckboxGroupInput(
          session,
          ns("select_comparisons_venn"),
          selected = tail(input$select_comparisons_venn, 5)
        )
      }
		})

    output$experiment_design <- renderText({
      experiment_design_txt(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model,
        select_block_factors_model = input$select_block_factors_model,
        select_interactions = input$select_interactions
      )
    })

    output$select_reference_levels <- renderUI({
      select_reference_levels_ui(
        sample_info = pre_process$sample_info(),
        select_factors_model = input$select_factors_model,
        data_file_format = pre_process$data_file_format(),
        counts_deg_method = input$counts_deg_method,
        id = id
      )
    })

    factor_reference_levels <- reactive(
      return(
        c(
          input$reference_level_factor_1,
				  input$reference_level_factor_2,
				  input$reference_level_factor_3,
				  input$reference_level_factor_4,
				  input$reference_level_factor_5,
				  input$reference_level_factor_6
        )
      )
    )

    deg <- reactiveValues(limma = NULL)
    observeEvent(
      input$submit_model_button, {
        req(!is.null(pre_process$raw_counts()))
      
        shinybusy::show_modal_spinner(
          spin = "orbit",
          text = "Running Analysis",
          color = "#000000"
        )

        deg$limma <- limma_value(
          data_file_format = pre_process$data_file_format(),
          counts_deg_method = input$counts_deg_method,
          raw_counts = pre_process$raw_counts(),
          limma_p_val = input$limma_p_val,
          limma_fc = input$limma_fc,
          select_model_comprions = input$select_model_comprions,
          sample_info = pre_process$sample_info(),
          select_factors_model = input$select_factors_model,
          select_interactions = input$select_interactions,
          select_block_factors_model = input$select_block_factors_model,
          factor_reference_levels = factor_reference_levels(),
          processed_data = pre_process$data(),
          counts_log_start = pre_process$counts_log_start(),
          p_vals = pre_process$p_vals()
        )

        shinybusy::remove_modal_spinner()
      }  
    )

    output$sig_gene_stats <- renderPlot({
      req(!is.null(deg$limma$results))
      sig_genes_plot(
        results = deg$limma$results
      )
    })

    output$sig_gene_stats_table <- renderTable({
      req(!is.null(deg$limma))
      
		  genes_stat_table(limma = deg$limma)
    },
      digits = 0,
      spacing = "s",
      include.rownames = FALSE,
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = T
    )

    output$list_comparisons_venn <- renderUI({
      req(!is.null(deg$limma))

      list_comp_venn(
        limma = deg$limma,
        up_down_regulated = input$up_down_regulated,
        id = id
      )
	  })

    output$venn_plot <- renderPlot({
      req(!is.null(deg$limma))
      req(!is.null(input$select_comparisons_venn))
      
		  plot_venn(
        limma = deg$limma,
        up_down_regulated = input$up_down_regulated,
        select_comparisons_venn = input$select_comparisons_venn
      )
    },
      height = 600,
      width = 600
    )

    # DEG STEP 2 --------
    output$list_comparisons <- renderUI({
      
      if (is.null(deg$limma$comparisons)) {
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
          choices = deg$limma$comparisons
	     )
      } 
	  })

    output$selectedHeatmap <- renderPlot({
	   
     bar = selectedHeatmap.data()$bar +2;
		 bar[bar==3] =2

	  	 myheatmap2( selectedHeatmap.data()$genes,bar,200,mycolor=input$heatColors1,c("Down","Up") )
	
    }, height = 400, width = 500)

  })
}

## To be copied in the UI
# mod_06_deg_ui("06_deg_ui_1")

## To be copied in the server
# mod_06_deg_server("06_deg_ui_1")
