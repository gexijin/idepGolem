#' mod_07_genome UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_07_genome_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "Genome",
    sidebarLayout(
      sidebarPanel(
        htmlOutput(outputId = ns("list_comparisons_genome")),
        tags$style(
          type = "text/css",
          "#genome-list_comparisons_genome{ width:100%;   margin-top: 10px}"
        ),
        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = ns("limma_p_val_viz"),
              label = "Genes: FDR ",
              value = 0.1,
              min = 1e-5,
              max = 1,
              step = .05
            ),
            tippy::tippy_this(
              ns("limma_p_val_viz"),
              "Set the FDR cutoff for showing genes on the chromosome plot.",
              theme = "light-border"
            )
          ),
          column(
            width = 6,
            numericInput(
              inputId = ns("limma_fc_viz"),
              label = "Min. Fold Change",
              value = 2,
              min = 1,
              max = 100,
              step = 0.5
            ),
            tippy::tippy_this(
              ns("limma_fc_viz"),
              "Only plot genes changing at least this many fold.",
              theme = "light-border"
            )
          )
        ),
        fluidRow(
          column(
            width = 6,
            checkboxInput(
              inputId = ns("label_gene_symbol"),
              label = "Label Genes",
              value = FALSE
            ),
            tippy::tippy_this(
              ns("label_gene_symbol"),
              "Show gene symbols next to highlighted loci.",
              theme = "light-border"
            )
          ),
          column(
            width = 6,
            checkboxInput(
              inputId = ns("ignore_non_coding"),
              label = "Coding genes only",
              value = TRUE
            ),
            tippy::tippy_this(
              ns("ignore_non_coding"),
              "Keep only protein-coding genes in the plot.",
              theme = "light-border"
            )
          )
        ),
        fluidRow(
          column(
            width = 6,
            checkboxInput(
              inputId = ns("hide_patches"),
              label = "Hide Patch Chr. ",
              value = TRUE
            ),
            tippy::tippy_this(
              ns("hide_patches"),
              "Hide alternate patch chromosomes from the view.",
              theme = "light-border"
            )
          ),
          column(
            width = 6,
            checkboxInput(
              inputId = ns("hide_chr"),
              label = "Hide Sparse Chrs.",
              value = FALSE
            ),
            tippy::tippy_this(
              ns("hide_chr"),
              "Hide chromosomes that contain four or fewer genes in your dataset.",
              theme = "light-border"
            )
          )
        ),
        HTML(
          "<hr style='height:1px;border:none;
          color:#333;background-color:#333;' />"
        ),
        h5("To identify enriched genomic loci:"),
        fluidRow(
          column(
            width = 6,
            selectInput(
              inputId = ns("ma_window_size"),
              label = "Window Size (Mb)",
              selected = 6,
              choices = c(1, 2, 4, 6, 8, 10, 15, 20),
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("ma_window_size"),
              "Set the width of the sliding window used to detect enriched regions.",
              theme = "light-border"
            )
          ),
          column(
            width = 6,
            selectInput(
              inputId = ns("ma_window_steps"),
              label = "Steps",
              selected = 2,
              choices = c(1, 2, 3, 4),
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("ma_window_steps"),
              "Choose the spacing between adjacent sliding windows.",
              theme = "light-border"
            )
          )
        ),
        selectInput(
          inputId = ns("ch_region_p_val"),
          label = "FDR cutoff for window",
          selected = 0.0001,
          choices = c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001),
          selectize = FALSE
        ),
        tippy::tippy_this(
          ns("ch_region_p_val"),
          "Threshold for calling genomic windows significantly enriched.",
          theme = "light-border"
        ),
        actionButton(
          inputId = ns("chr_data_popup"),
          label = "Download Plot Data",
          icon = icon("download")
        ),
        tippy::tippy_this(
          ns("chr_data_popup"),
          "Download the gene-level data behind the chromosome plot.",
          theme = "light-border"
        )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Chromosomes",
            plotly::plotlyOutput(
              outputId = ns("genome_plotly"),
              height = "900px"
            ),
            p("Select a region to zoom in. Mouse over the points to
            see more information on the gene. Enriched regions are
            highlighted by blue or red line segments paralell to the chromosomes."),
            p("Mouse over the figure to see more options on the top right, including download.")
          ),
          tabPanel(
            icon("info-circle"),
            includeHTML(app_sys("app/www/genome.html"))
          )
        )
      )
    )
  )
}

#' mod_07_genome Server Functions
#'
#' @noRd
mod_07_genome_server <- function(id, pre_process, deg, idep_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$list_comparisons_genome <- renderUI({
      if (is.null(deg$limma()$comparisons)) {
        tagList(
          selectInput(
            inputId = ns("select_contrast"),
            label = NULL,
            choices = list("All" = "All"),
            selected = "All",
            selectize = FALSE
          ),
          tippy::tippy_this(
            ns("select_contrast"),
            "Choose which comparison to display on the genome plots.",
            theme = "light-border"
          )
        )
      } else {
        tagList(
          selectInput(
            inputId = ns("select_contrast"),
            label = "Select a comparison:",
            choices = deg$limma()$comparisons,
            selectize = FALSE
          ),
          tippy::tippy_this(
            ns("select_contrast"),
            "Choose which comparison to display on the genome plots.",
            theme = "light-border"
          )
        )
      }
    })
    
    genome_plot <- reactive({
      req(!is.null(deg$limma()))
      req(!is.null(pre_process$all_gene_info()))
      req(!is.null(chr_data()))
      req(
        input$select_contrast,
        input$limma_p_val_viz,
        input$limma_fc_viz,
        input$ma_window_size,
        input$ma_window_steps,
        input$ch_region_p_val
      )
      withProgress(message = "Generating plot", {
        incProgress(0.2)
        chromosome_plotly(
          limma = deg$limma(),
          select_contrast = input$select_contrast,
          all_gene_info = pre_process$all_gene_info(),
          label_gene_symbol = input$label_gene_symbol,
          ma_window_size = input$ma_window_size,
          ma_window_steps = input$ma_window_steps,
          ch_region_p_val = input$ch_region_p_val,
          hide_patches = input$hide_patches,
          hide_chr = input$hide_chr,
          x = chr_data()$chr_data,
          x0 = chr_data()$other,
          moving_average = chr_data()$enriched_regions,
          ch_length_table = chr_data()$ch_length
        )
      })
    })

    # visualizing fold change on chrs.
    output$genome_plotly <- plotly::renderPlotly({
      genome_plot()
    })
    
    # Popup for chromosome data download options
    observeEvent(input$chr_data_popup, {
      req(!is.null(chr_data()))
      showModal(
        modalDialog(
          title = "Chromosome Data Options",
          # Dataset selection
          selectInput(
            inputId = ns("data_type"),
            label = "Select Dataset",
            choices = c("Enriched Genes" = "enriched_genes",
                        "Chromosome Data" = "chr_data",
                        "Enriched Region Boundaries" = "enriched_regions"),
                        selectize = FALSE
          ),
          # Chromosome data filtering options
          conditionalPanel(
            condition = "input.data_type == 'chr_data'",
            # Gene regulation selection
            selectInput(
              inputId = ns("gene_regulation"),
              label = "Select Gene Regulation",
              choices = c("All", "Up", "Down"),
              selected = "All",
              selectize = FALSE
            ),
            # Chromosome selection
            selectizeInput(
              inputId = ns("chr_select"),
              label = "Select Chromosomes",
              multiple = TRUE,
              choices = "All",
              selected = "All"
            ),
            ns = ns
          ),
          downloadButton(
            ns("download_chr_data"), 
            "Download Data"
            ),
          easyClose = TRUE,
          size = "s",
          footer = modalButton("Close")
        )
      )
    })
    
    # Get chromosome data using user-entered parameters
    chr_data <- reactive({
      req(!is.null(deg$limma()))
      req(!is.null(input$select_contrast))
      req(!is.null(pre_process$all_gene_info()))

      chromosome_data(
        limma = deg$limma(),
        select_contrast = input$select_contrast,
        all_gene_info = pre_process$all_gene_info(),
        ignore_non_coding = input$ignore_non_coding,
        limma_p_val_viz = input$limma_p_val_viz,
        limma_fc_viz = input$limma_fc_viz,
        ma_window_size = input$ma_window_size,
        ma_window_steps = input$ma_window_steps,
        ch_region_p_val = input$ch_region_p_val,
        hide_patches = input$hide_patches,
        hide_chr = input$hide_chr
      )
    })
    
    # Download config for data file download
    output$download_chr_data <- downloadHandler(
      filename = function(){
        req(!is.null(input$data_type))
        paste0(input$data_type, ".csv") #dynamic file name
      },
      content = function(file){
        req(!is.null(chr_data()))
        req(!is.null(input$gene_regulation))
        req(!is.null(input$chr_select))
        
        df <- chr_data()
        # Filter chromosome data
        if (input$data_type == "chr_data"){
          df$chr_data <- chr_filter(chr_data = df$chr_data,
                                    regulation = input$gene_regulation,
                                    chr_select = input$chr_select)
        }
        write.csv(df[[input$data_type]], file)
      }
    )
    
    observeEvent(input$chr_data_popup, {
      # Update chromosome selection dynamically
      updateSelectizeInput(
        inputId = "chr_select",
        session = session,
        choices = c(unique(chr_data()$chr_data$chromosome_name), "All"),
        selected = "All"
      )
    })
  })
}

## To be copied in the UI
# mod_07_genome_ui("mod_07_genome_ui_1")

## To be copied in the server
# mod_07_genome_server("mod_07_genome_ui_1")
