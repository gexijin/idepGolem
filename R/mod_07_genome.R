#' mod_07_genome UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_07_genome_ui <- function(id){
  ns <- NS(id)
  tabPanel(
    "Genome",
    sidebarLayout(  
      sidebarPanel( 
        htmlOutput(outputId = ns("list_comparisons_genome")),
        tags$style(
          type = "text/css",
          "#genome-list_comparisons_pathway{ width:100%;   margin-top:-12px}"
        ),
        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = ns("limma_p_val_viz"), 
              label = h5("Genes: FDR "), 
              value = 0.1,
              min = 1e-5, 
              max = 1,
              step = .05
            )
          ),
          column(
            width = 6,
            numericInput(
              inputId = ns("limma_fc_viz"), 
              label = h5("Fold change"), 
              value = 2,
              min = 1,
              max = 100,
              step = 0.5
            )
          )
        ),
        tags$style(
          type = "text/css",
          "#genome-limma_p_val_viz{ width:100%;   margin-top:-12px}"
        ),
        tags$style(
          type = "text/css",
          "#genome-limma_fc_viz{ width:100%;   margin-top:-12px}"
        ),
        fluidRow( 
          column(
            width = 6,
            checkboxInput(
              inputId = ns("label_gene_symbol"),
              label = "Label Genes",
              value = FALSE
            )
          ),
          column(
            width = 6,
            checkboxInput(
              inputId = ns("ignore_non_coding"),
              label = "Coding genes only",
              value = TRUE
            )
          )  
        ),
        HTML(
          "<hr style='height:1px;border:none;
          color:#333;background-color:#333;' />"
        ),
        fluidRow(
          column(
            width = 6,
            selectInput(
              inputId = ns("ma_window_size"),
              label = h5("Window Size (Mb)"),
              selected = 6,
              choices = c(1, 2, 4, 6, 8, 10, 15, 20)
            )
          ),
          column(
            width = 6,
            selectInput(
              inputId = ns("ma_window_steps"),
              label = h5("Steps"),
              selected = 2,
              choices = c(1, 2, 3, 4)
            )
          )
        ),
        selectInput(
          inputId = "chRegionPval", 
          label = h5("FDR cutoff for window"),
          selected = 0.0001,
          choices = c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001)
        ),
        HTML(
          "<hr style='height:1px;border:none;
          color:#333;background-color:#333;' />"
        ),
        actionButton(
          inputId = ns("run_preda"),
          label = "Run PREDA (5 mins)"
        ),
        br(),
        div(
          style = "display:inline-block; float:right",
          actionButton(
            inputId = ns("genome_tab_button"),
            label = "Info"
          )
        )
      ),
      mainPanel(                   
                  bsModal("genomeTabButtonID", "View your genes on chromosomes", "genomeTabButton", size = "large"
                          ,h3("Where are your differentially expressed genes (DEGs) located on the genome?" )
                          ,p("Red and blue dots represent significantly up- and down-regulated genes, respectively, according to the criteria on the side panel. These criteria could differ from the one in DEG1 tab. The distance of the dots from the closest chromosome is proportional to the log2 fold-change (FC).")

                         ,h3("Are there regions of the genome where genes are coherently up- or down-regulated?")
                         ,p("To answer this question, we scan the genome with sliding windows. Within each window we take several steps to slide forward. For example if you choose a window size = 6Mbps and steps = 3, the first window is from 0 to 6 Mbps, the 2nd  from 2 to 8Mbps, and the third from 4 to 10 Mbps, and so on.")
                         ,p("For all genes in a window/region, we test whether the mean of FC of these genes is zero using a t-test. All genes analyzed by DESeq2 or limma, significant or otherwise, are included in this analysis. Hence this result is indepdent of our DEG cutoffs. P values from the test of the mean are adjusted to FDR. Essentially, we considered genes located in a genomic region as a gene set or pathway, and we performed simple pathway analysis by asking whether these genes are behaving non-randomly.")

                         ,p("Based on an FDR cutoff for the windows, red and blue segments indicate genomic regions with genes coherently up- or down-regulated, respectively.  Below you can adjust the window size, and steps in a window, and FDR cutoff for windows.  Mouse over to see gene symbols or IDs. Zoom in regions of interest. The chromosomes may be only partly shown as we use the last gene's location to draw the line.")

                         ,p("As an alternative approach, you can use "
                             ,a("PREDA.", href="https://doi.org/10.1093/bioinformatics/btr404",target="_blank"),
                             "Very slow (5 mins), but may be useful in studying cancer or
                             other diseases that might involve chromosomal gain or loss."  ) 
                   )
           

                  ,plotlyOutput("genomePlotly",height = "900px")

  
        ,bsModal("modalExample111", "Differentially expressed genomic loci", "runPREDA", size="large"
                 ,fluidRow( 
                           column(3, numericInput("RegionsPvalCutoff", 
                                                  label = h5("Min. FDR"), 
                                                  value = 0.01,
                                                  min   = 1e-20,
                                                  max   = 1,
                                                  step  = .05) ),
                           column(3, numericInput("StatisticCutoff",
                                                  label = h5("Min. Statistic"),
                                                  min   = .2, 
                                                  max   = 1.5, 
                                                  value = .5,
                                                  step  =.1) ),

                           column(3,actionButton("showRegions", "Significant Loci")),
                           column(3,actionButton("showGenesRegions", "Genes"))
                 ) #fluidRow
                 ,tags$style(type='text/css', "#showRegions { width:100%; margin-top: 40px;}")
                 ,tags$style(type='text/css', "#showGenesRegions { width:100%; margin-top: 40px;}")
                 ,plotOutput("genomePlot")
        ) #bsModal

        ,bsModal("modalExample15", "Diff. expressed regions", "showRegions", size = "large"
                 ,downloadButton('downloadRegions', 'Download')
                 ,dataTableOutput("chrRegions"))

        ,bsModal("modalExample16", "Diff. expressed regions", "showGenesRegions", size = "large"
                 ,downloadButton('downloadGenesInRegions', 'Download')
                 ,dataTableOutput("genesInChrRegions"))

      )
    )   
  )
}
    
#' mod_07_genome Server Functions
#'
#' @noRd 
mod_07_genome_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_07_genome_ui("mod_07_genome_ui_1")
    
## To be copied in the server
# mod_07_genome_server("mod_07_genome_ui_1")
