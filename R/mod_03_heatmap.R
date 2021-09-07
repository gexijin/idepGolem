#' 03_heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_03_heatmap_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Heatmap",
    sidebarLayout(

      # Heatmap Panel Sidebar ----------
      sidebarPanel(
        sliderInput(
          inputId = ns("n_genes"),
          label = h4("Most variable genes to include:"),
          min = 0,
          max = 12000,
          value = 1000,
          step = 100
        ),
        HTML(
          '<hr style="height:1px;border:none;
           color:#333;background-color:#333;" />'
        ),

        # Heatmap customizing features ----------
        strong("Customize hierarchical clustering (Default values work well):"),
        fluidRow(
          column(width = 3, h5("Color")),
          column(
            width = 9,
            selectInput(
              inputId = ns("heatmap_color_select"),
              label = NULL,
              choices = "green-black-red",
              width = "100%"
            )
          )
        ),
        fluidRow(
          column(width = 4, h5("Distance")),
          column(
            width = 8,
            selectInput(
              inputId = ns("dist_functions"),
              label = NULL,
              choices = "Correlation",
              width = "100%"
            )
          )
        ),
        fluidRow(
          column(width = 4, h5("Linkage")),
          column(
            width = 8,
            selectInput(
              inputId = ns("hclust_functions"),
              label = NULL,
              choices = "average",
              width = "100%"
            )
          )
        ),
        fluidRow(
          column(width = 8, h5("Cut-off Z score")),
          column(
            width = 4,
            numericInput(
              inputId = ns("heatmap_cutoff"),
              label = NULL,
              value = 4,
              min = 2,
              step = 1
            )
          )
        ),

        # Checkbox features ------------
        checkboxInput(
          inputId = ns("gene_centering"),
          label = "Center genes (substract mean)",
          value = TRUE
        ),
        checkboxInput(
          inputId = ns("gene_normalize"),
          label = "Normalize genes (divide by SD)",
          value = FALSE
        ),
        checkboxInput(
          inputId = ns("sample_centering"),
          label = "Center samples (substract mean)",
          value = FALSE
        ),
        checkboxInput(
          inputId = ns("sample_normalize"),
          label = "Normalize samples(divide by SD)",
          value = FALSE
        ),
        checkboxInput(
          inputId = ns("no_sample_clustering"),
          label = "Do not re-order or cluster samples",
          value = FALSE
        ),

        # Sample coloring bar -----------
        htmlOutput(ns("list_factors_heatmap")),
        br(),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/heatmap/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          id = ns("heatmap_panels"),

          # Main heatmap ---------
          tabPanel(
            title = "Heatmap",
            br(),
            plotOutput(outputId = ns("heatmap_main"))
          ),

          # Interactive heatmap panel ----------
          tabPanel(
            title = "Interactive Heatmap",
            br(),
            # INSERT PLOTLY OUTPUT
          ),

          # Gene Standard Deviation Distribution ----------
          tabPanel(
            title = "Gene SD Distribution",
            br(),
            # INSERT SD DISTRIBUTION PLOT
          ),

          # Correlation matrix panel ----------
          tabPanel(
            title = "Correlation Matrix",
            br(),
            # Insert Correlation Matrix
          ),

          # Sample Tree Plot ---------
          tabPanel(
            title = "Sample Tree",
            br(),
            # INSERT SAMPLE TREE PLOT
          )
        )
      )
    )
  )
}

#' 03_heatmap Server Functions
#'
#' @noRd
mod_03_heatmap_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Heatmap Colors ----------
    heatmap_colors <- get_heatmap_colors()
    observe({
      updateSelectInput(
        session = session,
        inputId = "heatmap_color_select",
        choices = heatmap_colors$color_choices
      )
    })

    # Distance functions -----------
    dist_funs <- dist_functions()
    dist_choices <- setNames(
      1:length(dist_funs),
      names(dist_funs)
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "dist_functions",
        choices = dist_choices
      )
    })

    # Hclust Functions ----------
    hclust_funs <- hcluster_functions()
    hclust_choices <- setNames(
      1:length(hclust_funs),
      names(hclust_funs)
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "hclust_functions",
        choices = hclust_choices
      )
    })

    # Sample color bar render ----------
    output$list_factors_heatmap <- renderUI({
      req(!is.null(pre_process$sample_info()))

      selectInput(
        inputId = ns("select_factors_heatmap"),
        label = "Sample color bar:",
        choices = c(colnames(read_sample_info()), "Sample_Name")
      )
    })

  })
}

## To be copied in the UI
# mod_03_heatmap_ui("03_heatmap_ui_1")

## To be copied in the server
# mod_03_heatmap_server("03_heatmap_ui_1")
