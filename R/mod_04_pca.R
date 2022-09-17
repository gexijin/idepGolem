#' 04_pca UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList

mod_04_pca_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    title = "PCA",
    sidebarLayout(
      sidebarPanel(
        # width of shaded part of screen
        # width = 3,
        conditionalPanel(
          condition = "input.PCA_panels == 'PCA'",
          fluidRow(
            column(
              width = 6,
              selectInput(
                inputId = ns("PCAx"),
                label = "X-axis",
                choices = 1:5,
                selected = 1
              )
            ),
            column(
              width = 6,
              selectInput(
                inputId = ns("PCAy"),
                label = "Y-axis",
                choices = 1:5,
                selected = 2
              )
            )
          ),
          ns = ns
        ),
        # select design elements dynamically
        conditionalPanel(
          condition = "input.PCA_panels != 'PCAtools Package'",
          fluidRow(
            column(
              width = 6,
              uiOutput(
                outputId = ns("listFactors2")
              )
            ),
            column(
              width = 6,
              uiOutput(
                outputId = ns("listFactors1")
              )
            )
          ),
          ns = ns
        ),
        conditionalPanel(
          condition = "input.PCA_panels == 't-SNE'",
          fluidRow(
            actionButton(inputId = ns("seedTSNE"), label = "Re-calculate from new seed")
          ),
          ns = ns
        ),
        # PCATools plot options
        conditionalPanel(
          condition = "input.PCA_panels == 'PCAtools Package'",
          fluidRow(
            column(
              width = 6,
              selectInput(
                inputId = ns("x_axis_pc"),
                label = "X-axis",
                choices = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                selected = "PC1"
              )
            ),
            column(
              width = 6,
              selectInput(
                inputId = ns("y_axis_pc"),
                label = "Y-axis",
                choices = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                selected = "PC2"
              )
            )
          ),

          # Dynamic Color and Shape options
          fluidRow(
            column(
              width = 6,
              uiOutput(
                outputId = ns("pcatools_shape")
              )
            ),
            column(
              width = 6,
              uiOutput(
                outputId = ns("pcatools_color")
              )
            )
          ),
          # plot customization
          checkboxInput(
            inputId = ns("showLoadings"),
            label = "Show Loadings",
            value = FALSE
          ),
          checkboxInput(
            inputId = ns("encircle"),
            label = "Encircle",
            value = FALSE
          ),
          checkboxInput(
            inputId = ns("pointLabs"),
            label = "Point Labels",
            value = TRUE
          ),
          numericInput(
            inputId = ns("pointSize"),
            label = "Point Size (1-10)",
            value = 3.0,
            min = 1,
            max = 15
          ),
          ns = ns
        ),
        # Download report button
        downloadButton(
          outputId = ns("report"),
          label = "Generate Report"
        ),
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/pca/",
          target = "_blank"
        )
      ),
      mainPanel(
        tabsetPanel(
          id = ns("PCA_panels"),
          tabPanel(
            title = "PCA",
            plotOutput(
              outputId = ns("pca_plot_obj"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(ns("download_pca")),
            br(),
            br(),
            shiny::textOutput(
              outputId = ns("pc_correlation")
            ),
            br(),
            shiny::verbatimTextOutput(ns("image_dimensions")),
          ),
          tabPanel(
            "PCAtools Package",
            br(),
            plotOutput(
              outputId = ns("pcatools_biplot"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(ns("download_biplot")),
            br(),
            br(),
            br(),
            plotOutput(
              outputId = ns("pcatools_scree"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(ns("download_scree")),
            br(),
            br(),
            br(),
            plotOutput(
              outputId = ns("pcatools_eigencor"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(ns("download_eigencor")),
            br(),
            br()
          ),
          tabPanel(
            "MDS",
            br(),
            plotOutput(
              outputId = ns("mds_plot_obj"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(ns("download_mds")),
          ),
          tabPanel(
            "t-SNE",
            br(),
            plotOutput(
              outputId = ns("t_sne"),
              width = "100%",
              height = "500px"
            ),
            br(),
            ottoPlots::mod_download_figure_ui(ns("download_t_sne")),
            br()
          )
          # tabPanel(
          #   "Pathway Analysis of PCA",
          #   NULL
          # )
        )
      )
    )
  )
}
#' 05_pca Server Functions
#'
#' @noRd
mod_04_pca_server <- function(id, pre_process, idep_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Store client info in a convenience variable
    cdata <- session$clientData
    # get pca image dimensions
    output$image_dimensions <- renderText({
      a_ratio <- cdata[["output_pca-pca_plot_obj_width"]] / cdata[["output_pca-pca_plot_obj_height"]]
      paste(
        "Plot size (pixels): ",
        cdata[["output_pca-pca_plot_obj_width"]],
        " x ",
        cdata[["output_pca-pca_plot_obj_height"]],
        "\nAspect Ratio: ", cdata[["output_pca-pca_plot_obj_width"]] / cdata[["output_pca-pca_plot_obj_height"]],
        "\nPlot size (inches): ", "6.5 x ", round(6.5 / a_ratio, 3)
      )
    })

    # PCA plot ------------
    # reactive part -----
    pca_plot <- reactive({
      req(!is.null(pre_process$data()))

      p <- PCA_plot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info(),
        PCAx = input$PCAx,
        PCAy = input$PCAy,
        selected_shape = input$selectFactors2,
        selected_color = input$selectFactors1
      )
    })
    output$pca_plot_obj <- renderPlot({
      print(pca_plot())
    })

    # Download Button
    download_pca <- ottoPlots::mod_download_figure_server(
      id = "download_pca",
      filename = "pca_plot",
      figure = reactive({
        pca_plot()
      }),
      label = ""
    )

    # PC Factor Correlation ---------
    output$pc_correlation <- renderText({
      req(!is.null(pre_process$data()))
      pc_factor_correlation(
        data = pre_process$data(),
        sample_info = pre_process$sample_info()
      )
    })


    # t_SNE plot -----------------
    t_SNE_plot_obj <- reactive({
      req(!is.null(pre_process$data()))

      input$seedTSNE

      t_SNE_plot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info(),
        selected_shape = input$selectFactors2,
        selected_color = input$selectFactors1
      )
    })
    output$t_sne <- renderPlot({
      print(t_SNE_plot_obj())
    })
    # Download Button
    download_t_sne <- ottoPlots::mod_download_figure_server(
      id = "download_t_sne",
      filename = "t_sne_plot",
      figure = reactive({
        t_SNE_plot_obj()
      }),
      label = ""
    )

    # MDS plot ------------

    mds_plot <- reactive({
      req(!is.null(pre_process$data()))

      MDS_plot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info(),
        selected_shape = input$selectFactors2,
        selected_color = input$selectFactors1
      )
    })
    output$mds_plot_obj <- renderPlot({
      print(mds_plot())
    })
    # Download Button
    download_mds <- ottoPlots::mod_download_figure_server(
      id = "download_mds",
      filename = "mds_plot",
      figure = reactive({
        mds_plot()
      }),
      label = ""
    )

    # PCAtools biplot  ---------------------
    biplot <- reactive({
      req(!is.null(pre_process$data()))
      withProgress(message = "Generating Plots", {
        incProgress(0.2)

        PCA_biplot(
          data = pre_process$data(),
          sample_info = pre_process$sample_info(),
          select_gene_id = pre_process$select_gene_id(),
          all_gene_names = pre_process$all_gene_names(),
          selected_x = input$x_axis_pc,
          selected_y = input$y_axis_pc,
          encircle = input$encircle,
          showLoadings = input$showLoadings,
          pointlabs = input$pointLabs,
          point_size = input$pointSize,
          ui_color = input$selectColor,
          ui_shape = input$selectShape
        )
      })
    })

    output$pcatools_biplot <- renderPlot({
      print(biplot())
    })


    # Download Button
    download_biplot <- ottoPlots::mod_download_figure_server(
      id = "download_biplot",
      filename = "biplot",
      figure = reactive({
        biplot()
      }),
      label = ""
    )

    # PCAtools Scree Plot --------------------
    scree <- reactive({
      req(!is.null(pre_process$data()))

      PCA_Scree(
        processed_data = pre_process$data()
      )
    })
    output$pcatools_scree <- renderPlot({
      print(scree())
    })

    # Download Button
    download_scree <- ottoPlots::mod_download_figure_server(
      id = "download_scree",
      filename = "scree",
      figure = reactive({
        scree()
      }),
      label = ""
    )


    # PCAtools Eigencor Plot --------------------
    eigencor <- reactive({
      req(!is.null(pre_process$data()))

      p <- PCAtools_eigencorplot(
        processed_data = pre_process$data(),
        sample_info = pre_process$sample_info()
      )
      return(p)
    })
    output$pcatools_eigencor <- renderPlot({
      print(eigencor())
    })

    # Download Button
    download_eigencor <- ottoPlots::mod_download_figure_server(
      id = "download_eigencor",
      filename = "eigencor",
      figure = reactive({
        eigencor()
      }),
      label = ""
    )
    # select color
    output$listFactors1 <- renderUI({
      req(!is.null(pre_process$data()))

      if (is.null(pre_process$sample_info())) {
        return(NULL)
      } else {
        selectInput(
          inputId = ns("selectFactors1"),
          label = "Color ",
          choices = c(colnames(pre_process$sample_info()), "Names"),
          selected = "Names"
        )
      }
    })

    # select shape
    output$listFactors2 <- renderUI({
      req(!is.null(pre_process$data()))

      if (is.null(pre_process$sample_info())) {
        return(NULL)
      } else {
        tem <- c(colnames(pre_process$sample_info()), "Names")
        selectInput(
          inputId = ns("selectFactors2"),
          label = "Shape",
          choices = tem,
          selected = "Names"
        )
      }
    })

    # select color & shape for pcatools
    output$pcatools_color <- renderUI({
      req(!is.null(pre_process$data()))
      if (is.null(pre_process$sample_info())) {
        return(NULL)
      }

      selectInput(
        inputId = ns("selectColor"),
        label = "Color",
        choices = colnames(pre_process$sample_info()),
        selected = colnames(pre_process$sample_info())[1]
      )
    })
    output$pcatools_shape <- renderUI({
      req(!is.null(pre_process$data()))

      if (is.null(pre_process$sample_info())) {
        return(NULL)
      }
      selectInput(
        inputId = ns("selectShape"),
        label = "Shape",
        choices = colnames(pre_process$sample_info()),
        selected = colnames(pre_process$sample_info())[1]
      )
    })



    # Markdown report------------
    output$report <- downloadHandler(

      # For PDF output, change this to "report.pdf"
      filename = "pca_report.html",
      content = function(file) {
        withProgress(message = "Generating Report", {
          incProgress(0.2)


          # Copy the report file to a temporary directory before processing it, in
          # case we don't have write permissions to the current working dir (which
          # can happen when deployed).
          tempReport <- file.path(tempdir(), "pca_workflow.Rmd")
          # tempReport
          tempReport <- gsub("\\", "/", tempReport, fixed = TRUE)

          # This should retrieve the project location on your device:
          # "C:/Users/bdere/Documents/GitHub/idepGolem"
          wd <- getwd()

          markdown_location <- paste0(wd, "/vignettes/Reports/pca_workflow.Rmd")
          file.copy(from = markdown_location, to = tempReport, overwrite = TRUE)
          # Set up parameters to pass to Rmd document
          params <- list(
            pre_processed_data = pre_process$data(),
            pre_processed_descr = pre_process$descr(),
            sample_info = pre_process$sample_info(),
            pc_x = input$PCAx,
            pc_y = input$PCAy,
            color = input$selectFactors1,
            shape = input$selectFactors2,
            all_gene_names = pre_process$all_gene_names(),
            select_gene_id = pre_process$select_gene_id(),
            selected_x = input$x_axis_pc,
            selected_y = input$y_axis_pc,
            encircle = input$encircle,
            showLoadings = input$showLoadings,
            pointlabs = input$pointLabs,
            point_size = input$pointSize,
            ui_color = input$selectColor,
            ui_shape = input$selectShape
          )

          # stops report generation if params are missing
          req(params)
          # Knit the document, passing in the `params` list, and eval it in a
          # child of the global environment (this isolates the code in the document
          # from the code in this app).
          rmarkdown::render(
            input = tempReport, # markdown_location,
            output_file = file,
            params = params,
            envir = new.env(parent = globalenv())
          )
        })
      }
    )
  })
}

## To be copied in the UI
# mod_05_pca_ui("05_pca_ui_1")

## To be copied in the server
# mod_05_pca_server("05_pca_ui_1")
