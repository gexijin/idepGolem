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
        conditionalPanel(
          condition = "input.PCA_panels == '3D'",
          fluidRow(
            column(
              width = 6,
              selectInput(
                inputId = ns("PCAx3d"),
                label = "X-axis",
                choices = 1:5,
                selected = 1
              )
            ),
            column(
              width = 6,
              selectInput(
                inputId = ns("PCAy3d"),
                label = "Y-axis",
                choices = 1:5,
                selected = 2
              )
            ),
            column(
              width = 6,
              selectInput(
                inputId = ns("PCAz3d"),
                label = "Z-axis",
                choices = 1:5,
                selected = 3
              )
            )
          ),
          p("Use camera icon in figure to download image"),
          ns = ns
        ),
        # select design elements dynamically
        conditionalPanel(
          condition = "input.PCA_panels != 'PCAtools'",
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
            actionButton(
              inputId = ns("seedTSNE"),
              label = "Re-calculate"
            ),
            tags$br()
          ),
          br(),
          br(),
          ns = ns
        ),
        # PCATools plot options
        conditionalPanel(
          condition = "input.PCA_panels == 'PCAtools'",
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
        fluidRow(
          column(
            6,
            # Download report button
            downloadButton(
              outputId = ns("report"),
              label = "Report"
            )
          ),
          column(6,
            offset = 0,
            downloadButton(
              outputId = ns("pca_data"),
              label = "PCA data"
            )
          ),
          tippy::tippy_this(
            ns("report"),
            "Generate HTML report of PCA tab",
            theme = "light-border"
          ),
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
            #h5("This plot is interactive. Hover over the plot to see more details."),
            plotly::plotlyOutput(
              outputId = ns("interactive_pca_plot_obj"),
              width = "100%",
              height = "600px"
            ),
            ottoPlots::mod_download_figure_ui(ns("download_interactive_pca")),
            br(),
            br(),
            shiny::textOutput(
              outputId = ns("pc_correlation")
            ),
            br(),
          ),
          tabPanel(
            title = "3D",
            plotly::plotlyOutput(
              outputId = ns("pca_plot_obj_3d"),
              width = "100%",
              height = "700px"
            ),
            br(),
            br(),
            shiny::textOutput(
              outputId = ns("pc_correlation_3d")
            ),
          ),
          tabPanel(
            "PCAtools",
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
              outputId = ns("var_imp_x"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(ns("download_imp_x")),
            br(),
            br(),
            br(),
            plotOutput(
              outputId = ns("var_imp_y"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(ns("download_imp_y")),
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
            plotly::plotlyOutput(
              outputId = ns("mds_plot_obj"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(ns("download_mds")),
          ),
          tabPanel(
            "t-SNE",
            br(),
            plotly::plotlyOutput(
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
mod_04_pca_server <- function(id, load_data, pre_process, idep_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # Store client info in a convenience variable

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
        selected_color = input$selectFactors1,
        plots_color_select = load_data$plots_color_select()
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })
    output$interactive_pca_plot_obj <- plotly::renderPlotly({
      # remove legend titles for the interactive plot
      p <- pca_plot() +
        ggplot2::labs(color = NULL, shape = NULL)
      plotly::ggplotly(p, tooltip = "text")
    })
    # Download Button
    download_interactive_pca <- ottoPlots::mod_download_figure_server(
      id = "download_interactive_pca",
      filename = "pca_plot",
      figure = reactive({
        pca_plot()
      }),
      label = "",
      width = get_plot_width(
        client_data = session$clientData,
        plot_name = "interactive_pca_plot_obj",
        tab = id
      ),
      height = get_plot_height(
        client_data = session$clientData,
        plot_name = "interactive_pca_plot_obj",
        tab = id
      )
    )
    # PC Factor Correlation ---------
    output$pc_correlation <- renderText({
      req(!is.null(pre_process$data()))
      pc_factor_correlation(
        data = pre_process$data(),
        sample_info = pre_process$sample_info()
      )
    })

    # PCA plot 3D ------------
    # reactive part -----
    pca_plot_3d <- reactive({
      req(!is.null(pre_process$data()))

      p <- PCA_plot_3d(
        data = pre_process$data(),
        sample_info = pre_process$sample_info(),
        PCAx = input$PCAx3d,
        PCAy = input$PCAy3d,
        PCAz = input$PCAz3d,
        selected_shape = input$selectFactors2,
        selected_color = input$selectFactors1,
        plots_color_select = load_data$plots_color_select()
      )
    })
    output$pca_plot_obj_3d <- plotly::renderPlotly({
      pca_plot_3d()
    })

    # PC Factor Correlation ---------
    output$pc_correlation_3d <- renderText({
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

      p <- t_SNE_plot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info(),
        selected_shape = input$selectFactors2,
        selected_color = input$selectFactors1,
        plots_color_select = load_data$plots_color_select()
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })
    output$t_sne <- plotly::renderPlotly({
      req(t_SNE_plot_obj())
      plotly::ggplotly(
        t_SNE_plot_obj(),
        tooltip = "text"
      )
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

      p <- MDS_plot(
        data = pre_process$data(),
        sample_info = pre_process$sample_info(),
        selected_shape = input$selectFactors2,
        selected_color = input$selectFactors1,
        plots_color_select = load_data$plots_color_select()
      )
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })
    
    output$mds_plot_obj <- plotly::renderPlotly({
      req(mds_plot())
      plotly::ggplotly(mds_plot(), tooltip = "text")
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

        p <- PCA_biplot(
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
        #        refine_ggplot2(p, gridline = pre_process$plot_grid_lines())
      })
    })

    output$pcatools_biplot <- renderPlot({
      req(biplot())
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

      p <- PCA_Scree(
        processed_data = pre_process$data()
      )
      #      refine_ggplot2(p, gridline = pre_process$plot_grid_lines())
    })
    output$pcatools_scree <- renderPlot({
      req(scree())
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
      refine_ggplot2(
        p = p,
        gridline = pre_process$plot_grid_lines(),
        ggplot2_theme = pre_process$ggplot2_theme()
      )
    })
    output$pcatools_eigencor <- renderPlot({
      req(eigencor())
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
    
    # Variable importance for 1st selected PC
    var_plot1 <- reactive({
      req(!is.null(pre_process$data()))
      
      var_imp_plots(pre_process$data(),
                    pre_process$all_gene_names(),
                    input$x_axis_pc)
    })
    
    output$var_imp_x <- renderPlot({
      req(var_plot1())
      
      return(var_plot1())
    })
    
    # Download Button
    download_imp_x <- ottoPlots::mod_download_figure_server(
      id = "download_imp_x",
      filename = "var_imp_x",
      figure = reactive({
        var_plot1()
      }),
      label = "",
      width = 10,
      height = 6
    )
    
    # Variable importance for 2nd selected PC
    var_plot2 <- reactive({
      req(!is.null(pre_process$data()))
      
      var_imp_plots(pre_process$data(),
                    pre_process$all_gene_names(),
                    input$y_axis_pc)
    })
    
    output$var_imp_y <- renderPlot({
      req(var_plot2())
      
      return(var_plot2())
    })
    
    # Download Button
    download_imp_y <- ottoPlots::mod_download_figure_server(
      id = "download_imp_y",
      filename = "var_imp_y",
      figure = reactive({
        var_plot2()
      }),
      label = "",
      width = 10,
      height = 6
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

    # Download PCA table --------------
    pca_data <- reactive({
      get_pc(
        data = pre_process$data(),
        sample_info = pre_process$sample_info()
      )
    })

    output$pca_data <- downloadHandler(
      filename = function() {
        paste(gsub("-", "_", Sys.Date()), "_pca_data.csv", sep = "")
      },
      content <- function(file) {
        write.csv(pca_data(), file = file)
      }
    )


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

          markdown_location <- app_sys("app/www/RMD/pca_workflow.Rmd")
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
            ui_shape = input$selectShape,
            plots_color_select = load_data$plots_color_select()
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
