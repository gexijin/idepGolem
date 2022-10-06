#' 01_load_data UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_01_load_data_ui <- function(id) {
  ns <- shiny::NS(id)
  tabPanel(
    title = "Load Data",
    sidebarLayout(

      ##################################################################
      #       Load Data sidebar panel ----
      ##################################################################

      sidebarPanel(
        fluidRow(
          column(
            width = 9,
            p("Load a demo file below.
            Click the tabs to see some magic!")
          ),
          column(
            width = 3,
            # Reset Button -----------
            p(htmltools::HTML(
              "<div align=\"right\"><A HREF=\"javascript:history.go(0)\"
              >Reset</A></div>"
            )),
          )
        ),


        # Species Match Drop Down ------------
        # strong("1. Optional: Select or search for species"),
        fluidRow(
          column(
            width = 9,
            selectInput(
              inputId = ns("select_org"),
              label = strong("1. Optional: Select or search for species"),
              choices = " ",
              multiple = FALSE,
              selectize = TRUE
            )
          ),
          column(
            width = 3,
            # Species list and genome assemblies ----------
            actionButton(
              inputId = ns("genome_assembl_button"),
              label = "Info"
            )
          ),
          tippy::tippy_this(
            ns("genome_assembl_button"),
            "Additional info on annotated species",
            theme = "light-border"
          )
        ),

        # Conditional .GMT file input bar ----------
        conditionalPanel(
          condition = 'input.select_org == "NEW"',
          fileInput(
            inputId = ns("gmt_file"),
            label =
              "Upload a geneset .GMT file for enrichment analysis (optional)",
            accept = c(
              "text/csv",
              "text/comma-separated-values",
              "text/tab-separated-values",
              "text/plain",
              ".csv",
              ".tsv"
            )
          ),
          ns = ns
        ),

        # Buttons for data file format ----------
        fluidRow(
          column(
            width = 9,
            selectInput(
              inputId = ns("data_file_format"),
              label = strong("2. Choose data type"),
              choices = list(
                "Read counts data (recommended)" = 1,
                "Normalized expression values (RNA-seq FPKM, microarray, etc.)" = 2,
                "Fold-changes and corrected P values from CuffDiff or any other
                program" = 3
              ),
              selected = 1
            )
          ),
          column(
            width = 3,
            actionButton(
              inputId = ns("data_format_help"),
              label = "Info"
            ),
            tippy::tippy_this(
              ns("data_format_help"),
              "Additional info on accepted data types",
              theme = "light-border"
            )
          )
        ),
        # Conditional panel for fold changes data file ----------
        conditionalPanel(
          condition = "input.data_file_format == 3",
          checkboxInput(
            inputId = ns("no_fdr"),
            label = "Fold-changes only",
            value = FALSE
          ),
          ns = ns
        ),

        # Load expression data options ----------
        # Includes load demo action button, demo data dropdown, and expression
        # file upload box
        uiOutput(ns("load_data_ui")),

        # alternative UI output message for once expression data is loaded
        uiOutput(ns("load_data_alt")),

        # Experiment design file input ----------
        uiOutput(ns("design_file_ui")),
        fluidRow(
          column(
            width = 8,
            # Yes or no to converting IDs -------------
            checkboxInput(
              inputId = ns("no_id_conversion"),
              label = "Do not convert gene IDs",
              value = FALSE
            )
          ),
          column(
            width = 4,
            actionButton(ns("customize_button"), "Plot Options")
          )
        ),
        tippy::tippy_this(
          ns("customize_button"),
          "Customize plots throughout app",
          theme = "light-border"
        ),
        shinyBS::bsModal(
          id = ns("modalExample"),
          title = "Plot options",
          trigger = ns("customize_button"),
          size = "small",
          selectInput(
            inputId = ns("heatmap_color_select"),
            label = "Heatmap Color scheme:",
            choices = c(
              "Green-Black-Red",
              "Red-Black-Green",
              "Blue-White-Red",
              "Green-Black-Magenta",
              "Blue-Yellow-Red",
              "Blue-White-Brown",
              "Orange-White-Blue"
            ),
            selected = "Green-Black-Red",
            width = "100%"
          ),
          selectInput(
            inputId = ns("select_gene_id"),
            label = "Gene ID type for plots",
            choices = c("symbol", "ensembl_ID", "User_ID"),
            selected = "symbol"
          )
        ),

        # Link to public RNA-seq datasets ----------
        a(
          h4("Public RNA-seq datasets"),
          href = "http://bioinformatics.sdstate.edu/reads/"
        ),
        # Table output for species loading progress -----------
        tableOutput(ns("species_match")),
        a(
          h5("Citation information", align = "right"),
          href = "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6#citeas",
          target = "_blank"
        ),

        # Action button for Gene ID examples -----------
        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/data-format/",
          target = "_blank"
        )
      ),


      ##################################################################
      #       Load Data panel main ----
      ##################################################################
      mainPanel(
        shinyjs::useShinyjs(),

        # Table output for sample tissue type ----------
        DT::dataTableOutput(ns("sample_info_table")),
        br(),
        br(),

        # Display first 20 rows of the data ----------
        DT::dataTableOutput(ns("sample_20")),
        div(
          id = ns("load_message"),
          h4("Loading R packages, please wait ... ... ...")
        ),

        # Hide welcome screen after data is loaded -----
        uiOutput(ns("welcome_ui")),

        # Display file format help html document when prompted ----
        uiOutput(ns("format_help_ui"))
      )
    )
  )
}








#' 01_load_data Server Functions
#'
#' @noRd
### testing something, come back to later
mod_01_load_data_server <- function(id, idep_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # increase max input file size
    options(shiny.maxRequestSize = 2001024^2)

    welcome_modal <- shiny::modalDialog(
      title = "Welcome to iDEP!",
      tags$p("We hope that you find our application useful. If iDEP is used,
      even for preliminrary analysis, please cite:"),
      tags$h4("Ge, Son & Yao, iDEP: an integrated web application for
      differential expression and pathway analysis of RNA-Seq data,
      BMC Bioinformatics 19:1-24, 2018.", style = "color:#6B1518"),
      tags$p("By citing the iDEP paper, you will help this service remain
      available in the future. Visit the 'About' tab for more information."),
      tags$a(
        h5("Download citation"),
        href = "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6#citeas",
        target = "_blank"
      ),
      tags$h4("How-to videos coming soon!"),
      easyClose = T,
      size = "l"
    )

    shiny::showModal(welcome_modal)

    # Pop-up modal for gene assembl information ----
    observeEvent(input$genome_assembl_button, {
      shiny::showModal(
        shiny::modalDialog(
          size = "l",
          p("Search annotated species by common or scientific names,
          or NCBI taxonomy id. If your species cannot be found here,
          you can still use iDEP without pathway analysis."),
          DT::renderDataTable({
            df <- idep_data$org_info[, c("ensembl_dataset", "name", "totalGenes")]
            colnames(df) <- c("Ensembl/STRING-db ID", "Name (Assembly)", "Total Genes")
            row.names(df) <- NULL
            DT::datatable(
              df,
              options = list(
                pageLength = 20,
                scrollY = "400px"
              ),
              rownames = FALSE
            )
          })
        )
      )
    })


    # UI elements for load demo action button, demo data drop down, and -----
    # expression file upload
    output$load_data_ui <- renderUI({
      req(
        (is.null(input$go_button) || input$go_button == 0) &&
          is.null(input$expression_file)
      )
      req(input$data_file_format)

      # get demo data files based on specified format
      files <- idep_data$demo_file_info
      files <- files[files$type == input$data_file_format, ]
      choices <- setNames(as.list(files$ID), files$name)

      tagList(
        fluidRow(
          column(
            width = 6,
            actionButton(
              inputId = ns("go_button"),
              label = "Load demo:"
            ),
            tags$head(tags$style(
              "#load_data-go_button{color: red;
              font-size: 16px;}"
            )),
            tippy::tippy_this(
              ns("go_button"),
              "Load the selected demo data",
              theme = "light-border"
            )
          ),
          column(
            width = 6,
            selectInput(
              inputId = ns("select_demo"),
              label = NULL,
              choices = choices,
              selected = choices[[1]]
            ),
          )
        ),

        # Expression data file input
        fileInput(
          inputId = ns("expression_file"),
          label = strong("3. Upload expression data (CSV or text)"),
          accept = c(
            "text/csv",
            "text/comma-separated-values",
            "text/tab-separated-values",
            "text/plain",
            ".csv",
            ".tsv"
          )
        )
      )
    })

    # Alternate ui message to reset app once data is loaded ----
    output$load_data_alt <- renderUI({
      req(!(
        (is.null(input$go_button) || input$go_button == 0) &&
          is.null(input$expression_file)
      ))

      # reset message and action button
      tagList(
        h5("To load new files, reset the application from above.")
      )
    })

    # UI element for design file upload ----
    output$design_file_ui <- renderUI({
      req(is.null(input$go_button) || input$go_button == 0)

      tagList(
        fileInput(
          inputId = ns("experiment_file"),
          label = strong("4. Optional: Upload an experiment design file(CSV or text)"),
          accept = c(
            "text/csv",
            "text/comma-separated-values",
            "text/tab-separated-values",
            "text/plain",
            ".csv",
            ".tsv"
          )
        )
      )
    })


    # Provide species list for dropdown selection -----------
    observe({
      updateSelectizeInput(
        session = session,
        inputId = "select_org",
        choices = idep_data$species_choice,
        selected = idep_data$species_choice[1],
        server = TRUE
      )
    })

    # Show messages when on the Network tab or button is clicked ----
    observe({
      req(is.null(loaded_data()$data) && (
        tab() != "Load Data" || tab() != "About"
      ))

      showNotification(
        ui = paste("Pleaes load a demo file or your own data first."),
        id = "load_data_first",
        duration = NULL,
        type = "error"
      )
    })

    # Remove messages if the tab changes --------
    observe({
      req(!is.null(loaded_data()$data) ||
        tab() == "Load Data" || tab() == "About")
      removeNotification("load_data_first")
    })

    # Message for the status of the app ---------
    output$file_format <- renderUI({
      shinyjs::hideElement(id = "load_message")
      i <- "<h4>Ready to load data files.</h4>"
      htmltools::HTML(paste(i, collapse = "<br/>"))
    })

    # Change demo data based on selected format ----
    # returns a vector with file names  c(data, design)
    demo_data_file <- reactive({
      req(input$select_demo)

      files <- idep_data$demo_file_info
      ix <- which(files$ID == input$select_demo)
      return(c(
        files$expression[ix],
        files$design[ix]
      ))
    })

    # Reactive element to load the data from the user or demo data ---------
    loaded_data <- reactive({
      req(!is.null(input$go_button))
      input_data(
        expression_file = input$expression_file,
        experiment_file = input$experiment_file,
        go_button = input$go_button,
        demo_data_file = demo_data_file()[1],
        demo_metadata_file = demo_data_file()[2]
      )
    })

    # Sample information table -----------
    output$sample_info_table <- DT::renderDataTable({
      req(!is.null(loaded_data()$sample_info))

      DT::datatable(
        isolate({
          tem <- t(loaded_data()$sample_info)
          tem <- cbind(rownames(tem), tem)
          colnames(tem)[1] <- "Study_design"
          tem
        }),
        options = list(
          pageLength = 10,
          scrollX = "400px",
          dom = "t",
          ordering = F
        ),
        rownames = FALSE
      )
    })

    # First 20 rows of dataset table -----------
    output$sample_20 <- DT::renderDataTable({
      req(!is.null(conversion_info()$converted_data))

      DT::datatable(
        conversion_info()$converted_data[1:20, ],
        options = list(
          pageLength = 10,
          scrollX = "400px",
          dom = "t"
        ),
        rownames = TRUE
      )
    })

    # Get converted IDs ----------
    conversion_info <- reactive({
      req(!is.null(loaded_data()$data))

      # test data for correct format
      if (
        min(loaded_data()$data, na.rm = TRUE) < 0 & input$data_file_format == 1
      ) {
        showModal(modalDialog(
          title = "Somthing seems incorrect...",
          tags$p("Negative values were detected in this dataset. This is not
                 correct for the selected data type (Read Counts).
                 Please double check data type or the data.
                 You will not be able to continue until one of these
                 are resolved."),
          tags$br(),
          tags$p("Upon clicking okay, the application will reset."),
          size = "m",
          footer = actionButton(ns("reset_app"), "Okay, I will check my inputs!")
        ))
      } else {
        shinybusy::show_modal_spinner(
          spin = "orbit",
          text = "Loading Data",
          color = "#000000"
        )

        converted <- convert_id(
          rownames(loaded_data()$data),
          idep_data = idep_data,
          select_org = input$select_org,
          max_sample_ids = 200
        )

        all_gene_info <- gene_info(
          converted = converted,
          select_org = input$select_org,
          idep_data = idep_data
        )

        converted_data <- convert_data(
          converted = converted,
          no_id_conversion = input$no_id_conversion,
          data = loaded_data()$data
        )

        all_gene_names <- get_all_gene_names(
          mapped_ids = converted_data$mapped_ids,
          all_gene_info = all_gene_info
        )

        gmt_choices <- gmt_category(
          converted = converted,
          converted_data = converted_data$data,
          select_org = input$select_org,
          gmt_file = input$gmt_file,
          idep_data = idep_data
        )

        shinybusy::remove_modal_spinner()

        return(list(
          converted = converted,
          all_gene_info = all_gene_info,
          converted_data = converted_data$data,
          all_gene_names = all_gene_names,
          gmt_choices = gmt_choices
        ))
      }
    })

    observeEvent(input$reset_app, {
      session$reload()
    })

    # Species match table ----------
    output$species_match <- renderTable(
      {
        req(!is.null(input$go_button))
        if (is.null(input$expression_file) && input$go_button == 0) {
          return(NULL)
        }
        isolate({
          if (is.null(conversion_info()$converted)) {
            return(as.data.frame("ID not recognized."))
          }
          tem <- conversion_info()$converted$species_match
          if (nrow(tem) > 50) { # show only 50
            tem <- tem[1:50, , drop = FALSE]
          }
          if (is.null(tem)) {
            as.data.frame("ID not recognized.")
          } else {
            data.frame(
              "Species(genes matched)" = tem[, 1],
              check.names = FALSE
            )
          }
        })
      },
      digits = -1,
      spacing = "s",
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = TRUE
    )

    # Species match message ----------
    observe({
      req(
        tab() == "Load Data" &&
          !is.null(conversion_info()$converted) &&
          input$select_org == idep_data$species_choice[[1]] # species not selected
      )

      tem <- conversion_info()$converted$species_match
      showNotification(
        ui = paste0("Matched species is '", tem[1, ], ".' If that is not your
                    species, please click Reset and use the dropdown to select
                    the correct species first."),
        id = "species_match",
        duration = NULL,
        type = "warning"
      )
    })


    # Remove message if the tab changes --------
    observe({
      req(tab() != "Load Data")

      removeNotification("species_match")
    })

    output$welcome_ui <- renderUI({
      req(
        input$go_button == 0 &
          !is.null(input$go_button) &
          input$data_format_help == 0
      )

      tagList(
        fluidRow(
          column(
            width = 5,
            h3("Welcome to iDEP!")
          ),
          column(
            width = 6,
            img(
              src = "www/idep_logo.png",
              width = "43",
              height = "50"
            )
          )
        ),
        htmlOutput(ns("file_format")),
        includeHTML("inst/app/www/messages.html"),
        br(),
        img(
          src = "www/flowchart.png",
          align = "center",
          width = "562",
          height = "383"
        ),
        br(),
        img(
          src = "www/figs.gif",
          align = "center",
          width = "640",
          height = "480"
        )
      )
    })

    output$format_help_ui <- renderUI({
      req(input$data_format_help != 0)

      includeHTML("inst/app/www/format.html")
    })


    # Return data used in the following panels --------
    list(
      data_file_format = reactive(input$data_file_format),
      no_fdr = reactive(input$no_fdr),
      select_org = reactive(input$select_org),
      gmt_file = reactive(input$gmt_file),
      sample_info = reactive(loaded_data()$sample_info),
      all_gene_info = reactive(conversion_info()$all_gene_info),
      converted_data = reactive(conversion_info()$converted_data),
      all_gene_names = reactive(conversion_info()$all_gene_names),
      matched_ids = reactive(conversion_info()$converted$ids),
      gmt_choices = reactive(conversion_info()$gmt_choices),
      converted = reactive(conversion_info()$converted),
      no_id_conversion = reactive(input$no_id_conversion),
      heatmap_color_select = reactive(input$heatmap_color_select),
      select_gene_id = reactive(input$select_gene_id)
    )
  })
}

## To be copied in the UI
# mod_01_load_data_ui("load_data") # nolint

## To be copied in the server
# mod_01_load_data_server("load_data") # nolint
