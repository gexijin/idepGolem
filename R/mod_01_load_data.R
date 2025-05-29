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
    # move notifications and progress bar to the center of screen
    tags$head(
      tags$style(
        HTML(".shiny-notification {
              width: 300px;
              position:fixed;
              top: calc(85%);
              left: calc(5%);
               }
              
             .dis_gray { background-color: gray; }
              "
        )
      )
    ),
    sidebarLayout(
      
      ##################################################################
      #       Load Data sidebar panel ----
      ##################################################################
      sidebarPanel(
        # alternative UI output message for once expression data is loaded
        uiOutput(ns("load_data_alt")),
        # Species Match Drop Down ------------
        selectInput(
          inputId = ns("select_org"),
          label = NULL,
          choices = setNames(c(99, 999), c("Human", "None")), 
          multiple = FALSE,
          selectize = TRUE,
          selected = setNames(999, "None")
        ),
        
        fluidRow(
          column(
            width = 4, 
            strong("1. Species"),
          ),
          column(
            width = 8,
            textOutput(ns("selected_species"))
          )
        ),
        tags$head(tags$style("#load_data-selected_species{color: blue;
                                 font-size: 12px;
                                 font-style: italic;
                                 }"
        )
        ),
        
        # .GMT file input bar ----------
        fluidRow(
          column(
            width = 6,
            align = "center",
            # Species list and genome assemblies ----------
            actionButton(
              inputId = ns("genome_assembl_button"),
              label = strong("Select")
            )
          ),
          column(
            width = 6,
            align = "center",
            # Species list and genome assemblies ----------
            actionButton(
              inputId = ns("upload_gmt_button"),
              label = "Custom"
            ),
            tippy::tippy_this(
              ns("upload_gmt_button"),
              "Upload a custom pathway .GMT file to perform pathway analysis for any genome.",
              theme = "light-border"
            )
          )
        ),
        
        br(),
        # Dropdown for data file format ----------
        strong("2. Data type"),
        
        selectInput(
          inputId = ns("data_file_format"),
          label = NULL,
          choices = list(
            "..." = 0,
            "Read counts data (recommended)" = 1,
            "Normalized Expression data" = 2,
            "Fold-changes & adjusted P-values" = 3
          ),
          selected = 0,
          selectize = FALSE
        ),
        tippy::tippy_this(
          ns("data_file_format"),
          "We recommend uploading Read Counts Data, which can be analyzed using DESeq2.
          Choose Normalized Expression Data if your data is derived from 
          RNA-Seq(FPKM, RPKM, TPM), DNA microarray, proteomics data, etc. 
          Select Fold Change and Adjusted P-value, if you have already conducted D.E.G. analysis.
            ",
          theme = "light-border"
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
        br(),
        # Load expression data options ----------
        # Includes load demo action button, demo data dropdown, and expression
        # file upload box
        conditionalPanel(
          condition = "input.data_file_format != 0 && input.select_org != 999",
          uiOutput(ns("load_data_ui")),
          ns = ns
        ),
        # tags$style(
        #   HTML("
        #     #load_data-ui {
        #       display: block !important;
        #     }
        #   ")
        # ),
        br(),
        
        # Experiment design file input ----------
        conditionalPanel(
          condition = "input.data_file_format != 0 && input.select_org != 999",
          uiOutput(ns("design_file_ui")),
          ns = ns
        ),
        uiOutput(ns("example_genes_ui")),
        br(),
        checkboxInput(
          inputId = ns("customize_button"),
          label = strong("Global Settings"),
          value = FALSE
        ),
        selectInput(
          inputId = ns("multiple_map"),
          label = "Multiple mapped IDs:",
          choices = list(
            "Sum" = "sum",
            "Average" = "mean",
            "Median" = "median",
            "Max" = "max",
            "Max SD" = "max_sd"
          ),
          selected = "Sum",
          selectize = FALSE
        ),
        tippy::tippy_this(
          ns("multiple_map"),
          "When multiple IDs map to the same gene, we can summerize
          data in a certain way (sum, mean, median, max),
          or just keep the rows with the the most variation (max SD).
          When uploading transcript level counts, choose \"sum\" to aggregate gene level counts. ",
          theme = "light-border"
        ),
        selectInput(
          inputId = ns("plots_color_select"),
          label = "Plots Color scheme:",
          choices = c(
            "Set1",
            "Set2",
            "Set3",
            "Paired",
            "Dark2",
            "Accent",
            "Pastel1",
            "Pastel2",
            "Spectral"
          ),
          selected = "Set1",
          width = "100%"
        ),
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
          label = "Gene ID type for plots:",
          choices = c("symbol", "ensembl_ID", "User_ID"),
          selected = "symbol"
        ),
        selectInput(
          inputId = ns("ggplot2_theme"),
          label = "ggplot2 theme:",
          choices = c(
            "default", # no change
            "gray",
            "bw",
            "light",
            "dark",
            "classic",
            "minimal",
            "linedraw"
          ),
          selected = "default",
          width = "100%",
          selectize = FALSE
        ),
        tippy::tippy_this(
          ns("ggplot2_theme"),
          "Changes the styles for all 20 ggplot2 plots.",
          theme = "light-border"
        ),
        checkboxInput(
          inputId = ns("plot_grid_lines"),
          label = "Add grid lines to plots",
          value = FALSE
        ),
        checkboxInput(
          inputId = ns("no_id_conversion"),
          label = "Do not convert gene IDs",
          value = FALSE
        ),
        tippy::tippy_this(
          ns("no_id_conversion"),
          "If selected, uploaded gene IDs will not be converted to ENSEMBL gene IDs,
          which is used as a central id type in pathway databases.",
          theme = "light-border"
        ),
        br(),
        fluidRow(
          column(
            width = 4,
            align = "left",
            # Link to public RNA-seq datasets ----------
            a(
              "Public Data",
              href = "http://bioinformatics.sdstate.edu/reads/"
            ),
          ),
          column(
            width = 4,
            align = "center",
            a(
              "Cite iDEP",
              href = "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6#citeas",
              target = "_blank"
            ),
          ),
          column(
            width = 4,
            align = "right",
            # Action button for Gene ID examples -----------
            a(
              "Questions?",
              href = "https://idepsite.wordpress.com/data-format/",
              target = "_blank"
            )
          )
        ),
        # Table output for species loading progress -----------
        #br(),
        #br(),
        #tableOutput(ns("species_match"))
      ),
      
      
      ##################################################################
      #       Load Data panel main ----
      ##################################################################
      mainPanel(
        shinyjs::useShinyjs(),
        # connection issue button
        #actionButton(ns("server_connection"), "Server Connection Tips"),
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
        uiOutput(ns("show_gmt")),        
        # Display file format help html document when prompted ----
        uiOutput(ns("format_help_ui")),
        # Hide welcome screen after data is loaded -----
        uiOutput(ns("welcome_ui"))
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
    
    observe({
      shinyjs::toggle(id = "plots_color_select", condition = input$customize_button)
      shinyjs::toggle(id = "heatmap_color_select", condition = input$customize_button)
      shinyjs::toggle(id = "select_gene_id", condition = input$customize_button)
      shinyjs::toggle(id = "multiple_map", condition = input$customize_button)
      shinyjs::toggle(
        id = "no_id_conversion", 
        condition = (input$customize_button && !grepl("STRING", input$clicked_row))
      )
      shinyjs::toggle(id = "plot_grid_lines", condition = input$customize_button)
      shinyjs::toggle(id = "ggplot2_theme", condition = input$customize_button)
    })


    welcome_modal <- shiny::modalDialog(
      title = "iDEP: Empower all scientists!",
  
      tags$br(),
      tags$p(" If iDEP is used,
      even for preliminrary analysis, please cite: ",
        a(
          " BMC Bioinformatics 19:1-24, 2018, ",
          href = "https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6",
          target = "_blank"
        ),
        "  which has been cited ",
        a("890 times.",
          href = "https://scholar.google.com/scholar?oi=bibs&hl=en&cites=6502699637682046008,17999801138713500070,11001860275874506471",
          target = "_blank"
        )
      ),
      tags$h5("By citing the iDEP paper properly, you will help make this service
      available in the future. Just including the URL is not enough.",
        style = "color:#6B1518"
      ),
      easyClose = TRUE,
      size = "l"
    )

    #shiny::showModal(welcome_modal)

    # Pop-up modal for gene assembl information ----
    observeEvent(input$genome_assembl_button, {
      shiny::showModal(
        shiny::modalDialog(
          size = "l",
          h3("Click on a row to select. Then close this window."),
          p("Search annotated species by common or scientific names,
          or NCBI taxonomy id. Click on a row to select. 
          Use ENSEMBL annotation if available. Use STRING-db annotation as a last resort.  
           If your species cannot be found here,
          you can still use iDEP without pathway analysis."),
          easyClose = TRUE,
          DT::renderDataTable({
            df <- idep_data$org_info[, 
              c("ensembl_dataset", "name", "academicName", "taxon_id", "group")
            ]
            colnames(df) <- c(
              "Ensembl/STRING-db ID",
              "Name (Assembly)",
              "Academic Name",
              "Taxonomy ID",
              "Source"
            )
            row.names(df) <- NULL
            DT::datatable(
              df,
              selection = "single",
              options = list(
                lengthChange = FALSE,
                pageLength = 10,
                scrollY = "400px"
              ),
              callback = DT::JS(
                paste0(
                 "table.on('click', 'tr', function() {
                    var data = table.row(this).data();
                    if (data) {
                      Shiny.setInputValue('", id, "-clicked_row', data[0]);
                    }
                  });"
                )
              ),
              rownames = FALSE
            )
          })
        )
      )
    })

    # Pop-up modal for uploading GMT file ----
    observeEvent(input$upload_gmt_button, {
      shiny::showModal(
        shiny::modalDialog(
          size = "l",
          p("Upload a custom pathway .GMT file to perform pathway analysis. 
          This enables you to analyze data for species not annotated in our database."),
          fileInput(
            inputId = ns("gmt_file"),
            label =
              "Upload a custom pathway .GMT file",
            accept = c(
              "text/csv",
              "text/comma-separated-values",
              "text/tab-separated-values",
              "text/plain",
              ".csv",
              ".tsv",
              ".gmt"
            ),
            placeholder = ""
          ),
          easyClose = TRUE
        )
      )
    })


    observeEvent(input$gmt_file, {
      req(!is.null(input$gmt_file))
      updateSelectizeInput(
        session = session,
        inputId = "select_org",
        choices = "NEW",
        selected = "NEW",
        server = TRUE
      )
      updateCheckboxInput(
        inputId = "no_id_conversion",
        label = "Do not convert gene IDs",
        value = TRUE
      )
    })

    output$show_gmt <- renderUI({
      req(!is.null(input$gmt_file))

      in_file <- input$gmt_file
      in_file <- in_file$datapath
      lines <- scan(in_file, what = "", sep = "\n")
      total <- length(lines)
      if(length(lines) > 3){
        lines <- lines[1:3]
      }
      tagList(
        h4("The uploaded GMT file has ", total, " gene-sets. 
        Data will be analyzed in a custom mode, not using our database.  
        The gene IDs in the GMT file must match those in the expression data."),
        br(),
        tags$pre(
          style = "font-size: 12px;",
          lines[1],
          ifelse(length(lines) > 1, lines[2], ""),
          ifelse(length(lines) > 2, lines[2], ""),
        )
      )


    })
    shinyjs::hideElement(id = "select_org")
    
    selected_species_name <- reactiveVal("None")
    
    observeEvent(input$clicked_row, {
      
      # find species ID from ensembl_dataset
      selected <- find_species_id_by_ensembl(
        input$clicked_row, 
        idep_data$org_info
      )
      # assign name
      selected <- setNames(
        selected,
        find_species_by_id_name(selected, idep_data$org_info)
      )

      updateSelectizeInput(
        session = session,
        inputId = "select_org",
        choices = selected,
        selected = selected,
        server = TRUE
      )

      selected_species_name(
        find_species_by_id_name(selected, idep_data$org_info)
      )
    })

    output$selected_species <- renderText({
      if(is.null(input$gmt_file)){
        selected_species_name() 
      } else {
        return("Custom")
      }
    })
    
    observeEvent(input$data_file_format, {
      
      req(input$data_file_format != 0)
      updateSelectInput(
        session = session,
        inputId = "data_file_format",
        choices = list("Read counts data (recommended)" = 1,
                       "Normalized Expression data" = 2,
                       "Fold-changes & adjusted P-values" = 3),
        selected = input$data_file_format
        )
    })

    # UI elements for load demo action button, demo data drop down, and -----
    # expression file upload
    output$load_data_ui <- renderUI({
      req(
        (is.null(input$go_button) || input$go_button == 0)
      )
      req(input$data_file_format)
      
      if (input$data_file_format > 0){
      # get demo data files based on specified format
      files <- idep_data$demo_file_info
      files <- files[files$type == input$data_file_format, ]
      choices <- setNames(as.list(files$ID), files$name)
      } else {
        files <- NULL
        choices <- NULL
      }
      tagList(
        strong("3. Expression matrix (CSV, text, or xlsx)"),
        fluidRow(
          column(
            width = 9,
            # Expression data file input
            fileInput(
              inputId = ns("expression_file"),
              label = NULL,
              accept = c(
                "text/csv",
                "text/comma-separated-values",
                "text/tab-separated-values",
                "text/plain",
                ".csv",
                ".tsv",
                ".xlsx",
                ".xls"
              ),
              placeholder = ""
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
        fluidRow(
          column(
            width = 4,
            actionButton(
              inputId = ns("go_button"),
              label = "Demo"
            ),
            align = "right",
            tippy::tippy_this(
              ns("go_button"),
              "Load the selected demo file",
              theme = "light-border"
            )
          ),
          column(
            width = 8,
            align = "left",
            selectInput(
              inputId = ns("select_demo"),
              label = NULL,
              choices = choices,
              selected = choices[[1]],
              selectize = FALSE
            ),
            tippy::tippy_this(
              ns("select_demo"),
              "Select a demo file then click the \"Demo\" button",
              theme = "light-border"
            ),
            tags$style(
              type = "text/css",
              "#load_data-go_button { margin-top:-25px; color: red;}"
            ),
            tags$style(
              type = "text/css",
              "#load_data-select_demo { margin-top:-20px}"
            )
          )
        )
      )
    })

    # Alternate ui message to reset app once data is loaded ----
    output$load_data_alt <- renderUI({
      if (
        is.null(input$go_button) || input$go_button == 0 && is.null(input$expression_file)

      ) {
        tagList(
          p(
            "Watch a ",
            a(
              "video.",
              href = "https://youtu.be/Hs5SamHHG9s",
              target = "_blank"
            ),
            " Or click ",
            tags$span("Demo", id = "load-demo"),
            " below to try!"
          ),
          tags$script("
            document.getElementById('load-demo').style.color = 'red';
          ")
        )
      } else {
        # reset message and action button
        tagList(
          fluidRow(
            column(
              width = 12,
              align = "right",
              actionButton(
                inputId = ns("reset_app_new_data"),
                label = strong("Reset"),
                align = "right"
              )
            )
          ),
          br(),
          br()
        )
      }
    })

    observeEvent(input$reset_app_new_data, {
      session$reload()
    })

    # UI element for design file upload ----
    output$design_file_ui <- renderUI({
      req(is.null(input$go_button) || input$go_button == 0)
      
      tagList(

        strong("4. Optional: Exp. Design (CSV or text)"),
        fluidRow(
          column(
            width = 9,
            fileInput(
              inputId = ns("experiment_file"),
              label = NULL,
              accept = c(
                "text/csv",
                "text/comma-separated-values",
                "text/tab-separated-values",
                "text/plain",
                ".csv",
                ".tsv",
                ".xlsx",
                ".xls"
              ),
              placeholder = ""
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
        )
      )
    })
    
    #name variables like this: gene_ids_example_popup, not GeneIDsExamplePopup.
    output$example_genes_ui <- renderUI({
      actionButton(ns("gene_ids_example_popup"), "Gene IDs")
    })


    # Disables expression_file input to prevent multiple uploads
    observeEvent(input$expression_file, {
      shinyjs::disable("expression_file")
    })

    # Define the content of the basic modal
    observeEvent(input$gene_ids_example_popup, {
      shiny::showModal(
        shiny::modalDialog(
          title = "What do the gene IDs in our database look like?",
          tags$style(
            HTML(
              "#DataTables_Table_0_wrapper #DataTables_Table_0_filter label{
                width: 400px;
                float: left;
              }"
            )
          ),
          selectizeInput(
            inputId = ns("gene_id_examples"),
            label = "Select or search for species",
            choices = c("--Select species--", names(idep_data$species_choice))
          ),
          DT::dataTableOutput(ns("showGeneIDs4Species")),
          size = "l", # size is large
          easyClose = FALSE   # disabled: click outside the modal to close
        )
      )
    })
    
    geneIDs <- reactiveVal(NULL)
    
    observeEvent(input$gene_id_examples, {
      req(input$gene_id_examples != "--Select species--")
      ix <- which(idep_data$org_info$name2 == input$gene_id_examples)
      dbase <- idep_data$org_info$file[ix]
      geneIDs(
        showGeneIDs(
          species = input$gene_id_examples,
          db = dbase,
          nGenes = 10
        )
      )
    })
    
    # Render Gene ID example table in gene example Modal
    output$showGeneIDs4Species <- DT::renderDataTable({
      req(!is.null(geneIDs()))
      req(input$gene_id_examples != "--Select species--")
      removeNotification("ExampleIDDataQuery")
      
      DT::datatable(
        geneIDs(),
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        )
      )
    })

    # Disables experiment_file input to prevent multiple uploads
    observeEvent(input$experiment_file, {
      shinyjs::disable("experiment_file")
    })
    
    # Notification for data type and species selection
    observe({
      req(tab() == "Load Data")
      req(input$select_org == 999 || input$data_file_format == 0)
        
        showNotification(
          "Select a species and a data type before uploading data.",
          duration = 30,
          type = "error",
          id = "select_first"
          )
    })
    
    observe({
      req((input$select_org != 999 && input$data_file_format != 0) ||
          (tab() != "Load Data"))
      
      removeNotification("select_first")
      
    })

    # Show messages when on the Network tab or button is clicked ----
    observe({
      req(is.null(loaded_data()$data) && (
        tab() != "Load Data" || tab() != "About"
      ))

      showNotification(
        ui = paste("Please load a demo file or your own data first."),
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

    observe({
      req(!is.null(loaded_data()$data) && any(apply(loaded_data()$data, 2, function(col) all(col == 0))))

      showNotification(
        ui = paste("A sample has all values as zero. it is recommended to remove that sample."),
        id = "sample_remove_error",
        duration = NULL,
        type = "error"
      )
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
        #conversion_info()$converted_data[1:20, ],
        if(nrow(loaded_data()$data) > 20){
          loaded_data()$data[1:20, ]
        } else (loaded_data()$data),
        options = list(
          pageLength = 10,
          scrollX = "400px",
          dom = "t"
        ),
        rownames = TRUE
      )
    })

    observeEvent(input$expression_file, {
      # test data for correct format
      if (
        min(loaded_data()$data, na.rm = TRUE) < 0 && input$data_file_format == 1
      ) {
        showModal(modalDialog(
          title = "Data type does not match data format",
          tags$p("Negative values were detected in this dataset. This is not
                 correct for the selected data type (Read Counts).
                 Please double check data type or the data."),
          tags$br(),
          size = "m",
          easyClose = TRUE
          #footer = actionButton(ns("reset_app"), "Start over")
        ))
      }
    })

    # Get converted IDs ----------
    conversion_info <- reactive({
      req(!is.null(loaded_data()$data))


      req(!(min(loaded_data()$data, na.rm = TRUE) < 0 && input$data_file_format == 1))

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
          data = loaded_data()$data,
          multiple_map = input$multiple_map
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
      
    })

    observeEvent(input$reset_app, {
      session$reload()
    })

    # download database for selected species
    observeEvent(input$select_org, {
      ix <- which(idep_data$org_info$id == input$select_org)
      db_file <- idep_data$org_info[ix, "file"]
      dbname <- paste0(DATAPATH, "db/", db_file)
      if (!file.exists(dbname)) {
        withProgress(
          message = paste(
            "Download database for",
            idep_data$org_info[ix, "name2"],
           "(~5 minutes)"
           ), {
          incProgress(0.2)
          # download org_info and demo files to current folder
          options(timeout = 3000)
          download.file(
            url = paste0(db_url, db_ver, "/db/", db_file, ".gz"),
            destfile = paste0(dbname, ".gz"),
            mode = "wb",
            quiet = FALSE
          )
          incProgress(0.7)
          R.utils::gunzip(
            paste0(dbname, ".gz"), 
            remove = TRUE
          ) # untar and unzip the files
        })

      }
    })

    # Species match table ----------
    output$species_match <- renderTable({
      species_match_data()
      },
      digits = -1,
      spacing = "s",
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = TRUE
    )

    species_match_data <- reactive({
      req(!is.null(input$go_button))
      req(input$select_org)
      if (is.null(input$expression_file) && input$go_button == 0) {
        return(NULL)
      }
      req(input$data_file_format)
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
    })
    # Species match message ----------
    observe({
      req(
        tab() == "Load Data" &&
          #!is.null(conversion_info()$converted)
          species_match_data()[1,1] == "ID not recognized." &&
          input$select_org != "NEW"
      )


      showModal(modalDialog(
        title = "Please double check the selected species",
        tags$p("None of the gene IDs are recognzied. Possible causes: 1. Wrong species is selected. 
        2. Correct species is selected but we cannot map your gene IDs to Ensembl gene IDs. 
        3. Your species is not included in our database.  
        You can still run many analyses except pathway and enrichment."),
        size = "s",
        easyClose = TRUE
      ))
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
            width = 9,
            h4("iDEP: integrated Differential Expression & Pathway analysis")
          ),
          column(
            width = 3,
            img(
              src = "www/idep_logo.png",
              width = "43",
              height = "50"
            )
          )
        ),
        htmlOutput(ns("file_format")),
        includeHTML(app_sys("app/www/messages.html")),
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

      includeHTML(app_sys("app/www/format.html"))
    })

observeEvent(input$server_connection, {
  showModal(
    modalDialog(
      title = "Server Connection Tips",
      p("The iDEP webserver uses a load balancer to distribute incoming traffic across 125 virtual 
      servers, which are hosted on three physical servers(SDSU, Azure, & JetStream2). Each virtual server is capable of 
      serving multiple users simultaneously, provided it is not operating at full CPU capacity. 
      The assignment of users to specific virtual servers is determined by 
      their browser session and IP address."),
      tags$h4("Slow Performance:"),
      tags$ul(
        tags$li("Trying using iDEP in a new browser window (not just a new tab)."),
        tags$li("Using a different browser (Chrome, Edge, Safari, Firefox, etc.) to be assigned a different virtual server.")
      ),
      
      tags$h4("Connection Issues:"),
      tags$p("If you can't connect, try accessing iDEP directly using these mirror server URLs like this:"),
      tags$ul(
        tags$li(tags$a(href="http://149.165.173.123:55011/idep/", target="_blank", "JetStream2 (http://149.165.173.123:55011/idep/)")),
        tags$li(tags$a(href="http://4.236.179.243:55011/idep/", target="_blank", "Azure (http://4.236.179.243:55011/idep/)"))
      ),
      tags$p("Note: The port number (55011 in this example) can vary between 55001 and 55050, each pointing to a different virtual server."),
      
      tags$h4("Frequent Crashing:"),
      tags$p("There are several reasons iDEP might crash:"),
      tags$ul(
        tags$li("Problems with your data."),
        tags$li("Running a large analysis."),
        tags$li("The assigned virtual server overload (100% CPU).")
      ),
      tags$p("Try from a new browser window. Try again later. If issues persist, email us."),
      
      tags$h4("Insecure Connection (http):"),
      tags$p("We are transitioning to a secure https protocol. For now, you can:"),
      tags$ul(
        tags$li("Trust the current http connection."),
        tags$li("Download and run iDEP on your laptop.")
      ),
      easyClose = TRUE
    )
  )
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
      plots_color_select = reactive(input$plots_color_select),
      heatmap_color_select = reactive(input$heatmap_color_select),
      select_gene_id = reactive(input$select_gene_id),
      multiple_map = reactive(input$multiple_map),
      plot_grid_lines = reactive(input$plot_grid_lines),
      ggplot2_theme = reactive(input$ggplot2_theme)
    )
  })
}

## To be copied in the UI
# mod_01_load_data_ui("load_data") # nolint

## To be copied in the server
# mod_01_load_data_server("load_data") # nolint

