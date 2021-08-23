#' 01_load_data UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
### testing something, come back to later
mod_01_load_data_ui <- function(id) {
  ns <- shiny::NS(id)
  tabPanel("Load Data",
    sidebarLayout(

      # Load Data Panel Sidebar -----------
      sidebarPanel(

        # Button to load demo dataset ----------
        # Manually namespace the goButton in tag with id in module call
        actionButton(
          inputId = ns("go_button"),
          label = "Click here to load demo data"
        ),
        tags$head(tags$style(
          "#load_data-go_button{color: red;
          font-size: 16px;
          font-style: italic;}"
          )
        ),
        h5(" and just click the tabs for some magic!", style = "color:red"),

        # Reset Button -----------
        p(HTML(
          "<div align=\"right\"><A HREF=\"javascript:history.go(0)\"
           >Reset</A></div>"
          )
        ),

        # Species Match Drop Down ------------
        strong("1. Select or search for your species."),
        selectizeInput(
          inputId = ns("select_org"),
          label = NULL,
          choices = " ",
          multiple = TRUE,
          options = list(
            maxItems = 1,
            placeholder = "Best matching species",
            onInitialize = I('function() { this.setValue(""); }')
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
        radioButtons(
          inputId = ns("data_file_format"),
          label = "2. Choose data type",
          choices = list(
            "Read counts data (recommended)" = 1,
            "Normalized expression values (RNA-seq FPKM, microarray, etc.)" = 2,
            "Fold-changes and corrected P values from CuffDiff or any other
             program" = 3
          ),
          selected = 1
        ),

        # Conditional panel for fold changes data file ----------
        conditionalPanel(
          condition = "input.data_file_format == 3",
          checkboxInput(
            inputId = ns("no_fdr"),
            label = "Fold-changes only, no corrected P values",
            value = FALSE
          ),
          ns = ns
        ),

        # Expression data file input ----------
        fileInput(
          inputId = ns("expression_file"),
          label = "3. Upload expression data (CSV or text)",
          accept = c(
            "text/csv",
            "text/comma-separated-values",
            "text/tab-separated-values",
            "text/plain",
            ".csv",
            ".tsv"
          )
        ),

        # Link to public RNA-seq datasets ----------
        a(
          h4("Analyze public RNA-seq datasets for 9 species"),
          href = "http://bioinformatics.sdstate.edu/reads/"
        ),

        # Experiment design file input ----------
        fileInput(
          inputId = ns("experiment_file"),
          label = h5("Optional: Upload an experiment design file(CSV or text)"),
          accept = c(
            "text/csv",
            "text/comma-separated-values",
            "text/tab-separated-values",
            "text/plain",
            ".csv",
            ".tsv"
          )
        ),

        # Table output for species loading progress -----------
        tableOutput(ns("species")),

        # Action button for Gene ID examples -----------
        h5(
          "Check this out if you want example of our gene ids, or download gene
          mapping."
        ),

        # ADD GENE ID EXAMPLE CODE FOR BUTTON ----------
        actionButton(
          inputId = ns("gene_id_button"),
          label =  "Optional: COMING SOON!"
        ),
        ##################################################

        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/data-format/",
          target = "_blank"
        )


      ),

      # Load Data panel main -----------
      mainPanel(

        shinyjs::useShinyjs(),

        # Table output for sample tissue type ----------
        tableOutput(ns("sample_info_table")),

        # Display first 20 rows of the data ----------
        tableOutput(ns("sample_20")),

        # Instructions and flowchart ------------
        div(
          id = "load_message",
          h4("Loading R packages, please wait ... ... ...")
        ),
        htmlOutput("file_format"),
        h3(
          "We found an issue with the Gene Onotology database derived from
           Ensembl Release 103, which is used in iDEP 0.93. While we are fixing
           this issue, we have reverted the database to a previous version used
           in iDEP 0.92. "
        ),
        h4("Postdoc and GRA positions available!"),
        h4(
          "If your gene IDs are not recognized, please let us know. We might be
           able to add customized gene mappings to Ensembl gene IDs."
        ),
        h3(
          "New version 0.93 released on 5/23/2021 includes upgrades to R 4.05,
           Bioconductor 3.12, larger database (5000+ species) from Ensembl
           Release 103 and STRING-db v11. Massive, manually-collected pathway
           database for 20 model organisms. Fixed KEGG pathway chart and gene
           plot.", style = "color:red"
        ),
        h4(
          "We recently hired Jenny Qi for database updates and user support.",
          a(
            "Email Jenny for questions.",
            href = "mailto:gelabinfo@gmail.com?Subject=iDEP"
          )
        ),
        h5(
          "iDEP has not been thoroughly tested. Please let us know if you find
           any issue/bug."
        ),
        h5("We will be happy to help prepare your data for iDEP."),
        br(),
        img(
          src = "www/flowchart.png",
          align = "center",
          width = "562",
          height = "383"),

      )
    )
  )
}

#' 01_load_data Server Functions
#'
#' @noRd
### testing something, come back to later
mod_01_load_data_server <- function(id, idep_data, pre_process) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Provide species list for dropdown selection -----------
    observe({
      updateSelectizeInput(
        session = session,
        inputId = "select_org",
        choices = idep_data$species_choice,
        selected = idep_data$species_choice[1]
      )
    })

    # Sample information table -----------
    output$sample_info_table <- renderTable({
      if (is.null(read_sample_info())) {
        return(NULL)
      }
      isolate({
        tem <- t(read_sample_info())
        tem <- cbind(rownames(tem), tem)
        colnames(tem)[1] <- "Study_design"
        return(tem)
      })
      },
      include.rownames = FALSE,
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = TRUE
    )

    # First 20 rows of dataset table -----------
    output$sample_20 <- renderTable({
      in_file <- input$expression_file
      in_file <- in_file$datapath
      if (is.null(input$expression_file) && input$go_button == 0) {
        return(NULL)
      }
      if (is.null(input$expression_file) && input$go_button > 0) {
        in_file <- idep_data$demo_data_file
      }
      tem <- input$select_org
      isolate({
        x <- read.csv(in_file)
        if (dim(x)[2] <= 2 ) {
          x <- read.table(
            in_file,
            sep = "\t",
            header = TRUE)
        }
        x[1:20, ]
      })
      },
      include.rownames = FALSE,
      striped = TRUE,
      bordered = TRUE,
      width = "auto",
      hover = TRUE
    )

    read_data <- shiny::reactive({
      kurtosis_log <- 50

      in_file <- input$expression_file
      in_file <- in_file$datapath

      if (is.null(input$expression_file) && input$go_button == 0) {
        return(NULL)
      }

      if (is.null(input$expression_file) && input$go_button > 0) {
        in_file <- idep_data$demo_data_file
      }

      tem = input$data_file_format
      tem = pre_process$missing_value()

      # These are needed to make it responsive to changes
      if (!is.null(input$data_file_format))
        if (input$data_file_format == 1) {
          tem = pre_process$min_counts()
          tem = pre_process$n_min_samples_count()
          tem = pre_process$counts_log_start()
          tem = pre_process$counts_transform()
        }

      if (!is.null(input$data_file_format))
        if (input$data_file_format == 2) {
          tem = pre_process$log_transform_fpkm();
          tem = pre_process$log_start_fpkm();
          tem = pre_process$low_filter_fpkm()
          tem = pre_process$n_min_samples_fpkm()
        }

      isolate({
        withProgress(message = "Reading and pre-processing ", {
          # these packages moved here to reduce loading time
          library(edgeR, verbose = FALSE) # count data D.E.
          library(DESeq2, verbose = FALSE) # count data analysis

          if (is.null(input$data_file_format)) {
            return(NULL)
          }
          data_type_warning = 0
          data_type = c(TRUE)

          # Read file ----------
          # CSV attempt
          data <- read.csv(in_file, quote = "", comment.char = "")	

          # Tab-delimented if not CSV
          if (dim(data)[2] <= 2) {
            data <- read.table(
              in_file,
              sep = "\t",
              header = TRUE,
              quote = "",
              comment.char = ""
            )
          }

          # Remove non-numeric vars, except the first column -------
          for (i in 2:dim(data)[2]) {
            data_type <- c(data_type, is.numeric(data[, i]))
          }
          if (sum(data_type) <= 2) {
            return(NULL)
          }
          data <- data[, data_type]  # only keep numeric columns

          # rows with all missing values
          missing_filter = which(
            apply(data[, -1],
            1,
            function(y)
            sum(is.na(y))) != dim(data)[2] - 1
          )
          data <- data[missing_filter, ]

          data_size_original <- dim(data)
          data_size_original[2] <- data_size_original[2] - 1
          data[, 1] <- toupper(data[, 1])
          data[, 1] <- gsub(" |\"|\'", "", data[, 1])
          # remove spaces in gene ids
          # remove " in gene ids, mess up SQL query
          # remove ' in gene ids
          # remove one or two digits after "." at the end.
          # A35244.1 -> A35244  or A35244.23 -> A35244, but not more than two.
          # GLYMA.18G52160 stays the same.

          # Sort by standard deviation -----------
          data = data[order(-apply(
            data[, 2:dim(data)[2]],
            1,
            sd
            )
            )
          , ]

          # Remove duplicated genes ----------
          data <- data[!duplicated(data[, 1]), ]
          rownames(data) <- data[, 1]
          data <- as.matrix(data[, c(-1)])

          # remove "-" or "." from sample names ----------
          colnames(data) <- gsub("-", "", colnames(data))
          colnames(data) <- gsub("\\.", "", colnames(data))

          # Missng value in data ----------
          if (sum(is.na(data)) > 0) {
            if (pre_process$missing_value() == "geneMedian") {
              row_medians <- apply(data, 1, function(y)  median(y, na.rm = T))
              for(i in 1:dim(data)[2]) {
                val_miss_row = which(is.na(data[, i]))
                data[val_miss_row, i] <- row_medians[val_miss_row]
              }
            } else if (pre_process$missing_value() == "treatAsZero") {
              data[is.na(data)] <- 0
            } else if (pre_process$missing_value() == "geneMedianInGroup") {
              sample_groups = detect_groups(colnames(data))
              for (group in unique(sample_groups)) {
                samples = which(sample_groups == group)
                row_medians <- apply(
                  data[,samples, drop = F],
                  1,
                  function(y)  median(y, na.rm = T)
                )
                for (i in  samples) {
                  missing = which(is.na(data[, i]))
                  if(length(mssing) > 0)
                  data[missing, i]  <- row_medians[misssing]
                  }
              }
              if (sum(is.na(data)) > 0) {
                row_medians <- apply(
                  data,
                  1,
                  function(y)  median(y, na.rm = T)
                )
                for(i in 1:dim(data)[2]) {
                  missing = which(is.na(data[,i]))
                  data[missing, i] <- row_medians[missing]
                }
              }
            }
          }

          # Compute kurtosis ---------
          mean_kurtosis <- mean(apply(data, 2, e1071::kurtosis), na.rm = T)
          raw_counts <- NULL
          pvals <- NULL
          if (input$data_file_format == 2) {
            incProgress(1 / 3, "Pre-processing data")
            if (is.integer(data)) data_type_warning = 1;

            # Filters ----------

            # Not enough counts
            data <- data[which(apply(
              data,
              1,
              function(y) sum(y >= pre_process$low_filter_fpkm)
              ) >= input$n_min_samples_fpkm), ]

            # Same levels in every entry
            data <- data[which(apply(
              data,
              1,
              function(y) max(y)- min(y)) > 0), ]

            # Takes log if log is selected OR kurtosis is bigger than 50
            if (
              (pre_process$log_transform_fpkm == TRUE) |
              (mean_kurtosis > kurtosis_log)) {
                data <- log(data + abs(pre_process$log_start_fpkm()), 2)
            }

            std_dev <- apply(data, 1, sd)
            data <- data[order(-tem), ]

          } else if (input$data_file_format == 1) {
            incProgress(1 / 3, "Pre-processing counts data")
            tem = pre_process$counts_deg_method();
            tem = pre_process$counts_transform()

            if (!is.integer(data) & mean_kurtosis < kurtosis_log) {
              data_type_warning = -1
            }

            data <- round(data, 0)

            data <- data[which(apply(
              edgeR::cpm(edgeR::DGEList(counts = data)),
              1,
              function(y) sum(y >= pre_process$min_counts())) >=
              pre_process$n_min_samples_count()), ]

            raw_counts = data; # ???
            browser()

            # Construct DESeqExpression Object
            tem = rep("A", dim(data)[2]); tem[1] <- "B"
            col_data = cbind(colnames(data), tem)
            colnames(col_data) <- c("sample", "groups")
            dds <- DESeq2::DESeqDataSetFromMatrix(
              countData = data,
              colData = col_data,
              design = ~ groups
            )
            dds <- DESeq2::estimateSizeFactors(dds)

            incProgress(1 / 2, "transforming raw counts")

            # Counts Transformation ------------
            if(pre_process$counts_transform() == 3) {
              {rlog_data <- DESeq2::rlog(dds, blind = TRUE);
              rlog_data <- assay(rlog_data)}
            } else {
              if (pre_process$counts_transform() == 2) {
                vst_data <- vst(dds, blind=TRUE)
                vst_data <- assay(vst_data)
              } else {
                log_2_data <- log2(BiocGenerics::counts(
                  dds,
                  normalized = TRUE
                  ) + pre_process$counts_log_start()
                )
              }
            }
          } else if (input$data_file_format == 3) {
            n2 = (dim(data)[2] %/% 2)
            if (!input$no_fdr) {
              pvals = data[, 2 * (1:n2), drop = FALSE]
              data = data[, 2 * (1:n2) - 1, drop = FALSE]
              if (dim(data)[2] == 1) {
                placeholder <- rep(1, dim(data)[1])
                pvals <- cbind(pvals, placeholder)
                zero_placeholder <- rep(0, dim(data)[1])
                data <- cbind(data, zero_placeholder)
              }
            }
          }

          data_size = dim(data);

          sample_choice <- stats::setNames(
            as.list(1:(dim(data)[2])), colnames(data)
          )

          # observe({updateSelectInput(session, "scatterX", choices = sampleChoice, selected = sampleChoice[1]) })
          # bserve({updateSelectInput(session, "scatterY", choices = sampleChoice, selected = sampleChoice[2]) })
          validate(
            need(
              dim(data)[1] > 5 & dim(data)[2] >= 1, 
              "Data file not recognized. Please double check."
            )
          )

          incProgress(1, "Done.")

          sample_info_demo = NULL

          if (input$go_button > 0) {
            sample_info_demo <- t(read.csv(
              idep_data$demo_metadata_file,
              row.names = 1,
              header = T,
              colClasses = "character"
            ))
          }

          final_result <- list(
            data = as.matrix(data),
            mean_kurtosis = mean_kurtosis,
            raw_counts = raw_counts,
            data_type_warning = data_type_warning,
            data_size = c(data_size_original, data_size),
            sample_info_demo = sample_info_demo,
            pvals = pvals
          )

        return(final_result)
        })
      })
    })

    read_sample_info <- reactive ({
      if (is.null(input$expression_file) &&
      !is.null(read_data()$sample_info_demo)) {
        return(read_data()$sample_info_demo)
      }
      in_file <- input$experiment_file
      in_file <- in_file$datapath

      if (is.null(input$experiment_file) && input$go_button == 0) {
        return(NULL)
      }
      if (is.null(read_data())) {
        return(NULL)
      }

      isolate({
        if (is.null(input$data_file_format)) {
          return(NULL)
        }

        data_type_warning <- 0
        data_type <- c(TRUE)

        # Read experiment file ----------
        expr <- read.csv(
          in_file,
          row.names = 1,
          header = TRUE,
          colClasses = "character"
        )
        if (dim(expr)[2] <= 2) {
          expr <- read.table(
            in_file,
            row.names = 1,
            sep = "\t",
            header = TRUE,
            colClasses = "character"
          )
        }

        # remove "-" or "." from sample names ----------
        colnames(expr) = gsub("-", "", colnames(expr))
        colnames(expr) = gsub("\\.", "", colnames(expr))

        # Matching with column names of expression file ----------
        matches = match(
          toupper(colnames(read_data()$data)), toupper(colnames(expr))
        )
        matches = matches[which(!is.na(matches))] # remove NA

        validate(need(
          length(unique(matches)) == dim(read_data()$data)[2] &
                 dim(expr)[1] >= 1 & dim(expr)[1] < 500,
          "Error!!! Sample information file not recognized. Sample names
           must be exactly the same. Each row is a factor. Each column
           represent a sample.  Please see documentation on format."
          )
        )

        # Check factor levels, change if needed ----------
        for (i in 1:dim(expr)[1]) {
          expr[i, ] = gsub("-", "", expr[i, ])
          expr[i, ] = gsub("\\.", "", expr[i,])
        }

        # Factor levels match ----------
        if (length(unique(matches)) == dim(read_data()$data)[2]) {
          expr = expr[, matches]

          if(
            sum(apply(expr, 1, function(y) length(unique(y)))) >
            length(unique(unlist(expr)))) {
              factor_names = apply(
                expr,
                2,
                function(y) paste0(names(y), y)
              )
              rownames(factor_names) = rownames(expr)
              expr <- factor_names
            }
          return(t(expr))
        } else {
          return(NULL)
        }
      })
    })


    list(
      data_file_format = reactive(input$data_file_format)
    )
  })
}

## To be copied in the UI
# mod_01_load_data_ui("load_data") # nolint

## To be copied in the server
# mod_01_load_data_server("load_data") # nolint
