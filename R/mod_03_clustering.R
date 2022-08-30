#' 03_heatmap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_03_clustering_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Clustering",
    sidebarLayout(

      # Heatmap Panel Sidebar ----------
      sidebarPanel(
        width = 3,
        conditionalPanel(
          condition = "input.cluster_panels == 'Hierarchical' | 
          input.cluster_panels == 'Gene SD Distribution' ",
          
          numericInput(
            inputId = ns("n_genes"), 
            label = h4("Top n most variable genes to include:"), 
            min = 10, 
            max = 12000, 
            value = 1000, 
            step = 10
          ), 
          ns = ns
        ),
        
        conditionalPanel(
          condition = "(input.cluster_panels == 'Hierarchical' | 
            input.cluster_panels == 'sample_tab') &&  input.cluster_meth == 2",
          
          # k- means slidebar -----------
            
            sliderInput(
              inputId = ns("k_clusters"),
              label = "Number of Clusters:",
              min   = 2,
              max   = 20,
              value = 4,
              step  = 1
            ),
            
            # Re-run k-means with a different seed
            actionButton(
              inputId = ns("k_means_re_run"),
              label = "Re-Run"
            ),
            
            # Elbow plot pop-up 
            actionButton(
              inputId = ns("elbow_pop_up"),
              label = "How many clusters?"
            ),
            # Line break ---------
            HTML(
              '<hr style="height:1px;border:none;
           color:#333;background-color:#333;" />'
            ),
          ns = ns
        ),

      

        # Select Clustering Method ----------
        conditionalPanel(
          condition = "input.cluster_panels == 'Hierarchical' | 
            input.cluster_panels == 'sample_tab'",
          
          selectInput(
            inputId = ns("cluster_meth"),
            label = "Select Clustering Method:",
            choices = list(
              "Hierarchical" = 1,
              "k-Means" = 2
            ),
            selected = 1
          ),

          ns = ns
        ),

        # Heatmap customizing features ----------
        conditionalPanel(
          condition = "input.cluster_panels == 'Hierarchical' ",
          
          # Gene ID Selection -----------
          selectInput(
            inputId = ns("select_gene_id"),
            label = "Gene ID type on Zoomed heatmap:",
            choices = NULL,
            selected = NULL
          ),
          
          # Sample coloring bar -----------
          htmlOutput(ns("list_factors_heatmap")),

          strong("Customize heatmap:"),
          fluidRow(
            br(),
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

          ns = ns
        ),

        # Clustering methods for hierarchical ----------
        conditionalPanel(
          condition = "input.cluster_meth == 1 && 
            (input.cluster_panels == 'Hierarchical' | 
            input.cluster_panels == 'sample_tab')",
          fluidRow(
            column(width = 4, h5("Distance")),
            column(
              width = 8,
              selectInput(
                inputId = ns("dist_function"),
                label = NULL,
                choices = NULL,
                width = "100%"
              )
            )
          ),
          fluidRow(
            column(width = 4, h5("Linkage")),
            column(
              width = 8,
              selectInput(
                inputId = ns("hclust_function"),
                label = NULL,
                choices = c(
                  "average", "complete", "single",
                  "median", "centroid", "mcquitty"
                ),
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
                value = 3,
                min = 2,
                step = 1
              )
            )
          ),
          ns = ns
        ),

        # Checkbox features ------------
        conditionalPanel(
          condition = "input.cluster_panels == 'Hierarchical' | 
            input.cluster_panels == 'sample_tab' ",

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
          ns = ns
        ),

        conditionalPanel(
          condition = "input.cluster_panels == 'Hierarchical' ",

          checkboxInput(
            inputId = ns("no_sample_clustering"),
            label = "Do not cluster samples",
            value = TRUE
          ),
          checkboxInput(
            inputId = ns("show_row_dend"),
            label = "Show Row Dendogram",
            value = TRUE
          ),
          br(),
          downloadButton(
            outputId = ns("download_heatmap_data"),
            label = "Heatmap data"
          ),
          ns = ns
        ),

        downloadButton(
          outputId = ns("report"),
          label = "Generate Report"
        ), 
        

        a(
          h5("Questions?", align = "right"),
          href = "https://idepsite.wordpress.com/heatmap/",
          target = "_blank"
        )
      ),





      #########################################################################
      # Main Panel
      #########################################################################

      mainPanel(
        tabsetPanel(
          id = ns("cluster_panels"),

          # Heatmap panel ----------
          tabPanel(
            title = "Hierarchical",
            br(),

            fluidRow(
              column(
                width = 4,
                plotOutput(
                  outputId = ns("heatmap_main"),
                  height = "500px",
                  width = "100%",
                  brush = ns("ht_brush")
                ),
                uiOutput(
                  outputId = ns("ht_click_content")
                )
              ),
              column(
                width = 8,
                plotOutput(
                  outputId = ns("sub_heatmap"),
                  height = "650px",
                  width = "100%",
                  click = ns("ht_click")
                )
              )
            ),
            checkboxInput(
              inputId = ns("cluster_enrichment"), 
              label = strong("Enrichment analysis on 
                selected genes or k-means clusters"),
              value = TRUE
            ),
            conditionalPanel(
              condition = "input.cluster_enrichment == 1 ",
              fluidRow(
                column(
                  width = 4,
                  htmlOutput(outputId = ns("select_go_selector"))
                ),
                column(
                  width = 4,
                  checkboxInput(
                    inputId = ns("filtered_background"),
                    label = "Use filtered genes as background.",
                    value = FALSE
                  )
                ),
                column(
                  width = 4,
                  checkboxInput(
                    inputId = ns("remove_redudant"),
                    label = "Remove Redudant Gene Sets",
                    value = FALSE
                  )
                ),
                tags$style(
                  type = 'text/css',
                  "#clustering-min_set_size {width:100%; margin-top:-12px}"
                ),
                tags$style(
                  type ='text/css',
                  "#clustering-max_set_size {width:100%; margin-top:-12px}"
                )
              ),
              mod_11_enrichment_ui(ns("enrichment_table_cluster")),
              ns = ns
            )

          ),

          # Gene Standard Deviation Distribution ----------
          tabPanel(
            title = "Gene SD Distribution",
            br(),
            plotOutput(
              outputId = ns("sd_density_plot"),
              width = "100%",
              height = "500px"
            ),
            ottoPlots::mod_download_figure_ui(ns("dl_gene_dist"))
          ),
          # Sample Tree -----------------
          tabPanel(
            title = "Sample Tree",
            value = "sample_tab", 
            h5(
              "Using genes with maximum expression level at the top 75%.
               Data is transformed and clustered as specified in the sidebar."
            ),
            br(),
            plotOutput(
              outputId = ns("sample_tree"),
              width = "100%",
              height = "400px"
            ),
            ottoPlots::mod_download_figure_ui(ns("dl_sample_tree"))
          )
        )
      )
    )
  )
}








#########################################################################
# Server function
#########################################################################

#' 03_heatmap Server Functions
#'
#' @noRd
mod_03_clustering_server <- function(id, pre_process, idep_data, tab) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Interactive heatmap environment
    shiny_env <- new.env()

    # Update Slider Input ---------
    observe({
      req(tab() == "Clustering")
      req(!is.null(pre_process$data()))
      if (nrow(pre_process$data()) > 12000) {
        max_genes <- 12000
      } else {
        max_genes <- round(nrow(pre_process$data()), -2)
      }
      updateNumericInput(
        inputId = "n_genes", 
        max = max_genes
      )
    })

    # Heatmap Colors ----------
    heatmap_colors <- list(
      "Green-Black-Red" = c("green", "black", "red"),
      "Red-Black-Green" = c("red", "black", "green"), 
      "Blue-White-Red" = c("blue", "white", "red"),
      "Green-Black-Magenta" = c("green", "black", "magenta"),
      "Blue-Yellow-Red" = c("blue", "yellow", "red"),
      "Blue-White-Brown" = c("blue", "white", "brown"), 
      "Orange-White-Blue" = c("orange", "white", "blue")
    )
    heatmap_choices <- c(
      "Green-Black-Red",
      "Red-Black-Green", 
      "Blue-White-Red",
      "Green-Black-Magenta",
      "Blue-Yellow-Red",
      "Blue-White-Brown", 
      "Orange-White-Blue"
    )
    observe({
      updateSelectInput(
        session = session,
        inputId = "heatmap_color_select",
        choices = heatmap_choices
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
        inputId = "dist_function",
        choices = dist_choices
      )
    })

    # Hclust Functions ----------
    hclust_funs <- hcluster_functions()

    # Gene ID Name Choices ----------
    observe({
      req(!is.null(pre_process$all_gene_names()))

      updateSelectInput(
        session = session,
        inputId = "select_gene_id",
        choices = colnames(pre_process$all_gene_names()),
        selected = "symbol"
      )
    })

    # Sample color bar selector ----------
    output$list_factors_heatmap <- renderUI({
      selectInput(
        inputId = ns("select_factors_heatmap"),
        label = "Sample Color Bar:",
        choices = c("Sample_Name", colnames(pre_process$sample_info()))
      )
    })

    # GMT choices for enrichment ----------
    output$select_go_selector <- renderUI({
	    req(!is.null(pre_process$gmt_choices()))
      selected <- "GOBP"
      if("KEGG" %in% pre_process$gmt_choices()) {
        selected <- "KEGG"
      }
	    selectInput(
        inputId = ns("select_go"),
        label = NULL,
        choices = pre_process$gmt_choices(),
        selected = selected
      )
    })

    # Standard Deviation Density Plot ----------
    sd_density_plot <- reactive({
      req(!is.null(pre_process$data()))
      
      sd_density(
        data = pre_process$data(),
        n_genes_max = input$n_genes
      )
    })
    
    output$sd_density_plot <- renderPlot({
      print(sd_density_plot())
    })
    
    dl_gene_dist <- ottoPlots::mod_download_figure_server(
      id = "dl_gene_dist", 
      filename = "sd_density_plot", 
      figure = reactive({ sd_density_plot() })
    )
    
    

    # Heatmap Data -----------
    heatmap_data <- reactive({
      req(!is.null(pre_process$data()))

      process_heatmap_data(
        data = pre_process$data(),
        n_genes_max = input$n_genes,
        gene_centering = input$gene_centering,
        gene_normalize = input$gene_normalize,
        sample_centering = FALSE,
        sample_normalize = FALSE,
        all_gene_names = pre_process$all_gene_names(),
        select_gene_id = input$select_gene_id
      )
    })

    # HEATMAP -----------
    # Information on interactivity
    # https://jokergoo.github.io/2020/05/15/interactive-complexheatmap/
    output$heatmap_main <- renderPlot({
      req(!is.null(heatmap_data()))
      req(!is.null(input$select_factors_heatmap))

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Creating Heatmap",
        color = "#000000"
      )

      # Assign heatmap to be used in multiple components
      shiny_env$ht <- heatmap_main(
        data = heatmap_data(),
        cluster_meth = input$cluster_meth,
        heatmap_cutoff = input$heatmap_cutoff,
        sample_info = pre_process$sample_info(),
        select_factors_heatmap = input$select_factors_heatmap,
        dist_funs = dist_funs,
        dist_function = input$dist_function,
        hclust_function = input$hclust_function,
        no_sample_clustering = input$no_sample_clustering,
        heatmap_color_select = heatmap_colors[[input$heatmap_color_select]],
        row_dend = input$show_row_dend,
        k_clusters = input$k_clusters,
        re_run = input$k_means_re_run
      )

      # Use heatmap position in multiple components
      shiny_env$ht_pos_main <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht)

      shinybusy::remove_modal_spinner()

      return(shiny_env$ht)
    })

    # Heatmap Click Value ---------
    output$ht_click_content <- renderUI({
      if (is.null(input$ht_click)) { 
        "Click on zoomed heatmap"
      } else {
        cluster_heat_click_info(
          click = input$ht_click,
          ht_sub = shiny_env$ht_sub,
          ht_sub_obj = shiny_env$ht_sub_obj,
          ht_pos_sub = shiny_env$ht_pos_sub,
          sub_groups = shiny_env$sub_groups,
          group_colors = shiny_env$group_colors,
          cluster_meth = input$cluster_meth,
          click_data = shiny_env$click_data
        )
      }
    })

    # Subheatmap creation ---------
    output$sub_heatmap <- renderPlot({
      if (is.null(input$ht_brush)) {
        grid::grid.newpage()
        grid::grid.text("Select a region on the heatmap to zoom in. 
        Gene IDs shows up when less than 60 genes are selected.", 0.5, 0.5)
      } else {
        submap_return <- heat_sub(
          ht_brush = input$ht_brush,
          ht = shiny_env$ht,
          ht_pos_main = shiny_env$ht_pos_main,
          heatmap_data = heatmap_data(),
          sample_info = pre_process$sample_info(),
          select_factors_heatmap = input$select_factors_heatmap,
          cluster_meth = input$cluster_meth
        )

        # Objects used in other components ----------
        shiny_env$ht_sub_obj <- submap_return$ht_select
        shiny_env$submap_data <- submap_return$submap_data
        shiny_env$sub_groups <- submap_return$sub_groups
        shiny_env$group_colors <- submap_return$group_colors
        shiny_env$click_data <- submap_return$click_data
        
        shiny_env$ht_sub <- ComplexHeatmap::draw(
          shiny_env$ht_sub_obj,
          annotation_legend_list = submap_return$lgd,
          annotation_legend_side = "top"
        )

        shiny_env$ht_pos_sub <- InteractiveComplexHeatmap::htPositionsOnDevice(shiny_env$ht_sub)

        return(shiny_env$ht_sub)
      }
    })

    # Enrichment Analysis ----------
    # Gene sets reactive
    pathway_table <- reactive({
      req(!is.null(input$select_gene_id))
      req(!is.null(input$ht_brush) || input$cluster_meth == 2)

      shinybusy::show_modal_spinner(
        spin = "orbit",
        text = "Running Analysis",
        color = "#000000"
      )
      
      pathway_info <- list()
      
      if (input$cluster_meth == 1) {
        gene_names <- merge_data(
          all_gene_names = pre_process$all_gene_names(),
          data = shiny_env$submap_data,
          merge_ID = input$select_gene_id
        )
      
        # Only keep the gene names and scrap the data
        gene_names_query <- dplyr::select_if(gene_names, is.character)

        req(!is.null(pre_process$all_gene_names()))
        req(!is.null(input$select_go))

        gene_sets <- read_pathway_sets(
          all_gene_names_query = gene_names_query,
          converted = pre_process$converted(),
          go = input$select_go,
          select_org = pre_process$select_org(),
          gmt_file = pre_process$gmt_file(),
          idep_data = idep_data,
          gene_info = pre_process$all_gene_info()
        )

        pathway_info[["Hierarchical_Selection"]] <- find_overlap(
          pathway_table = gene_sets$pathway_table,
          query_set = gene_sets$query_set,
          total_genes = gene_sets$total_genes,
          processed_data = pre_process$data(),
          gene_info = pre_process$all_gene_info(),
          go = input$select_go,
          idep_data = idep_data,
          select_org = pre_process$select_org(),
          sub_pathway_files = gene_sets$pathway_files,
          use_filtered_background = input$filtered_background,
          reduced = input$remove_redudant
        )
         # k-means-----------------------------------------------------
      } else if (input$cluster_meth == 2) {
        # Get the cluster number and Gene 
        
        req(heatmap_data())
        req(input$k_clusters)
        req(input$select_gene_id)
        req(shiny_env$ht)
        
        row_ord <- ComplexHeatmap::row_order(shiny_env$ht)
        
        req(!is.null(names(row_ord)))
        
        for (i in 1:length(row_ord)) {
          if (i == 1) {
          clusts <- data.frame(
            "cluster" = rep(names(row_ord[i]), length(row_ord[[i]])),
            "row_order" = row_ord[[i]]
          )
          } else {
          tem <- data.frame(
            "cluster" = rep(names(row_ord[i]), length(row_ord[[i]])),
            "row_order" = row_ord[[i]]
          )
          clusts <- rbind(clusts, tem)
          }
        }
        clusts$id <- rownames(heatmap_data()[clusts$row_order, ]) 

        # disregard user selection use clusters for enrichment
        for (i in 1:input$k_clusters) {
          cluster_data <- subset(clusts, cluster == i)
          row.names(cluster_data) <- cluster_data$id

          gene_names <- merge_data(
            all_gene_names = pre_process$all_gene_names(),
            data = cluster_data,
            merge_ID = input$select_gene_id
          )
      
          # Only keep the gene names and scrap the data
          gene_names_query <- dplyr::select_if(gene_names, is.character)

          req(!is.null(pre_process$all_gene_names()))
          req(!is.null(input$select_go))
          
          gene_sets <- read_pathway_sets(
            all_gene_names_query = gene_names_query,
            converted = pre_process$converted(),
            go = input$select_go,
            select_org = pre_process$select_org(),
            gmt_file = pre_process$gmt_file(),
            idep_data = idep_data,
            gene_info = pre_process$all_gene_info()
          )

          pathway_sub_info <- find_overlap(
            pathway_table = gene_sets$pathway_table,
            query_set = gene_sets$query_set,
            total_genes = gene_sets$total_genes,
            processed_data = pre_process$data(),
            gene_info = pre_process$all_gene_info(),
            go = input$select_go,
            idep_data = idep_data,
            select_org = pre_process$select_org(),
            sub_pathway_files = gene_sets$pathway_files,
          use_filtered_background = input$filtered_background,
          reduced = input$remove_redudant
          )

          pathway_info[[paste0("Cluster ", i)]] <- pathway_sub_info
        }
      }

      shinybusy::remove_modal_spinner()

      return(pathway_info)
    })

    # Sample Tree ----------
    sample_tree <- reactive({
      req(!is.null(pre_process$data()), input$cluster_meth == 1)
      
      draw_sample_tree(
        tree_data = pre_process$data(),
        gene_centering = input$gene_centering,
        gene_normalize = input$gene_normalize,
        sample_centering = FALSE,
        sample_normalize = FALSE,
        hclust_funs = hclust_funs,
        hclust_function = input$hclust_function,
        dist_funs = dist_funs,
        dist_function = input$dist_function
      )  
      p <- recordPlot() 
      return(p)
    })
    
    output$sample_tree <- renderPlot({
      print(sample_tree())
    })
    
    dl_sample_tree <- ottoPlots::mod_download_figure_server(
      id = "dl_sample_tree", 
      filename = "sample_tree", 
      figure = reactive({ sample_tree() })
    )
    
    observeEvent(input$cluster_meth, {
      if (input$cluster_meth == 1){
        showTab(
          inputId = "cluster_panels", 
          target = "sample_tab"
        ) 
      }
    })
    
    observeEvent(input$cluster_meth, {
      if(input$cluster_meth == 2){
        hideTab(
          inputId = "cluster_panels",
          target = "sample_tab"
        )
      }
    })
    
    # k-Cluster elbow plot ----------
    output$k_clusters <- renderPlot({
      req(!is.null(heatmap_data()))

      k_means_elbow(
        heatmap_data = heatmap_data()
      )
    })
    # pop-up modal 
    observeEvent(input$elbow_pop_up, {
      showModal(modalDialog(
        plotOutput(ns("k_clusters")), 
        footer = NULL, 
        easyClose = TRUE, 
        title = tags$h5(
          "Following the elbow method, one should choose k so that adding 
          another cluster does not substantially reduce the within groups sum of squares.",
          tags$a(
            "Wikipedia",
            href = "https://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set",
            target = "_blank"
          )
        ),  
      ))
    })

     # Heatmap Download Data -----------
    heatmap_data_download <- reactive({
      req(!is.null(pre_process$all_gene_names()))
      req(!is.null(heatmap_data()))

      merged_data <- merge_data(
        pre_process$all_gene_names(),
        heatmap_data(),
        merge_ID = input$select_gene_id
      )
    })

    output$download_heatmap_data <- downloadHandler(
      filename = function() {
        "heatmap_data.csv"
      },
      content = function(file) {
        write.csv(heatmap_data_download(), file)
      }
    )

  enrichment_table_cluster <- mod_11_enrichment_server(
    id = "enrichment_table_cluster",
    results = reactive({ pathway_table() }) # does not update?
  )
    
    # Markdown report------------
    output$report <- downloadHandler(
      
      # For PDF output, change this to "report.pdf"
      filename ="clustering_report.html",
      content = function(file) {
        #Show Loading popup
        shinybusy::show_modal_spinner(
          spin = "orbit",
          text = "Generating Report",
          color = "#000000"
        )
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "clustering_workflow.Rmd")
        #tempReport
        tempReport<-gsub("\\", "/",tempReport,fixed = TRUE)
        
        #This should retrieve the project location on your device:
        #"C:/Users/bdere/Documents/GitHub/idepGolem"
        wd <- getwd()
        
        markdown_location <-paste0(wd, "/vignettes/Reports/clustering_workflow.Rmd")
        file.copy(from=markdown_location,to = tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(
          pre_processed_data = pre_process$data(),
          sample_info = pre_process$sample_info(),
          all_gene_names = pre_process$all_gene_names(),
          n_genes = input$n_genes,
          k_clusters = input$k_clusters,
          cluster_meth = input$cluster_meth,
          select_gene_id = input$select_gene_id,
          list_factors_heatmap = input$list_factors_heatmap,
          heatmap_color_select = heatmap_colors[[input$heatmap_color_select]],
          dist_function = input$dist_function,
          hclust_function = input$hclust_function,
          heatmap_cutoff = input$heatmap_cutoff,
          gene_centering = input$gene_centering,
          gene_normalize = input$gene_normalize,
          no_sample_clustering = input$no_sample_clustering,
          show_row_dend = input$show_row_dend
            
          
        )
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(
          input = tempReport,#markdown_location, 
          output_file = file,
          params = params,
          envir = new.env(parent = globalenv())
        )
        shinybusy::remove_modal_spinner()
        
      }
      
    )
    

  })
}

## To be copied in the UI
# mod_03_heatmap_ui("03_heatmap_ui_1")

## To be copied in the server
# mod_03_heatmap_server("03_heatmap_ui_1")
