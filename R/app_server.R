#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # file size is 5MB by default. This changes it to 30MB
  # options(shiny.maxRequestSize = 30*1024^2)
  options(warn = -1) # turn off warning
  pdf(file = NULL)

  # define where database is located
  db_ver <<- "data113"
  db_url <<- "http://bioinformatics.sdstate.edu/data/"

  # if environmental variable is not set, use relative path
  DATAPATH <<- Sys.getenv("IDEP_DATABASE")[1]
  # if not defined in the environment, use too levels above
  if (nchar(DATAPATH) == 0) {
    DATAPATH <<- paste0("../../data/")
  }
  #Add version
  DATAPATH <<- paste0(DATAPATH, "/", db_ver, "/")
  org_info_file <<- paste0(DATAPATH, "demo/orgInfo.db")
  if(!file.exists(org_info_file)) {
    DATAPATH <<- paste0("./", db_ver, "/")
    org_info_file <<- paste0(DATAPATH, "demo/orgInfo.db")
  }

  # load static data files such as list of species, gmt files, etc
  # This could be moved to run_app as global variable, as in global.R
  # see https://github.com/ThinkR-open/golem/issues/6
  idep_data <- get_idep_data()

  # Tab Variable to control reactivity
  tab <- reactive(input$navbar)

  load_data <- mod_01_load_data_server(
    id = "load_data",
    idep_data = idep_data,
    tab = tab
  )
  pre_process <- mod_02_pre_process_server(
    id = "pre_process",
    load_data = load_data,
    tab = tab
  )
  mod_03_clustering_server(
    id = "clustering",
    pre_process = pre_process,
    load_data = load_data,
    idep_data = idep_data,
    tab = tab
  )
  mod_04_pca_server(
    id = "pca",
    load_data = load_data,
    pre_process = pre_process,
    idep_data = idep_data
  )
  deg <- mod_05_deg_server(
    id = "deg",
    pre_process = pre_process,
    idep_data = idep_data,
    load_data = load_data,
    tab = tab
  )
  mod_06_pathway_server(
    id = "pathway",
    pre_process = pre_process,
    deg = deg,
    idep_data = idep_data,
    tab = tab
  )
  mod_07_genome_server(
    id = "genome",
    pre_process = pre_process,
    deg = deg,
    idep_data = idep_data
  )
  mod_08_bicluster_server(
    id = "bicluster",
    pre_process = pre_process,
    idep_data = idep_data,
    tab = tab
  )
  mod_09_network_server(
    id = "network",
    pre_process = pre_process,
    idep_data = idep_data,
    tab = tab
  )
  mod_10_doc_server(
    id = "doc",
    pre_process = pre_process,
    idep_data = idep_data,
    tab = tab
  )
}
