#' fct_database.R: A fie with all functions to deal with database
#'
#'
#' @section fct_database.R functions:
#' \code{connect_convert_db} connects to the converID.db 
#' and returns the objects.
#'
#'
#'
#' @name fct_database.R
NULL


DATAPATH <- "../data/data103/"


#' connect_convert_db connects to the converID.db and returns the objects.
#'
#'
#' @description
#'
#'
#' @param datapath
#'
#'
#' @return
#'
#'
#' @examples
connect_convert_db <- function(datapath = DATAPATH) {
  return(DBI::dbConnect(
    drv = RSQLite::dbDriver("SQLite"),
    dbname = paste0(datapath, "convertIDs.db"),
    flags = RSQLite::SQLITE_RO
  ))
}


#' get_idep_data
#'
#'
#' @description
#'
#'
#' @param datapath
#'
#'
#' @return
#'
#'
#' @examples
get_idep_data <- function(datapath = DATAPATH) {
  kegg_species_id <- read.csv(paste0(datapath, "data_go/KEGG_Species_ID.csv"))

  gmt_files <- list.files(
    path = paste0(datapath, "pathwayDB"),
    pattern = ".*\\.db"
  )
  gmt_files <- paste(datapath,
    "pathwayDB/", gmt_files,
    sep = ""
  )

  gene_info_files <- list.files(
    path = paste0(datapath, "geneInfo"),
    pattern = ".*GeneInfo\\.csv"
  )
  gene_info_files <- paste(datapath,
    "geneInfo/", gene_info_files,
    sep = ""
  )

  demo_data_file <- paste0(datapath, "data_go/BcellGSE71176_p53.csv")
  demo_metadata_file <- paste0(
    datapath,
    "data_go/BcellGSE71176_p53_sampleInfo.csv"
  )

  conn_db <- connect_convert_db()

  quotes <- DBI::dbGetQuery(
    conn = conn_db,
    statement = "select * from quotes"
  )
  quotes <- paste0(
    "\"", quotes$quotes,
    "\"", " -- ", quotes$author, ".       "
  )

  string_species_go_data <- read.csv(paste0(
    datapath,
    "data_go/STRING11_species.csv"
  ))

  org_info <- DBI::dbGetQuery(
    con = conn_db,
    statement = "select distinct * from orgInfo"
  )

  annotated_species_count <- sort(table(org_info$group))

  go_levels <- DBI::dbGetQuery(
    conn = conn_db,
    statement = "select distinct id, level from GO
         WHERE GO = 'biological_process'"
  )

  go_level_2_terms <- go_levels[which(go_levels$level %in% c(2, 3)), 1]

  id_index <- DBI::dbGetQuery(
    conn = conn_db,
    statement = "select distinct * from idIndex"
  )

  DBI::dbDisconnect(conn = conn_db)

  return(list(
    kegg_species_id = kegg_species_id,
    gmt_files = gmt_files,
    gene_info_files = gene_info_files,
    demo_data_file = demo_data_file,
    demo_metadata_file = demo_metadata_file,
    quotes = quotes,
    string_species_go_data = string_species_go_data,
    org_info = org_info,
    annotated_species_count = annotated_species_count,
    go_levels = go_levels,
    go_level_2_terms = go_level_2_terms,
    id_index = id_index
  ))
}

# find idType based on index
find_id_type_by_id <- function(id, id_index) { # find
  return(id_index$idType[as.numeric(id)])
}


# find species name use id
find_species_by_id <- function(species_id, org_info) {
  return(org_info[which(org_info$id == speciesID), ])
}
# just return name
find_species_by_id_name <- function(species_id, org_info) { # find species name use id
  return(org_info[which(org_info$id == species_id), 3])
}
