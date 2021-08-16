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


DATAPATH <- "/mnt/f/data/data103/"


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
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param species_id DESCRIPTION.
#' @param org_info DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
find_species_by_id <- function(species_id, org_info) {
  return(org_info[which(org_info$id == speciesID), ])
}


# just return name
#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param species_id DESCRIPTION.
#' @param org_info DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
find_species_by_id_name <- function(species_id, org_info) {
  # find species name use id
  return(org_info[which(org_info$id == species_id), 3])
}


convert_id <- function(query, species_choice, annotated_species_counts,
                       select_org = "BestMatch") {
  query <- gsub(pattern = "\"|\'", "", x = query)
  # remove " in gene ids, mess up SQL query
  # remove ' in gene ids
  # |\\.[0-9] remove anything after A35244.1 -> A35244
  #  some gene ids are like Glyma.01G002100

  query_set <- clean_gene_set(unlist(strsplit(toupper(query), "\t| |\n|\\,")))
  conn_db <- connect_convert_db()

  if (select_org == "BestMatch") { # query all species
    result <- DBI::dbGetQuery(
      conn = conn_db,
      statement = "select distinct id,ens,species,idType
     from mapping where id in (?)",
      params = query_set
    )
  } else { # organism has been selected query specific one
    result <- DBI::dbGetQuery(
      conn = conn_db,
      statement = "select distinct id,ens,species,idType
     from mapping where species = ? and id in (?)",
      params = list(select_org, query_set)
    )
  }
  DBI::dbDisconnect(conn = conn_db)

  if (nrow(result) == 0) {
    return(NULL)
  }

  if (select_org == species_choice[[1]]) { # if best match species
    combination <- paste(result$species, result$idType)
    sorted_counts <- sort(table(combination), decreasing = T)

    # Try to use Ensembl instead of STRING-db genome annotation
    # when  one species matched, it is a number, not a vector
    # if the #1 species and #2 are close
    # 1:3 are Ensembl species
    # and #2 come earlier (ensembl) than #1
    tmp <- sum(annotated_species_counts[1:3])
    if (!is.integer(sorted_counts) &&
      sorted_counts[1] <= sorted_counts[2] * 1.1 &&
      as.numeric(gsub(
        pattern = " .*", "",
        x = names(sorted_counts[1])
      )) > tmp &&
      as.numeric(gsub(
        pattern = " .*", "",
        x = names(sorted_counts[2])
      )) < tmp) {
      tem <- sorted_counts[2]
      sorted_counts[2] <- sorted_counts[1]
      names(sorted_counts)[2] <- names(sorted_counts)[1]
      sorted_counts[1] <- tem
      names(sorted_counts)[1] <- names(tem)
    }


    recognized <-
      result <- result[which(combination == names(sorted_counts[1])), ]


    species_matched <- sorted_counts
    tmp <- as.numeric(gsub(pattern = " .*", "", x = names(sorted_counts)))
    names(species_matched) <- sapply(X = tmp, FUN = find_species_by_id_name)
    species_matched <- as.data.frame(species_matched)

    if (length(sorted_counts) == 1) { # if only  one species matched
      species_matched[1, 1] <- paste(rownames(species_matched),
        "(", species_matched[1, 1], ")",
        sep = ""
      )
    } else { # if more than one species matched
      species_matched[, 1] <- as.character(species_matched[, 1])
      species_matched[, 1] <- paste(species_matched[, 1],
        " (", species_matched[, 2], ")",
        sep = ""
      )
      species_matched[1, 1] <- paste(
        species_matched[1, 1],
        "   ***Used in mapping***  To change, select from above and resubmit query."
      )
      species_matched <- as.data.frame(species_matched[, 1])
    }
  } else { # if species is selected
    result <- result[which(result$species == select_org), ]
    if (nrow(result) == 0) {
      return(NULL)
    } # stop("ID not recognized!")
    species_matched <- as.data.frame(paste(
      "Using selected species ",
      find_species_by_id_name(select_org)
    ))
  }

  # remove duplicates in query gene ids
  result <- result[which(!duplicated(result[, 1])), ]
  # remove duplicates in ensembl_gene_id
  result <- result[which(!duplicated(result[, 2])), ]
  colnames(species_matched) <- c("Matched Species (genes)")
  conversion_table <- result[, 1:2]
  colnames(conversion_table) <- c("User_input", "ensembl_gene_id")
  conversion_table$Species <- sapply(result[, 3], findSpeciesByIdName)
  return(list(
    originalIDs = query_set,
    ids = unique(result[, 2]),
    species = findSpeciesById(result$species[1]),
    species_matched = species_matched,
    conversion_table = conversion_table
  ))
}
