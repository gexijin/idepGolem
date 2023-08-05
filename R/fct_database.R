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



#' Connect to the convertIDs database and return the
#' objects.
#'
#' Create a database connection with the DBI package.
#'
#' @param datapath Folder path to the data file
#'
#' @export
#' @return Database connection.
connect_convert_db <- function(datapath = DATAPATH) {
  if (!file.exists(org_info_file)) {
    # download org_info and demo files to current folder
    withProgress(message = "Download demo data and species database", {
      incProgress(0.2)
      file_name <- paste0(db_ver, ".tar.gz")
      options(timeout = 3000)
      download.file(
        url = paste0(db_url, db_ver, "/", file_name),
        destfile = file_name,
        mode = "wb",
        quiet = FALSE
      )
      untar(file_name) # untar and unzip the files
      file.remove(file_name) # delete the tar file to save storage
    })
  }

  return(DBI::dbConnect(
    drv = RSQLite::dbDriver("SQLite"),
    dbname = org_info_file,
    flags = RSQLite::SQLITE_RO
  ))
}


#' Connect to the convertIDs database for the species and return the
#' objects.
#'
#' Create a database connection with the DBI package.
#'
#' @param datapath Folder path to the data file
#' @param select_org The slected species
#' @param idep_data  Data object that includes org_info
#'
#' @export
#' @return Database connection.
connect_convert_db_org <- function(datapath = DATAPATH, select_org, idep_data) {
  ix <- which(idep_data$org_info$id == select_org)
  db_file <- idep_data$org_info[ix, "file"]
  return(try(
    DBI::dbConnect(
    drv = RSQLite::dbDriver("SQLite"),
    dbname = paste0(datapath, "db/", db_file),
    flags = RSQLite::SQLITE_RO
  )))
}


#' Load iDEP data
#'
#' Use this function call to load data that is
#' used in other functions.
#'
#'
#' @param datapath Folder path to the iDEP data
#' @export
#' @return Large list of the iDEP data.
#' 1. kegg_species_id:  KEGG species list
#' 2. gmt_files: list of pathway files
#' 3. gene_info_files: list of geneInfo files
#' 4. demo_data_file: demo data file
#' 5. demo_metadata_file: experimental design file for demo data
#' 6. quotes:  quotes
#' 7. string_species_go_data: List of STRING species
#' 8. org_info: orgInfo for Ensembl species
#' 9. annotated_species_count: total number of annotated species
#' 10. go_levels: GO levels
#' 11. go_level_2_terms: mapping of GO levels to terms
#' 12. id_index: idtype and index
#' 13. species_choice: list of species for populating selection.
#'
get_idep_data <- function(datapath = DATAPATH) {

  conn_db <- connect_convert_db()

  # demo_data_info.csv file
  # columns: ID, expression, design, type, name
  # ID column must be unique
  # The type column 1- read counts, 2- normalized , 3 LFC&Pval
  # The files must be available in the demo folder of the database.
  # Leave blank if design file is missing
  # Default file have smallest ID, within the same type.

  demo_file_info <- read.csv(
    paste0(datapath, "demo/demo_data_info.csv")
  )
  # add path for expression matrix
  demo_file_info$expression <- paste0(
    datapath,
    "demo/",
    demo_file_info$expression
  )

  # add path for design file if exist
  ix <- which(nchar(demo_file_info$design) > 2)
  demo_file_info$design[ix] <- paste0(
    datapath,
    "demo/",
    demo_file_info$design[ix]
  )

  org_info <- DBI::dbGetQuery(
    conn = conn_db,
    statement = "select * from orgInfo;"
  )
  annotated_species_count <- sort(table(org_info$group))

  species_choice <- setNames(as.list(org_info$id), org_info$name2)
  species_choice <- append(
    setNames("NEW", "**NEW SPECIES**"),
    species_choice
  )

  #species_choice <- append(
  #  setNames("BestMatch", "Best matching species"),
  #  species_choice
  #)
  top_choices <- c(
    #"Best matching species", 
    #"**NEW SPECIES**", 
    "Human", "Mouse", "Rat", "Cow",
    "Zebrafish", "Pig", "Chicken", "Macaque", "Dog", "Drosophila melanogaster",
    "Caenorhabditis elegans", "Saccharomyces cerevisiae",
    "Arabidopsis thaliana", "Zea mays", "Glycine max",
    "Oryza sativa Indica Group", "Oryza sativa Japonica Group", "Vitis vinifera"
  )

  other_choices <- names(species_choice)[
    !(names(species_choice) %in% top_choices)
  ]
  species_choice <- species_choice[c(top_choices, other_choices)]
  #org_info <- org_info[order(org_info$group), ]
  ix <- match(org_info$name2, top_choices)
  org_info <- org_info[order(ix), ]
  org_info <- org_info[order(org_info$group == "STRINGv11.5"), ]
  # GO levels
  go_levels <- DBI::dbGetQuery(
    conn = conn_db,
    statement = "select distinct id, level from GO
         WHERE GO = 'biological_process'"
  )
  go_level_2_terms <- go_levels[which(go_levels$level %in% c(2, 3)), 1]

  quotes <- DBI::dbGetQuery(
    conn = conn_db,
    statement = "select quotes from quotes;"
  )
  quotes <- quotes[, 1]

  DBI::dbDisconnect(conn = conn_db)


  # id_index <- DBI::dbGetQuery(
  #  conn = conn_db,
  #  statement = "select distinct * from idIndex"
  # )


  return(list(
#    kegg_species_id = kegg_species_id,
#    gmt_files = gmt_files,
#    gene_info_files = gene_info_files,
    demo_file_info = demo_file_info,
    quotes = quotes,
#    string_species_go_data = string_species_go_data,
    org_info = org_info,
    annotated_species_count = annotated_species_count,
    go_levels = go_levels,
    go_level_2_terms = go_level_2_terms,
#    id_index = id_index,
    species_choice = species_choice
  ))
}


#' Find the ID type
#'
#' Using an inputted ID, find the type of IDs
#'
#' @param id ID to find the type from
#' @param id_index Index of IDs to use
#'
#' @export
find_id_type_by_id <- function(id, id_index) {
  # find idType based on index
  return(id_index$idType[as.numeric(id)])
}

#' Find a species by ID
#'
#' Find a species in the iDEP database with an
#' ID.
#'
#' @param species_id Species ID to search the database with
#' @param org_info iDEP data org_info file
#'
#' @export
#' @return Species information in \code{org_info} from the
#'  matched ID.
find_species_by_id <- function(species_id, org_info) {
  return(org_info[which(org_info$id == species_id), ])
}


#' Find a species by ID
#'
#' Find a species in the iDEP database with an
#' ID.
#'
#' @param species_id Species ID to search the database with
#' @param org_info iDEP data org_info file
#'
#' @export
#' @return Only return the species name with this function.
find_species_by_id_name <- function(species_id, org_info) {
  # find species name use id
  return(org_info[which(org_info$id == species_id), 3])
}

#' Find a species id by ensembl dataset name
#'
#' Find a species in the iDEP database with an
#' ID.
#'
#' @param species_id Species ID to search the database with
#' @param org_info iDEP data org_info file
#'
#' @export
#' @return Only return the species name with this function.
find_species_id_by_ensembl <- function(ensembl_dataset, org_info) {
  # find species name use id
  return(org_info[which(org_info$ensembl_dataset == ensembl_dataset), "id"])
}

#'  Convert gene IDs to ensembl data.
#'
#' @param query A character vector of gene IDs
#' @param idep_data A instance of the output from \code{\link{get_idep_data}()}
#' @param select_org A character of the species that wants to be looked up,
#'  default to \code{"BestMatch"}
#' @param max_sample_ids The number of gene ids used to determine species
#'
#' @export
#' @return A large list of the conversion ID information that was gathered
#'  from querying the database with the original IDs.
convert_id <- function(query,
                       idep_data,
                       select_org = "BestMatch",
                       max_sample_ids = 100) {
  # Solves the issue of app shut down when species
  # is deleted after genes are uploaded.
  if (is.null(select_org)) {
    return(NULL)
  }

  query <- gsub(pattern = "\"|\'", "", x = query)
  # remove " in gene ids, mess up SQL query
  # remove ' in gene ids
  # |\\.[0-9] remove anything after A35244.1 -> A35244
  #  some gene ids are like Glyma.01G002100
  query_set <- clean_query(query_input = query)
  query_string <- paste0("('", paste(query_set, collapse = "', '"), "')")

  # use a small set of genes to guess species and idType; to improve speed
  # use a small set of gene ids, with the max #
  # when the query is small, use the quary

  # acutal number of samples, for calculating % later
  n_sample_ids <- length(query_set)
  test_query_string <- query_string
  if (length(query_set) > max_sample_ids) {
    n_sample_ids <- max_sample_ids
    test_query_set <- sample(query_set, max_sample_ids)
    test_query_string <- paste0(
      "('",
      paste(test_query_set, collapse = "', '"),
      "')"
    )
  }

  conn_db <- connect_convert_db_org(
    select_org = select_org,
    idep_data = idep_data
  )

  # if database connection error
  # see ? try
  if(inherits(conn_db, "try-error")) {
    return(NULL)
  }


  # if species is selected ---------------------------------------------------
  query_statement <- paste0(
    "select id,ens,idType from mapping where id IN ",
    query_string
  )

  result <- DBI::dbGetQuery(conn_db, query_statement)

  DBI::dbDisconnect(conn_db)
  if (nrow(result) == 0) {
    return(NULL)
  }

  # resolve multiple ID types, get the most matched
  best_id_type <- as.integer(
    names(
      sort(
        table(result$idType),
        decreasing = TRUE
      )
    )[1]
  )
  result <- result[result$idType == best_id_type, ]

  matched <- as.data.frame(paste(
    "Selected:",
    find_species_by_id_name(
      species_id = select_org,
      org_info = idep_data$org_info
    )
  ))


  # Needs review---------------------------------------
  # one to many, keep one ensembl id, randomly
  # remove duplicates in query gene ids
  result <- result[which(!duplicated(result[, 1])), ]

  # many user id to one ensembl id mapping, keep all
  # remove duplicates in ensembl_gene_id
  # result <- result[which(!duplicated(result[, 2])), ]

  colnames(matched) <- c("Matched Species")

  conversion_table <- result[, 1:2]
  colnames(conversion_table) <- c("User_input", "ensembl_gene_id")
  conversion_table$Species <- find_species_by_id_name(
    species_id = as.integer(select_org),
    org_info = idep_data$org_info
  )

  species <- find_species_by_id(
    species_id = as.integer(select_org),
    org_info = idep_data$org_info
  )

  return(list(
    origninal_ids = query_set,
    ids = unique(result[, 2]),
    species = species,
    species_match = matched,
    conversion_table = conversion_table
  ))
}

#' Read pathway sets for gene query
#'
#' This function provides the gene set information
#' to perform enrichment analysis on the provided
#' query.
#'
#' @param all_gene_names_query Subsetted data frame of genes names from
#'  \code{\link{get_all_gene_names}()} to query database with.
#' @param converted List of conversion information on the original
#'   IDs returned from the \code{\link{convert_id}()}.
#' @param go String designating the section of the database to query for pathway
#'   analysis. See \code{\link{gmt_category}()} for choices.
#' @param select_org String designating what organism is being analyzed.
#' @param gmt_file For NEW species the gmt file to use
#'   for the pathway analysis
#' @param idep_data Data built in to idep
#' @param gene_info The gene info from the converted IDs and
#'   the function gene_info()
#' 
#' @export
#' @return This function returns a list with values that are
#'   used in the find_overlap function. The list contains
#'   pathway_table which is the overlap and total genes for
#'   each pathway that is enriched in the query. The list
#'   also contains the query_set of genes, the total_genes
#'   number which is used in the calculation of the p-values
#'   and the pathway files that contain gmt information on
#'   the mathced species.
read_pathway_sets <- function(all_gene_names_query,
                              converted,
                              go,
                              select_org,
                              gmt_file,
                              idep_data,
                              gene_info) {
  id_not_recognized <- as.data.frame("ID not recognized!")

  if (select_org == "NEW" && is.null(gmt_file)) {
    return(as.data.frame("No GMT file provided!"))
  } else if (select_org == "NEW" && !is.null(gmt_file)) {
    in_file <- gmt_file
    in_file <- in_file$datapath

    return(read_gmt(in_file))
  }

  if (ncol(all_gene_names_query) == 1) {
    return(id_not_recognized)
  }

  query_set <- all_gene_names_query[, 2]

  #  if (!is.null(gene_info)) {
  #    if (dim(gene_info)[1] > 1) {
  #      gene_info <- gene_info[which(
  #        gene_info$gene_biotype == "protein_coding"
  #      ), ]
  #      query_set <- intersect(query_set, gene_info[, 1])
  #    }
  #  }

  if (length(query_set) == 0) {
    return(id_not_recognized)
  }

  pathway <- connect_convert_db_org(
    select_org = select_org,
    idep_data = idep_data
  )

  # if database connection error
  if(inherits(pathway, "try-error")) {
    return(id_not_recognized)
  }

  # does the db file has a categories table?
  ix <- grep(
    pattern = "pathway",
    x = DBI::dbListTables(pathway)
  )

  if (length(ix) == 0) {
    return(id_not_recognized)
  }

  if (is.null(go)) {
    go <- "GOBP"
  }

  sql_query <- build_pathway_query(go, query_set)

  result <- DBI::dbGetQuery(pathway, sql_query)

  if (dim(result)[1] == 0) {
    return(pathway_table <- NULL)
  }

  # List pathways and frequency of genes
  pathway_ids <- stats::aggregate(
    result$pathwayID,
    by = list(unique_values = result$pathwayID),
    FUN = length
  )
  colnames(pathway_ids) <- c("pathway_id", "overlap")

  if (dim(pathway_ids)[1] == 0) {
    return(pathway_table <- NULL)
  } else {

    pathway_info <- DBI::dbGetQuery(
      pathway,
      paste(
        "SELECT DISTINCT id,n,Description,memo FROM pathwayInfo WHERE id IN ('",
        paste(pathway_ids[, 1], collapse = "', '"), "') ",
        sep = ""
      )
    )

    # if duplicates pathway_info, remove
    # sorted so that if some GO categories are repeated, keep the one with URL
    pathway_info <- pathway_info[
      order(
        pathway_info$id,
        pathway_info$memo,
        decreasing = TRUE
      ),
    ]
    pathway_info <- pathway_info[!duplicated(pathway_info$id), ]

    pathway_merge <- merge(
      x = pathway_ids,
      y = pathway_info,
      by.x = "pathway_id",
      by.y = "id",
      all = TRUE
    )

    # Convert pathways into lists
    gene_sets <- lapply(
      pathway_merge$pathway_id,
      function(x) result[which(result$pathwayID == x), 1]
    )
    names(gene_sets) <- pathway_merge$description
    # note this might cause errors if length differs.
    # which could happen if 
    pathway_merge$gene_sets <- gene_sets
  }

  # list of query genes that have at least one pathway in db
  query_set_db <- unique(result$gene)
  query_set <- query_set_db # only use genes included in DB

  # total number of genes in pathway db
  sql_query <- "SELECT COUNT ( DISTINCT gene ) FROM pathway"

  if (go != "All") {
    sql_query <- paste(
      sql_query,
      " WHERE category='", go, "'",
      sep = ""
    )
  }
  total_genes_db <- DBI::dbGetQuery(
    pathway,
    sql_query
  )
  total_genes_db <- as.integer(total_genes_db)

  total_genes <- total_genes_db

  DBI::dbDisconnect(pathway)

  # Gene sets and info for the enrichment analysis
  return(list(
    pathway_table = pathway_merge,
    query_set = query_set,
    total_genes = total_genes
#    pathway_files = pathway_files
  ))
}

#' Background gene pathway sets
#'
#' This function reads the pathway sets for the filtered
#' gene IDs and performs pathway analysis for the query
#' with the filtered background.
#'
#' @param processed_data Matrix of gene data that has been through
#'   \code{\link{pre_process}()}
#' @param gene_info Dataframe of converted IDs information from the function
#'   \code{\link{gene_info}()}
#' @param sub_query Vector of IDs that the enrichment
#'   analysis should be performed on. This list is also returned from
#'   \code{\link{read_pathway_sets}()}.
#' @param go String designating the section of the database to query for pathway
#'   analysis. See \code{\link{gmt_category}()} for choices.
#' @param pathway_table Data frame of results from
#'  \code{\link{read_pathway_sets}()}. If this data frame is NULL or 0 rows
#'  there this function will return no significant enrichment.
#' @param idep_data List of data returned from \code{\link{get_idep_data}()}
#' @param sub_pathway_files String designating file location for GMT files in
#'   the database that contain information for the matched species. This string
#'   is returned from \code{\link{read_pathway_sets}()}.
#' @param select_org Species selected.
#'
#' @export
#' @return Pathway gene set table for the background genes. Used
#'   in find_overlap to calculate pvals for the filtered background
#'
#' @seealso This function is used internally in \code{\link{find_overlap}()}
#'
background_pathway_sets <- function(processed_data,
                                    gene_info,
                                    sub_query,
                                    go,
                                    pathway_table,
                                    idep_data,
                                    select_org,
                                    sub_pathway_files) {
  query_set <- rownames(processed_data)

  #  if (!is.null(gene_info)) {
  #    if (dim(gene_info)[1] > 1) {
  #      query_set <- intersect(
  #        query_set,
  #        gene_info[which(gene_info$gene_biotype == "protein_coding"), 1]
  #      )
  #    }
  #  }

  pathway <- connect_convert_db_org(
    select_org = select_org,
    idep_data = idep_data
  )
  # if database connection error
  if(inherits(pathway, "try-error")) {
    return(NULL)
  }

  if (length(intersect(query_set, sub_query)) == 0) {
    # If none of the selected genes are in background genes
    return(list(
      bg_result = as.data.frame(
        "None of the selected genes are in the background genes!"
      )
    ))
  }

  # Make sure the background set includes the query set
  query_set <- unique(c(query_set, sub_query))

  if (is.null(go)) {
    go <- "GOBP"
  }

  # this is slow when gene lists are big
  #sql_query <- build_pathway_query(go, query_set)
  #results <- DBI::dbGetQuery(pathway, sql_query)

  # retrieve all pathways for the category
  sql_query <- "SELECT gene, pathwayID FROM pathway "
  if (go != "All") {
    sql_query <- paste0(sql_query, " WHERE category = '", go, "'")
  }

  # since there are so many genes, this takes a long time
  # we are not using the genes, just query all the pathways for the category
  #sql_query <- build_pathway_query(go, query_set)

  results <- DBI::dbGetQuery(pathway, sql_query)

  # only keep genes in the query_set, this is faster than query using SQL
  results <- results[results$gene %in% query_set, ]

  if (dim(results)[1] == 0) {
    return(list(
      bg_result = as.data.frame("No matching species or gene ID file!")
    ))
  }
  bg_result <- table(results$pathwayID)
  bg_result <- as.data.frame(bg_result)
  colnames(bg_result) <- c("pathway_id", "overlap_bg")

  pathway_table_bg <- merge(
    pathway_table,
    bg_result,
    by = "pathway_id",
    all.x = TRUE
  )

  # list of background genes that have at least one pathway in db
  pathway_table_bg$total_genes_bg <- length(unique(results$gene))

  DBI::dbDisconnect(pathway)

  return(pathway_table_bg)
}

#' Database choices for the converted IDs
#'
#' This function provides a list of the portions of the
#' database to query that contain genes from the converted
#' IDs.
#'
#' @param converted List of converted gene information  from
#'   \code{\link{convert_id}()}
#' @param converted_data Data matrix element with the converted ensembl_IDs from
#'   the list returned from \code{\link{convert_data}()}
#' @param select_org String designating the organism being analyzed
#' @param gmt_file For NEW species the gmt file to use
#'   for the pathway analysis, otherwise \code{NULL}
#' @param idep_data List of data returned from \code{\link{get_idep_data}()}
#'
#' @export
#' @return A character vector of names of the section of data to perform
#'   pathway analysis on
gmt_category <- function(converted,
                         converted_data,
                         select_org,
                         gmt_file,
                         idep_data) {
  if (select_org == "NEW" && !is.null(gmt_file)) {
    return(list(custom_gene_set = "Custom"))
  }

  id_not_recognized <- as.data.frame("ID not recognized!")

  if (is.null(converted)) {
    return(id_not_recognized)
  }

  query_set <- rownames(converted_data)

  if (length(query_set) == 0) {
    return(id_not_recognized)
  }

  conn_db <- connect_convert_db_org(
    select_org = select_org,
    idep_data = idep_data
  )

  # if database connection error
  if(inherits(conn_db, "try-error")) {
    return(id_not_recognized)
  }

  # does the db file has a categories table?
  ix <- grep(
    pattern = "categories",
    x = DBI::dbListTables(conn_db)
  )

  if (length(ix) == 0) {
    return(id_not_recognized)
  }


  # Generate a list of geneset categories such as "GOBP", "KEGG" from file
  gene_set_category <- DBI::dbGetQuery(conn_db, "select distinct * from categories")
  DBI::dbDisconnect(conn_db)

  gene_set_category <- sort(gene_set_category[, 1])
  category_choices <- setNames(as.list(gene_set_category), gene_set_category)

  # Set order of popular elements
  top_choices <- c("GOBP", "GOCC", "GOMF", "KEGG")

  # if KEGG or others are not in the list, remove
  top_choices <- top_choices[top_choices %in% gene_set_category]

  other_choices <- names(category_choices)[
    !(names(category_choices) %in% top_choices)
  ]
  category_choices <- category_choices[c(top_choices, other_choices)]

  # Change names to the full description for display
  names(category_choices)[match("GOBP", category_choices)] <- "GO Biological Process"
  names(category_choices)[match("GOCC", category_choices)] <- "GO Cellular Component"
  names(category_choices)[match("GOMF", category_choices)] <- "GO Molecular Function"
  category_choices <- append(setNames("All", "All available gene sets"), category_choices)

  return(category_choices)
}

#' Read pathway gene sets
#'
#' Use the IDs from the converted database to find the
#' gene sets for all the pathways in the database. Returns
#' a list with each entry a vector of IDs corresponding to
#' a description of a pathway in the database.
#'
#' @param converted List of converted gene information  from
#'   \code{\link{convert_id}()}
#' @param all_gene_names Data frame of gene names from
#'   \code{\link{get_all_gene_names}()}
#' @param go String designating the section of the database to query for pathway
#'   analysis. See \code{\link{gmt_category}()} for choices.
#' @param select_org String designating with organism is being analyzed
#' @param idep_data List of data returned from \code{\link{get_idep_data}()}
#' @param my_range Vector of the (min_set_size, max_set_size)
#' 
#' @export
#' @return A list with each entry a list of gene IDs that correspond to
#'  a pathway.
read_gene_sets <- function(converted,
                           all_gene_names,
                           go,
                           select_org,
                           idep_data,
                           my_range) {
  id_not_recognized <- as.data.frame("ID not recognized!")
  if (is.null(converted)) {
    return(id_not_recognized)
  }
  if (!is.null(all_gene_names[, 2])) {
    query_set <- all_gene_names[, 2]
  }
  if (is.null(query_set) || length(query_set) == 0) {
    return(id_not_recognized)
  }

  pathway <- connect_convert_db_org(
    select_org = select_org,
    idep_data = idep_data
  )

# for testing
#pathway <- DBI::dbConnect(
#    drv = RSQLite::dbDriver("SQLite"),
#    dbname = "C:/work/iDEP_data/data104b/pathwayDB/Human__hsapiens_gene_ensembl.db",
#    flags = RSQLite::SQLITE_RO
#  )
# browser()
  # if database connection error
  if(inherits(pathway, "try-error")) {
    return(id_not_recognized)
  }

  # does the db file has a categories table?
  ix <- grep(
    pattern = "pathway",
    x = DBI::dbListTables(pathway)
  )

  if (length(ix) == 0) {
    return(id_not_recognized)
  }
  if (is.null(go)) {
    go <- "GOBP"
  }

  # retrieve all pathways for the category
  sql_query <- "SELECT  gene, pathwayID FROM pathway "
  if (go != "All") {
    sql_query <- paste0(sql_query, " WHERE category = '", go, "'")
  }

  # since there are so many genes, this takes a long time
  # we are not using the genes, just query all the pathways for the category
  #sql_query <- build_pathway_query(go, query_set)

  result <- DBI::dbGetQuery(pathway, sql_query)

  # only keep genes in the query_set, this is faster than query using SQL
  result <- result[result$gene %in% query_set, ]

  if (dim(result)[1] == 0) {
    return(list(x = as.data.frame("No matching species or gene ID file!")))
  }
  # List pathways and frequency of genes
  pathway_ids <- aggregate(
    result$pathwayID,
    by = list(unique.values = result$pathwayID),
    FUN = length
  )
  pathway_ids <- pathway_ids[which(pathway_ids[, 2] >= my_range[1]), ]
  pathway_ids <- pathway_ids[which(pathway_ids[, 2] <= my_range[2]), ]
  if (dim(pathway_ids)[1] == 0) {
    gene_sets <- NULL
  }

  # filter pathways
  result <- result[result$pathwayID %in% pathway_ids[, 1], ]
  # convert into lists using the split funciton.
  gene_sets <- split(result$gene, result$pathwayID)

  pathway_info <- DBI::dbGetQuery(
    pathway,
    paste(
      "select distinct id,Description,memo from pathwayInfo where id IN ('",
      paste(pathway_ids[, 1], collapse = "', '"),
      "') ",
      sep = ""
    )
  )
  # add pathway name to gene sets
  ix <- match(names(gene_sets), pathway_info[, 1])
  names(gene_sets) <- pathway_info[ix, 2]
  DBI::dbDisconnect(pathway)
  return(
    list(
      gene_lists = gene_sets,
      pathway_info = pathway_info
    )
  )
}

#' Convert IDs from ensembl to entrez
#'
#' Convert an ID qeury for a species from the ensembl
#' ID type to entrez type.
#'
#' @param query Vector of IDs to convert
#' @param species String designating the organism being analyzed
#' @param org_info org_info file from the list returned from
#'   \code{link{get_idep_data}()}
#' @param idep_data  Data object with species info for connecting to database
#'
#'
#' @export
#' @return The queried genes with converted IDs.
convert_ensembl_to_entrez <- function(query,
                                      species,
                                      org_info,
                                      idep_data) {
  query_set <- clean_gene_set(
    unlist(strsplit(toupper(names(query)), "\t| |\n|\\, "))
  )

  # Note uses species Identifying
  species_id <- org_info$id[which(org_info$ensembl_dataset == species)]
  # idType 6 for entrez gene ID
  convert <- connect_convert_db_org(
    select_org = species_id,
    idep_data = idep_data
  )
  # if database connection error
  if(inherits(convert, "try-error")) {
    return(NULL)
  }
  id_type_entrez <- DBI::dbGetQuery(
    convert,
    paste(
      "select distinct * from idIndex where idType = 'entrezgene_id'"
    )
  )
  if (dim(id_type_entrez)[1] != 1) {
    cat("Warning! entrezgene ID not found!")
  }
  id_type_entrez <- as.numeric(id_type_entrez[1, 1])

  result <- DBI::dbGetQuery(
    convert,
    paste0(
      "select  id,ens from mapping where ",
      " idType ='", id_type_entrez, "' ",
      " AND ens IN ('", paste(query_set, collapse = "', '"), "')"
    )
  )

  DBI::dbDisconnect(convert)
  if (dim(result)[1] == 0) {
    return(NULL)
  }
  ix <- match(result$ens, names(query))

  tem <- query[ix]
  names(tem) <- result$id
  return(tem)
}


#' Create SQL query statement for reading pathway db
#'
#' Convert an ID qeury for a species from the ensembl
#' ID type to entrez type.
#'
#' @param go  Pathway category "KEGG", "GOBP"
#' @param query_set  List of genes
#' @export
#' @return The SQL SELECT statement
build_pathway_query <- function(go, query_set) {

  sql_query <- "SELECT gene, pathwayID FROM pathway WHERE "

  # faster if category is first
  if (go != "All") {
    sql_query <- paste0(sql_query, " category = '", go, "' AND ")
  }

  # Get Gene sets
  sql_query <- paste(
    sql_query,
    " gene IN ('",
    paste(query_set, collapse = "', '"),
    "')",
    sep = ""
  )
  return(sql_query)
}
