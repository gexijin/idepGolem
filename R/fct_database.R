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


DATAPATH <- "D:/data/data103/"


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
  org_info <- org_info[order(org_info$name), ]

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

  species_choice <- setNames(as.list(org_info$id), org_info$name2)
  species_choice <- append(
    setNames("NEW", "**NEW SPECIES**"),
    species_choice
  )
  species_choice <- append(
    setNames("BestMatch", "Best matching species"),
    species_choice
  )
  top_choices <- c(
    "Best matching species", "**NEW SPECIES**", "Human", "Mouse", "Rat", "Cow",
    "Zebrafish", "Pig", "Chicken", "Macaque", "Dog", "Drosophila melanogaster",
    "Caenorhabditis elegans", "Saccharomyces cerevisiae",
    "Arabidopsis thaliana", "Zea mays", "Glycine max",
    "Oryza sativa Indica Group", "Oryza sativa Japonica Group", "Vitis vinifera"
  )
  other_choices <- names(species_choice)[
    !(names(species_choice) %in% top_choices)
  ]
  species_choice <- species_choice[c(top_choices, other_choices)]
  
  genome_assembl <- paste0(datapath, "data_go/orgInfo_animal_plant_metazoa.csv")
  genome_assembl <- read.csv(genome_assembl)
  genome_assembl <- genome_assembl |>
    dplyr::select(name2, academicName, assembly) |>
    dplyr::rename(
      "Species Name" = name2, "Academic Name" = academicName,
      "Genome Assembly" = assembly
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
    id_index = id_index,
    species_choice = species_choice,
    genome_assembl = genome_assembl
  ))
}

idep_data <- get_idep_data()


#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param id DESCRIPTION.
#' @param id_index DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
find_id_type_by_id <- function(id, id_index) {
  # find idType based on index
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
  return(org_info[which(org_info$id == species_id), ])
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



#' convert_id This function takes gene IDs and converts them to ensembl data.
#'
#' FUNCTION_DESCRIPTION
#'
#' @param query A character vector of gene IDs
#' @param idep_data A instance of the output from get_idep_data
#'  (link to documentation)
#' @param select_org A character of the species that wants to be looked up,
#'  default to \code{"BestMatch"}
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
convert_id <- function(query, idep_data,
                       select_org = "BestMatch") {
  query <- gsub(pattern = "\"|\'", "", x = query)
  # remove " in gene ids, mess up SQL query
  # remove ' in gene ids
  # |\\.[0-9] remove anything after A35244.1 -> A35244
  #  some gene ids are like Glyma.01G002100
  query_set <- clean_query(query_input = query)

  conn_db <- connect_convert_db()

  if (select_org == "BestMatch") { # query all species
    result <- DBI::dbGetQuery(
      conn = conn_db,
      statement = "select distinct id,ens,species,idType
     from mapping where id in (?)",
      params = list(query_set)
    )
  } else { # organism has been selected query specific one
    result <- DBI::dbGetQuery(
      conn = conn_db,
      statement = "select distinct id,ens,species,idType
     from mapping where species = ? and id in (?)",
      params = list(rep(select_org, length(query_set)), query_set)
    )
  }
  DBI::dbDisconnect(conn = conn_db)

  if (nrow(result) == 0) {
    return(NULL)
  }

  if (select_org == idep_data$species_choice[[1]]) { # if best match species
    combination <- paste(result$species, result$idType)
    sorted_counts <- sort(table(combination), decreasing = T)

    # Try to use Ensembl instead of STRING-db genome annotation
    # when  one species matched, it is a number, not a vector
    # if the #1 species and #2 are close
    # 1:3 are Ensembl species
    # and #2 come earlier (ensembl) than #1
    tmp <- sum(idep_data$annotated_species_counts[1:3])
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

    result <- result[which(combination == names(sorted_counts[1])), ]
    species_matched <- sorted_counts
    tmp <- as.numeric(gsub(pattern = " .*", "", x = names(sorted_counts)))
    names(species_matched) <- sapply(
      X = tmp,
      FUN = find_species_by_id_name,
      org_info = idep_data$org_info
    )
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
        "***Used in mapping***"
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
      find_species_by_id_name(
        species_id = select_org,
        org_info = idep_data$org_info
      )
    ))
  }

  # remove duplicates in query gene ids
  result <- result[which(!duplicated(result[, 1])), ]
  # remove duplicates in ensembl_gene_id
  result <- result[which(!duplicated(result[, 2])), ]
  colnames(species_matched) <- c("Matched Species (genes)")
  conversion_table <- result[, 1:2]
  colnames(conversion_table) <- c("User_input", "ensembl_gene_id")
  conversion_table$Species <- sapply(
    X = result[, 3],
    FUN = find_species_by_id_name,
    org_info = idep_data$org_info
  )

  species <- find_species_by_id(
    species_id = result$species[1],
    org_info = idep_data$org_info
  )

  return(list(
    origninal_ids = query_set,
    ids = unique(result[, 2]),
    species = species,
    species_matched = species_matched,
    conversion_table = conversion_table
  ))
}



#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param query DESCRIPTION.
#' @param species DESCRIPTION.
#' @param idep_date DESCRIPTION.
#' @param convert_type DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
convert_ensembl <- function(query, species, idep_date,
                            convert_type = "entrez") {
  query_set <- clean_query(query_input = query)
  # note uses species Identifying
  species_id <-
    idep_date$orgInfo$id[which(idep_date$orgInfo$ensembl_dataset == species)]

  conn_db <- connect_convert_db()
  if (convert_type == "entrez") {
    id_type_entrez <- DBI::dbGetQuery(
      conn = conn_db,
      statement = "select distinct * from idIndex
      where idType = 'entrezgene_id'"
    )
    id_type_entrez <- as.numeric(id_type_entrez[1, 1])
    result <- DBI::dbGetQuery(
      conn = conn_db,
      statement = "select id,ens,species from mapping
     where ens IN (?) AND idType = ?",
      params = list(query_set, id_type_entrez)
    ) # slow
    # idType 6 for entrez gene ID
  } else {
    id_type_kegg <- DBI::dbGetQuery(
      conn = conn_db,
      statement = "select distinct * from idIndex
      where idType = 'kegg'"
    )
    id_type_entrez <- as.numeric(id_type_entrez[1, 1])
    result <- DBI::dbGetQuery(
      conn = conn_db,
      statement = "select id,ens,species from mapping
     where ens IN (?) AND idType = ?",
      params = list(query_set, id_type_kegg)
    ) # slow
  }
  DBI::dbDisconnect(conn = conn_db)

  if (nrow(result) == 0) {
    return(NULL)
  }
  result <- subset(
    x = result,
    subset = species == species_id,
    select = -species
  )

  ix <- match(x = result$ens, table = names(query))

  tem <- query[ix]
  names(tem) <- result$id
  return(tem)
}

#' READ GENE SETS
read_pathway_sets <- function (
  all_gene_names,
  converted,
  go,
  select_org,
  gmt_range,
  gmt_file,
  idep_data,
  gene_info
) {
  browser()
	id_not_recognized = as.data.frame("ID not recognized!")

  if(select_org == "NEW" && !is.null(gmt_file)) {
    in_file <- gmt_file
    in_file <- in_file$datapath

    return(read_gmt(in_file))
  }

	if(ncol(all_gene_names) == 1) {
    return(id_not_recognized)
  }

  query_set <- all_gene_names[, 2]

  if(!is.null(gene_info)) {
    if(dim(gene_info)[1] > 1) {  
	    gene_info <- gene_info[which(
        gene_info$gene_biotype == "protein_coding"
      ), ]  
	    query_set <- intersect(query_set, gene_info[, 1])
	  }
  }

  if(length(query_set) == 0) {
    return(id_not_recognized)
  } 

	ix = grep(converted$species[1,1], idep_data$gmt_files)
  total_genes <- converted$species[1,7]

	# If selected species is not the default "bestMatch", use that species directly
	if(select_org != "BestMatch") {  
		ix = grep(find_species_by_id(select_org)[1, 1], idep_data$gmt_files)
		total_genes <- org_info[which(
      org_info$id == as.numeric(select_org)
    ), 7]
	}

  if (length(ix) == 0) {
    return(id_not_recognized )
  }

	pathway <- DBI::dbConnect(
    drv = RSQLite::dbDriver("SQLite"),
    dbname = idep_data$gmt_files[ix],
    flags= RSQLite::SQLITE_RO
  )
	
	if(is.null(go)) {
    go <- "GOBP"
  }

	sql_query = paste(
    "select distinct gene, pathwayID from pathway where gene IN ('",
    paste(query_set, collapse = "', '"), "')",
    sep = ""
  )
	# cat(paste0("\n\nhere:",GO,"There"))

	if(go != "All") {
    sql_query = paste0(sql_query, " AND category ='", go,"'")
  } 
	result <- DBI::dbGetQuery(pathway, sql_query)

	if(dim(result)[1] == 0) {
    return(list(x = as.data.frame("No matching species or gene ID file!")))
  }

	# list pathways and frequency of genes
	pathway_ids <- stats::aggregate(
    result$pathwayID,
    by = list(unique.values = result$pathwayID),
    FUN = length
  )
	pathway_ids <- pathway_ids[which(pathway_ids[, 2] >= gmt_range[1]), ]
	pathway_ids <- pathway_ids[which(pathway_ids[, 2] <= gmt_range[2]), ]
  colnames(pathway_ids) <- c("pathway_id", "overlap")

	if(dim(pathway_ids)[1] == 0) {
    gene_sets = NULL
  } else {
    # Convert pathways into lists
	  gene_sets <- lapply(
      pathway_ids[, 1],
      function(x)  result[which(result$pathwayID == x), 1]
    )
	  pathway_info <- DBI::dbGetQuery(
      pathway,
      paste(
        "select distinct id, n, Description from pathwayInfo where id IN ('", 
				paste(pathway_ids[, 1], collapse = "', '"), "') ", sep = ""
      ) 
    )
	  ix <- match(pathway_ids[, 1], pathway_info[, 1])
    pathway_merge <- merge(
      x = pathway_ids,
      y = pathway_info,
      by.x = "pathway_id",
      by.y = "id"
    )
	  names(gene_sets) <- pathway_info[ix, 2]  

    test <- total_genes - length(query_set)
	  if (test < 0) {
	    test <- 0
 	  }
  }
	
	DBI::dbDisconnect(pathway)

  # Gene sets and info for the enrichment analysis
	return(list(
    gene_sets = gene_sets,
    pathway_info = pathway_merge,
    total_gene_test = test
  ))
} 

#' Database choices for the converted IDs
gmt_category <- function(
  converted,
  converted_data,
  select_org,
  gmt_file,
  idep_data
) {
	if (select_org == "NEW" && !is.null(gmt_file)) {
    return(list(custom_gene_set ="Custom"))
  }
		
	id_not_recognized = as.data.frame("ID not recognized!")

	if(is.null(converted)) {
    return(id_not_recognized)
  }

	query_set <- rownames(converted_data)

	if(length(query_set) == 0){
    return(id_not_recognized )
  }

	ix = grep(converted$species[1,1], idep_data$gmt_files)

	if (length(ix) == 0) {
    return(id_not_recognized)
  }
	
	# If selected species is not the default "bestMatch", use that species directly
	if(select_org != idep_data$species_choice[[1]]) {  
		ix = grep(find_species_by_id(select_org)[1,1], idep_data$gmt_files)
		if (length(ix) == 0) {
      return(id_not_recognized)
    }
	}

	pathway <- DBI::dbConnect(
    drv = RSQLite::dbDriver("SQLite"),
    dbname = idep_data$gmt_files[ix],
    flags= RSQLite::SQLITE_RO
  )

	# Generate a list of geneset categories such as "GOBP", "KEGG" from file
	gene_set_category <-  DBI::dbGetQuery(pathway, "select distinct * from categories") 
	gene_set_category  <- sort(gene_set_category[, 1])
	category_choices <- setNames(as.list(gene_set_category ), gene_set_category )
	
	# Set order of popular elements
  top_choices <- c("GOBP", "GOCC", "GOMF", "KEGG")
  other_choices <- names(category_choices)[
    !(names(category_choices) %in% top_choices)
  ]
  category_choices <- category_choices[c(top_choices, other_choices)]

	# Change names to the full description for display
	names(category_choices)[match("GOBP", category_choices)] <- "GO Biological Process"
	names(category_choices)[match("GOCC", category_choices)] <- "GO Cellular Component"
	names(category_choices)[match("GOMF", category_choices)] <- "GO Molecular Function"
	category_choices <- append(setNames("All", "All available gene sets"), category_choices)
	
	DBI::dbDisconnect(pathway)

	return(category_choices )
} 