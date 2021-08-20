#' fct_01_load_data.R This file holds all of the main data analysis functions
#' associated with first tab of the iDEP website.
#'
#'
#' @section fct_01_load_data.R functions:
#'
#'
#'
#' @name fct_01_load_data.R
NULL

# retrieve detailed info on genes
gene_info <- function(converted, select_org, idep_data) {
  check <- check_object_state(
    check_exp = (is.null(converted)),
    true_message = as.data.frame("ID not recognized!")
  )
  if (check$bool) {
    return(check)
  }

  query_set <- converted$ids
  check <- check_object_state(
    check_exp = (length(query_set) == 0),
    true_message = as.data.frame("ID not recognized!")
  )
  if (check$bool) {
    return(check)
  }

  ix <- grep(
    pattern = converted$species[1, 1],
    x = idep_data$gene_info_files
  )
  check <- check_object_state(
    check_exp = (length(ix) == 0),
    true_message = as.data.frame("No matching gene info file found")
  )
  if (check$bool) {
    return(check)
  }
  
  # If selected species is not the default "bestMatch",
  # use that species directly
  if (select_org != idep_date$species_choice[[1]]) {
    ix <- grep(
      pattern = find_species_by_id(
        species_id = select_org,
        org_info = idep_data$org_info
      )[1, 1],
      x = idep_data$gene_info_files
    )
  }

  check <- check_object_state(
    check_exp = (length(ix) != 1),
    true_message = as.data.frame("Multiple geneInfo file found!")
  )
  if (check$bool) {
    return(check)
  } else {
    gene_info_csv <- read.csv(as.character(idep_data$geneInfoFiles[ix]))
    gene_info_csv[, 1] <- toupper(gene_info_csv[, 1])
  }

  set <- match(gene_info_csv$ensembl_gene_id, query_set)
  set[which(is.na(set))] <- "Genome"
  set[which(set != "Genome")] <- "List"
  return(cbind(gene_info, set))
}