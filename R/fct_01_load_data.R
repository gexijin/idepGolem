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




#' Load basic data information
#' 
#' This function does the immediate loading of the data and
#' sample info to present in the data preview table and the
#' sample info table. The data undergoes very basic filtering
#' and transformation before entering the table.
#' 
#' @param expression_file The data path for the input into the
#' expression file file input bar
#' @param experiment_file The data path for the input into the
#' experiment file file input bar
#' @param go_button Action button input that tells the app to
#' load the demo data files
#' @param demo_data_file Expression demo data path
#' @param demo_metadata_file Experiment demo data path
#' 
#' @return This returns a list that contains the expression data
#' and the sample information. If there is no experiment file it
#' only returns the expression data.
#'
load_info <- function(
  expression_file,
  experiment_file,
  go_button,
  demo_data_file,
  demo_metadata_file
) {

  in_file_data <- expression_file
  in_file_data <- in_file_data$datapath

  if (is.null(in_file_data) && go_button == 0) {
    return(NULL)
  } else if (go_button > 0) {
    in_file_data <- demo_data_file
  }

  isolate({
    # Read expression file -----------
    data <- read.csv(in_file_data, quote = "", comment.char = "")
    # Tab-delimented if not CSV
    if (dim(data)[2] <= 2) {
      data <- read.table(
        in_file_data,
        sep = "\t",
        header = TRUE,
        quote = "",
        comment.char = ""
      )
    }

    # Filter out non-numeric columns ---------
    num_col <- c(TRUE)
    for (i in 2:dim(data)[2]) {
      num_col <- c(num_col, is.numeric(data[, i]))
    }
    if (sum(num_col) <= 2) {
      return(NULL)
    }
    data <- data[, num_col]

    # Format gene ids --------
    data[, 1] <- toupper(data[, 1])
    data[, 1] <- gsub(" |\"|\'", "", data[, 1])

    # Remove duplicated genes ----------
    data <- data[!duplicated(data[, 1]), ]

    # Set gene ids as rownames and get rid of column ---------
    rownames(data) <- data[, 1]
    data <- as.matrix(data[, c(-1)])

    # Remove "-" or "." from sample names ----------
    colnames(data) <- gsub("-", "", colnames(data))
    colnames(data) <- gsub("\\.", "", colnames(data))
  })

  # Read experiment file ----------
  in_file_expr <- experiment_file
  in_file_expr <- in_file_expr$datapath
  if (is.null(in_file_expr) && go_button == 0) {
    return(list(
      data = data
    ))
  } else if (go_button > 0) {
    sample_info_demo <- t(read.csv(
    demo_metadata_file,
    row.names = 1,
    header = T,
    colClasses = "character"
  ))
  return(list(
    sample_info = sample_info_demo,
    data = data
  ))
  }

  isolate({
    # Read experiment file ----------
    expr <- read.csv(
      in_file_expr,
      row.names = 1,
      header = TRUE,
      colClasses = "character"
    )
    if (dim(expr)[2] <= 2) {
      expr <- read.table(
        in_file_expr,
        row.names = 1,
        sep = "\t",
        header = TRUE,
        colClasses = "character"
      )
    }
    # remove "-" or "." from sample names ----------
    colnames(expr) <- gsub("-", "", colnames(expr))
    colnames(expr) <- gsub("\\.", "", colnames(expr))

    # Matching with column names of expression file ----------
    matches <- match(
      toupper(colnames(data)), toupper(colnames(expr))
    )
    matches <- matches[which(!is.na(matches))] # remove NA
    validate(need(
      length(unique(matches)) == dim(data)[2] &
        dim(expr)[1] >= 1 & dim(expr)[1] < 500,
      "Error!!! Sample information file not recognized. Sample names
       must be exactly the same. Each row is a factor. Each column
       represent a sample.  Please see documentation on format."
    ))

    # Check factor levels, change if needed ----------
    for (i in 1:dim(expr)[1]) {
      expr[i, ] <- gsub("-", "", expr[i, ])
      expr[i, ] <- gsub("\\.", "", expr[i, ])
    }

    # Factor levels match ---------
    if (length(unique(matches)) == dim(data)[2]) {
      expr <- expr[, matches]
      if (
        sum(apply(expr, 1, function(y) length(unique(y)))) >
        length(unique(unlist(expr)))) {
          factor_names <- apply(
            expr,
            2,
            function(y) paste0(names(y), y)
          )
          rownames(factor_names) <- rownames(expr)
          expr <- factor_names
        }
        return(list(
          data = data,
          sample_info = t(expr)
        ))
    } else {
      return(list(
        data = data
      ))
    }
  })
}