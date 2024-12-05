#' fct_01_load_data.R This file holds all of the main data analysis functions
#' associated with first tab of the iDEP website.
#'
#'
#' @section fct_01_load_data.R functions:
#'
#'
#' @name fct_01_load_data.R
NULL

#' Retrieve detailed info on genes
#'
#' This function retrieves detailed gene information from the
#' database for the matched species.
#'
#' @param converted List returned \code{\link{convert_id}()}. Contains
#'  information about the gene IDs for the matched species.
#' @param select_org String indicating selected organism for the expression
#'  data.Default is "BestMatch."
#' @param idep_data Data files from the database obtained with
#'  \code{\link{get_idep_data}()}
#'
#' @export
#' @return A data frame containing information on all ensembl IDs for the
#'  matched species. Usually a very large data frame due to the amount
#'  of IDs that the data base contains.
#'
#' @family load data functions
gene_info <- function(converted,
                      select_org,
                      idep_data) {
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

  conn_db <- connect_convert_db_org(
    select_org = select_org,
    idep_data = idep_data
  )

  # if database connection error
  # see ? try
  if(inherits(conn_db, "try-error")) {
    return(check)
  }

  # does the db file has a geneInfo table?
  ix <- grep(
    pattern = "geneInfo",
    x = DBI::dbListTables(conn_db)
  )

  check <- check_object_state(
    check_exp = (length(ix) == 0),
    true_message = as.data.frame("No matching gene info file found")
  )
  if (check$bool) {
    return(check)
  }

  check <- check_object_state(
    check_exp = (length(ix) != 1),
    true_message = as.data.frame("Multiple geneInfo file found!")
  )
  if (check$bool) {
    return(check)
  } else {
    query_statement <- "select * from geneInfo;"
    gene_info_csv <- DBI::dbGetQuery(conn_db, query_statement)
    DBI::dbDisconnect(conn_db)

# This is removed to make it faster. Make sure the original database
# is upper case, and there is no extra space character in Symbol
#    gene_info_csv[, 1] <- toupper(gene_info_csv[, 1])
    # remove the two spaces in gene symbol
 #   gene_info_csv$symbol <- gsub(
 #     " ",
 #     "",
 #     gene_info_csv$symbol
  #  )

  }

  set <- match(gene_info_csv$ensembl_gene_id, query_set)
  set[which(is.na(set))] <- "Genome"
  set[which(set != "Genome")] <- "List"
  return(cbind(gene_info_csv, set))
}


# Safe conversion function   "11,231" -> 11231
safe_numeric_conversion <- function(x) {
  converted <- suppressWarnings(as.numeric(gsub(",", "", x)))
  ifelse(is.na(converted), NA, converted)
}


#' Load basic data information
#'
#' This function does the immediate loading of the data and
#' sample info to present in the data preview table and the
#' sample info table. The data undergoes very basic filtering
#' and transformation before entering the table.
#'
#' @param expression_file Data file path for the expression file, should be
#'  accessed with \code{expression_file$datapath}
#' @param experiment_file Data file path for the experiment file, should be
#'  accessed with \code{experiment_file$datapath}
#' @param go_button TRUE/FALSE to load the demo data files
#' @param demo_data_file Expression demo data path \code{idep_data$demo_data_file}
#' @param demo_metadata_file Experiment demo data path
#'  \code{idep_data$demo_metadata_file}
#'
#' @export
#' @return This returns a list that contains the expression data
#' and the sample information. If there is no experiment file it
#' only returns the expression data.
#'
#' @family load data functions
input_data <- function(expression_file,
                       experiment_file,
                       go_button,
                       demo_data_file,
                       demo_metadata_file) {
  in_file_data <- expression_file
  in_file_data <- in_file_data$datapath

  if (is.null(in_file_data) && go_button == 0) {
    return(NULL)
  } else if (go_button > 0) { # use demo data
    in_file_data <- demo_data_file
  }
  isolate({
    # Read expression file -----------

    file_extension <- tolower(tools::file_ext(in_file_data))
    if ( file_extension == "xlsx" || 
         file_extension == "xls"
      )  {
      data <- readxl::read_excel(in_file_data)
      data <- data.frame(data)
    } else {
       data <- read.csv(in_file_data,
         header = TRUE, stringsAsFactors = FALSE,
         quote = "\"", comment.char = "",
         blank.lines.skip = TRUE
       )
    }

    # Try tab-delimented if not CSV
    if (ncol(data) <= 1) {
      data <- read.table(
        in_file_data,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        quote = "\"",
        comment.char = "",
        blank.lines.skip = TRUE
      )
    }

    # try semicolon -delimented if not CSV
    if (ncol(data) <= 1) {
      data <- read.table(
        in_file_data,
        sep = ";",
        header = TRUE,
        stringsAsFactors = FALSE,
        quote = "\"",
        comment.char = "",
        blank.lines.skip = TRUE
      )
    }

    # try space-delimented if not CSV
    if (ncol(data) <= 1) {
      data <- read.table(
        in_file_data,
        sep = " ",
        header = TRUE,
        stringsAsFactors = FALSE,
        quote = "\"",
        comment.char = "",
        blank.lines.skip = TRUE
      )
    }

    
    # if more than one column
    if (ncol(data) > 1) {
      # Convert all columns after the first one to numeric using the safe function
      data[, -1] <- lapply(data[, -1], function(col) {
        if (is.character(col)) {
          return(sapply(col, safe_numeric_conversion))
        } else {
          return(col)
        }
      })

      # rows where all values after the first column are NA cause issues
      # if a row has all NA values down stream processing fails
      # Remove rows where all values after the first column are NA
      data <- data[!apply(is.na(data[, -1]), 1, all), ]

      # Identify rows where all values after the first column are NA
      #all_na_rows <- apply(is.na(data[, -1]), 1, all)
      # Change all NA values in those rows to 0, column by column
      #data[all_na_rows, -1] <- lapply(
      #  data[all_na_rows, -1],
      #  function(col) ifelse(is.na(col), 0, col)
      #)

      # remove columns where all values are NA
      data <- data[, !apply(is.na(data), 2, all)]

      # remove columns where all values are 0
      ix <- apply(data, 2, function(col) all(col == 0))
      if (sum(ix) > 0) {
        showNotification(
          ui = paste0("Warning!!! Columns with all zero values are deleted: ", 
                      paste0(names(data)[ix], collapse = ", ")),
          id = "zero_expression_file",
          duration = 10,
          type = "warning"
        )
        data <- data[, !ix]        
      }
    }

    # cannot parse file; only one or two column
    if (ncol(data) <= 1) {
      showNotification(
        ui = "Error!!! Expression file not recognized. 
        Click the Reset button, examine the file, and try again.",
        id = "error_expression_file",
        duration = NULL,
        type = "error"
      )
      return(NULL)
    }
    # Order by SD ----------
    #data <- data[order(-apply(
    #  data[, 2:ncol(data)],
    #  1,
    #  sd
    #)), ]

    # Format gene ids --------
    #iconv converts latin1 to UTF-8; otherwise toupper(ENSG00000267023Ê) causes error
    data[, 1] <- toupper(iconv(data[, 1], "latin1", "UTF-8"))
    data[, 1] <- gsub(" |\"|\'", "", data[, 1])

    # Remove duplicated genes ----------
    data <- data[!duplicated(data[, 1]), ]

    # Remove rows without genes IDs----------
    data <- data[!is.na(data[, 1]), ]

    # Set gene ids as rownames and get rid of column ---------
    rownames(data) <- data[, 1]
    data <- as.matrix(data[, c(-1)])

    # use janitor to clean up column names;  too slow!!!
    #data <- janitor::clean_names(data)

    # Remove "-" or "." from sample names ----------
    colnames(data) <- gsub("-", "", colnames(data))
    colnames(data) <- gsub("\\.", "", colnames(data))
  })

  # Read experiment file ----------
  in_file_expr <- experiment_file
  in_file_expr <- in_file_expr$datapath
  if (is.null(in_file_expr) && go_button == 0) {
    return(list(
      data = data,
      sample_info = NULL
    ))
  } else if (go_button > 0) {
    sample_info_demo <- NULL
    # if design file is not ""
    if (!is.null(demo_data_file)) {
      if (nchar(demo_metadata_file) > 2) {
        sample_info_demo <- t(read.csv(
          demo_metadata_file,
          row.names = 1,
          header = TRUE,
          colClasses = "character"
        ))
      }
    }
    return(list(
      sample_info = sample_info_demo,
      data = data
    ))
  }

  isolate({
    # Read experiment file ----------
    file_extension <- tolower(tools::file_ext(in_file_expr))
    if ( file_extension == "xlsx" || 
         file_extension == "xls"
      )  {
      expr <- readxl::read_excel(in_file_expr)
      expr <- data.frame(expr)
    } else {
      expr <- read.csv(
        in_file_expr,
        row.names = 1,
        header = TRUE,
        colClasses = "character",
        quote = "\"",
        comment.char = "",
        blank.lines.skip = TRUE
      )
    }

    # Try tab-delimented if not CSV
    if (ncol(expr) <= 1) {
      expr <- read.table(
        in_file_expr,
        row.names = 1,
        sep = "\t",
        header = TRUE,
        colClasses = "character",
        quote = "\"",
        comment.char = "",
        blank.lines.skip = TRUE
      )
    }


    # remove "-" or "." from sample names ----------
    colnames(expr) <- gsub("-", "", colnames(expr))
    colnames(expr) <- gsub("\\.", "", colnames(expr))

    rownames(expr) <- gsub("-", "", rownames(expr))
    rownames(expr) <- gsub("\\.", "", rownames(expr))
    rownames(expr) <- gsub(" ", "", rownames(expr))

    # Matching with column names of expression file ----------
    matches <- match(
      toupper(colnames(data)), toupper(colnames(expr))
    )
    matches <- matches[which(!is.na(matches))] # remove NA

    validate(need(
      length(unique(matches)) == ncol(data) &&
        nrow(expr) >= 1 && nrow(expr) < 500,
      "Error!!! Sample information file not recognized. Column names
       must be exactly the same. Each row is a factor. Each column
       represent a sample.  Please see documentation on format.
       Please click on the Reset button and try again."
    ))

    # Check factor levels, change if needed ----------
    ignored_factors <- c()
    for (i in 1:nrow(expr)) {
      expr[i, ] <- gsub("-", "", expr[i, ])
      expr[i, ] <- gsub("\\.", "", expr[i, ])
      # Remove whitespace
      expr[i, ] <- gsub(" ", "", expr[i, ])
      # Convert to upper case to avoid mixing
      expr[i, ] <- toupper(expr[i, ])
      if(length(unique(t(expr[i, ]))) == 1) {
        ignored_factors <- c(ignored_factors, i)
      }
    }

    if(length(ignored_factors) == nrow(expr)) {
      # ignore design file.
      matches <- 1
      showNotification(
        ui = "Ignoring design file. All rows have one unique levels.",
        id = "ignore_all_rows",
        duration = NULL,
        type = "default"
      )
    } else if(length(ignored_factors) > 0) {    
      # delete the row with only one level
      ignored_rows <- paste0(row.names(expr)[ignored_factors], collapse = ", ")
      showNotification(
        ui = paste("Ignoring rows with one unique levels:", ignored_rows),
        id = "ignore_some_rows",
        duration = NULL,
        type = "default"
      )
      expr <- expr[-1 * ignored_factors, ] 
    }

    # Factor levels match ---------
    if (length(unique(matches)) == ncol(data)) {
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

#' Convert expression data rownames to ensembl
#'
#' This function works with the converted ids to take in the
#' expression data and swap the rownames to the converted ids.
#' It returns the data exactly the same except for the changed
#' ids.
#'
#' @param converted Data from \code{\link{convert_id}()} function containing
#'  converted ids
#' @param no_id_conversion TRUE/FALSE for converting data ids or not
#' @param data Data from the input expression file
#' @param multiple_map String to designate how to handle values when multiple
#'   ids are matched to the same gene. Options are 'mean', 'median', 'sum',
#'   'max', and 'max_id'.
#'
#' @export
#' @return Returns original data with rownames converted to ensembl
#'
#' @family load data functions
convert_data <- function(converted,
                         data,
                         no_id_conversion,
                         multiple_map =
                           c("mean", "median", "sum", "max", "max_sd")) {
  if (is.null(converted) || no_id_conversion) {
    return(list(
      data = data,
      mapped_ids = rownames(data)
    ))
  } else {
    mapping <- converted$conversion_table

    rownames(data) <- toupper(rownames(data))
    merged <- merge(
      mapping[, 1:2],
      data,
      by.y = "row.names",
      by.x = "User_input",
      all.y = TRUE
    )
    no_match <- which(is.na(merged[, 2]))
    merged[no_match, 2] <- merged[no_match, 1]

    mapped_ids <- merged[, 1:2]

    # remove user id
    merged <- merged[, -1]

    if (multiple_map == "max_sd") {
      # Multiple matches use one with highest SD ----------
      tmp <- apply(merged[, 2:(ncol(merged))], 1, sd)
      merged <- merged[order(merged[, 1], -tmp), ]
      merged <- merged[!duplicated(merged[, 1]), ]
    } else {
      # aggregate multiple ids mapped to the same gene
      merged <- switch(multiple_map,
        "sum" = aggregate(
          . ~ ensembl_gene_id,
          data = merged,
          FUN = sum_fun,
          #na.rm = TRUE, # if used, missing values are replaced with NA.
          na.action = na.pass  # otherwise, missing values --> genes removed
        ),
        "mean" = aggregate(
          . ~ ensembl_gene_id,
          data = merged,
          FUN = mean_fun,
          na.action = na.pass
        ),
        "max" = aggregate(
          . ~ ensembl_gene_id,
          data = merged,
          FUN = max_fun,
          na.action = na.pass
        ),
        "median" = aggregate(
          . ~ ensembl_gene_id,
          data = merged,
          FUN = median_fun,
          na.action = na.pass
        )
      )
    }

    rownames(merged) <- merged[, 1]
    merged <- as.matrix(merged[, c(-1)])

    # Order by SD ----------
    merged <- merged[order(-apply(
      merged[, 1:ncol(merged)],
      1,
      function(x) sd(x, na.rm = TRUE)
    )), ]

    return(list(
      data = merged,
      mapped_ids = mapped_ids
    ))
  }
}

#' Get all matched gene names
#'
#' This function will create a data frame with all the matched
#' gene names from the database.
#'
#' @param mapped_ids Matched IDs from \code{\link{convet_id}()}
#' @param all_gene_info Gene information matched species in idep database from
#'  \code{\link{gene_info}()}
#'
#'
#' @export
#' @return Data frame containing all the matched ID names from idep
#' database. Three columns denotes a recognized species for which
#' idep had gene names for. Two columns means the IDs were converted
#' to ensembl format, but no species was found for the gene names.
#' One means no conversion occurred.
#'
#' @family load data functions
get_all_gene_names <- function(mapped_ids,
                               all_gene_info) {
  # not converted
  if (is.null(dim(mapped_ids))) {
    return(data.frame(
      "User_ID" = mapped_ids,
      "ensembl_ID" = mapped_ids, # dummy data
      "symbol" = mapped_ids # dummy data
    ))
  } else if (!is.null(all_gene_info$bool)) { # ensembl ID only, no symbol
    return(data.frame(
      "User_ID" = mapped_ids[, 1],
      "ensembl_ID" = mapped_ids[, 2],
      "symbol" = mapped_ids[, 1] # dummy data
    ))
  } else {
    mapped_ids <- data.frame(
      "User_ID" = mapped_ids[, 1],
      "ensembl_ID" = mapped_ids[, 2]
    )
    all_names <- merge(
      mapped_ids,
      all_gene_info[, c("ensembl_gene_id", "symbol")],
      by.x = "ensembl_ID",
      by.y = "ensembl_gene_id",
      all.x = T
    )
    all_names <- all_names[!duplicated(all_names$ensembl_ID), ]
    all_names$symbol[all_names$symbol == ""] <- NA
    all_names$symbol[is.na(all_names$symbol)] <- {
      all_names$ensembl_ID[is.na(all_names$symbol)]
    }
    duplicates <- all_names$symbol %in%
      all_names$symbol[duplicated(all_names$symbol)]
    all_names$symbol[duplicates] <- paste(
      all_names$symbol[duplicates],
      all_names$ensembl_ID[duplicates]
    )
    all_names <- dplyr::select(
      all_names,
      User_ID,
      ensembl_ID,
      symbol,
      tidyselect::everything()
    )

    return(all_names)
  }
}

#' Calculate sum of vector
#'
#' This function refines the sum function to handle missing value, used 
#' with the aggregate function
#'
#' @param x a vector 
#'
#' @export
#' @return a number, or NA if all the numbers are NA
#'
sum_fun <- function(x){
  if(sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    sum(x, na.rm = TRUE)
  }
}

#' Calculate max of vector
#'
#' This function refines the max function to handle missing value, used 
#' with the aggregate function
#'
#' @param x a vector 
#'
#' @export
#' @return a number, or NA if all the numbers are NA
#'
max_fun <- function(x){
  if(sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    max(x, na.rm = TRUE)
  }
}


#' Calculate mean of vector
#'
#' This function refines the mean function to handle missing value, used 
#' with the aggregate function
#'
#' @param x a vector 
#'
#' @export
#' @return a number, or NA if all the numbers are NA
#'
mean_fun <- function(x){
  if(sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    mean(x, na.rm = TRUE)
  }
}

#' Calculate median of vector
#'
#' This function refines the median function to handle missing value, used 
#' with the aggregate function
#'
#' @param x a vector 
#'
#' @export
#' @return a number, or NA if all the numbers are NA
#'
median_fun <- function(x){
  if(sum(is.na(x)) == length(x)) {
    return(NA)
  } else {
    median(x, na.rm = TRUE)
  }
}




#' Show sample gene IDs in modal
#'
#' Show example gene IDs for a selected species.
#'
#' @param species String indicating selected organism for the expression
#' @param db String indicating selected database
#' @param nGenes Number of genes to show for each idType
#'
#' @export
#' @return a data frame with example gene IDs
#'
showGeneIDs <- function(species, db, nGenes = 10){
  # Given a species ID, this function returns 10 gene ids for each idType
  if(species == "BestMatch")
    return(as.data.frame("Select a species above.") )

  converted <- NULL
  try(
    converted <- DBI::dbConnect(
      drv = RSQLite::dbDriver("SQLite"),
      dbname = paste0(DATAPATH, "/db/", db),
      flags = RSQLite::SQLITE_RO #read only mode
    ),
    silent = TRUE
  )

  if(is.null(converted)){
    showNotification(
      ui = paste("Selected database is not downloaded"),
      id = "db_notDownloaded",
      duration = 2.5,
      type = "error"
    )
    return()
  }
  removeNotification("db_notDownloaded")
  showNotification(
    ui = paste("Querying Data...  May take 5 minutes."),
    id = "ExampleIDDataQuery",
    duration = NULL,
    type = "message"
  )

  idTypes <- DBI::dbGetQuery(
    conn = converted,   
    paste0( 
      "WITH RandomIds AS (
      SELECT m.idType,
           m.id,
           ROW_NUMBER() OVER (PARTITION BY m.idType ORDER BY RANDOM()) AS rn
      FROM Mapping m
      )
      SELECT i.*, r.id AS RandomId
      FROM idIndex i
      LEFT JOIN RandomIds r ON i.id = r.idType AND r.rn <= ", nGenes, ";"
    )
  )
  DBI::dbDisconnect(converted) # suggested by GitHub Copilot

  result <- aggregate(
    RandomId ~ idType, 
    data = idTypes,
    FUN = function(x) paste(x, collapse = "; ")
  )
  colnames(result) <- c("ID Type", "Examples")
  
  # put symbols first, refseq next, followed by ensembls. Descriptions (long gnee names) last
  result <- result[ order( grepl("ensembl", result$'ID Type'), decreasing = TRUE), ]
  result <- result[ order( grepl("refseq", result$'ID Type'), decreasing = TRUE), ]
  result <- result[ order( grepl("symbol", result$'ID Type'), decreasing = TRUE), ]
  result <- result[ order( grepl("description", result$'ID Type'), decreasing = FALSE), ]

  return(result)
}

### Add Example Gene ID column to database

# -- Create a temporary table to store the concatenated Example IDs
# CREATE TEMPORARY TABLE TempExampleIds AS
# WITH RandomIds AS (
#   SELECT m.idType,
#   m.id,
#   ROW_NUMBER() OVER (PARTITION BY m.idType ORDER BY RANDOM()) AS rn
#   FROM Mapping m
# )
# SELECT i.id, GROUP_CONCAT(DISTINCT r.id) AS "ExampleIDs"
# FROM idIndex i
# LEFT JOIN RandomIds r ON i.id = r.idType AND r.rn <= 100
# GROUP BY i.id;
# 
# -- Add new column
# -- ALTER TABLE idIndex
# -- ADD ExampleIDs TEXT;
# 
# -- Fill column with example genes
# UPDATE idIndex AS i
# SET "ExampleIDs" = (
#   SELECT ExampleIDs
#   FROM TempExampleIds AS t
#   WHERE t.id = i.id
# );
# 
# -- Drop the temporary table
# DROP  TABLE TempExampleIds;
# 
# -- View results
# SELECT * FROM idIndex



### Policies' Content ###
terms_of_use_content <- function() {
  HTML('
    <div class="policy">
      <h1>iDEP Terms of Service</h1>
      <p style="font-weight: bold;">Updated Aug. 19th, 2024</p>
      <p>Welcome to Orditus LLC ("we", "us", or "our"). These Terms of Use ("Terms") govern your use of our website RTutor.ai and any services, content, or features provided by us (collectively, the "Services"). By accessing or using our Services, you agree to be bound by these Terms. If you do not agree to these Terms, you may not use our Services.</p>
      <h2>1. Description of Services</h2>
      <p>RTutor is designed to generate and run R and Python code. Using the power of ChatGPT, RTutor translates your natural language into scripts to analyze and interpret data.</p>
      <h2>2. Acceptance of Terms</h2>
      <p>By accessing and using our Services, you accept and agree to be bound by these Terms and our Privacy Policy. If you do not agree to these Terms, you should not use our Services.</p>
      <h2>3. Changes to Terms</h2>
      <p>We reserve the right to modify these Terms at any time. Any changes will be effective immediately upon posting of the revised Terms. Your continued use of the Services after the posting of any changes constitutes your acceptance of the revised Terms.</p>
      <h2>4. User Obligations and Restrictions</h2>
      <p>You agree to use the Services only for lawful purposes and in a manner that does not infringe the rights of, restrict, or inhibit anyone else\'s use and enjoyment of the Services. Prohibited behavior and use includes but is not limited to:</p>
      <p>You must be at least 13 years old to use this Service.<br>
      You must not use this Service commercially in any form without our written consent.<br>
      You agree not to transmit obscene or offensive content.<br>
      You agree not to use this Service to generate content that is misleading, harmful, or illegal.<br>
      You must not overuse the Services in a manner that adversely impacts the performance or availability of the Services for other users.<br>
      You must not upload data that is confidential or personally identifiable without proper anonymization.<br>
      You agree not to attempt to reverse-engineer or otherwise tamper with the AI models and algorithms.<br>
      You must not copy or scrape the content of the Service in any form without our written consent.</p>
      <h2>5. User Content</h2>
      <p>By submitting, posting, or displaying content, you grant us a worldwide, non-exclusive, royalty-free license to use, reproduce, adapt, modify, publish, and distribute such content in any and all media. You are responsible for any data you submit to our Services. By submitting data, you agree that:<br><br>
      You have the necessary rights and permissions to share the data.<br>
      You will not submit any data that violates the rights of others or any applicable laws.<br>
      We may use, store, and process basic information about your data as described in our Privacy Policy.</p>
      <h2>6. Intellectual Property</h2>
      <p>All content, trademarks, and data on this website, including but not limited to software, databases, text, graphics, icons, hyperlinks, private information, designs, and agreements, are the property of or licensed to Orditus LLC and as such are protected from infringement by local and international legislation and treaties. All content, AI models, algorithms, and data used or generated by our Services are the property of Orditus LLC or its licensors and are protected by intellectual property laws. By using our Services, you agree that:
      <br>You retain ownership of the data you submit to the AI assistant.
      <br>You grant us a license to use, reproduce, adapt, and modify your inputs as necessary to provide and improve our Services.
      <br>We retain ownership of all models, algorithms, and enhancements made to the Services based on aggregated and anonymized user data.</p>
      <h2>7. Limitation of Liability</h2>
      <p>To the fullest extent permitted by law, Orditus LLC shall not be liable for any indirect, incidental, special, consequential, or punitive damages arising from or related to your use of the Services, including but not limited to:<br><br>
      <br>Inaccuracies or errors in AI-generated results.
      <br>Loss or corruption of data.
      <br>Misuse of AI-generated insights for decision-making.
      <br>Unauthorized access to or use of your account or data.<br>
      In no event shall Orditus LLC nor its directors, employees, partners, or affiliates, be liable for any indirect, incidental, special, consequential, or punitive damages, including without limitation, loss of profits, data, use, goodwill, or other intangible losses, resulting from:<br>
      Your use or inability to use the Services.<br>
      Any unauthorized access to or use of our servers and/or any personal information stored therein.<br>
      Any breaks or stoppages in our service that affect the sending or receiving of data to or from the Services.</p>
      <h2>8. Termination</h2>
      <p>We may terminate or suspend your access to our Services immediately, without prior notice or liability, for any reason whatsoever, including without limitation if you breach the Terms. Upon termination, your right to use the Services will immediately cease.</p>
      <h2>9. Privacy Policy</h2>
      <p>Your privacy is important to us. Please review our Privacy Policy to understand how we collect, use, and protect your information. We may use this information to improve the Service. By using our Services, you consent to the data practices described in our Privacy Policy.</p>
      <h2>10. Governing Law</h2>
      <p>These Terms shall be governed and construed in accordance with the laws of South Dakota, USA without regard to its conflict of law provisions.</p>
      <h2>11. Contact Us</h2>
      <p style="padding-bottom: 90px;">If you have any questions about these Terms, please contact us at ge@orditus.com.</p>
    </div>
  ')
}

privacy_policy_content <- function() {
  HTML('
    <div class="policy">
      <h1>iDEP Privacy Policy</h1>
      <p style="font-weight: bold;">Updated Aug. 19th, 2024</p>
      <p>We are committed to protecting your privacy. This Privacy Policy explains how ("we", "us", or "our"), Orditus LLC, collects, uses, discloses, and safeguards your information when you use our website RTutor.ai and our services (collectively, the "Services"). By accessing or using the Services, you agree to the terms of this Privacy Policy. If you do not agree with the terms of this Privacy Policy, please do not access or use the Services.</p>
      <h2>1. Information We Collect</h2>
      <h3>1.1 Personal Information</h3>
      <p>We may collect personal information that you provide directly to us when you use our Services or communicate with us. This information may include:<br><br>
      Name<br>
      Email address<br>
      Company name<br><br></p>
      <h3>1.2 Usage Information</h3>
      <p>We may collect information from third parties that your browser sends whenever you visit our website or use our Services. This data may include:<br><br>
      Browser type,<br>
      Browser version,<br>
      Pages of our Service that you visit,<br>
      City of user,<br>
      Time and date of your visit,<br>
      Time spent on those pages,<br>
      Other diagnostic data.<br><br>
      We do not collect IP addresses or demographic data.</p>
      <h3>1.3 Upload Information and Prompts</h3>
      <p>When you upload data to our Services, it is stored temporarily for the duration of your session. Once your session ends or you close your browser, all uploaded data is automatically deleted. We do not track or store any of the data you upload.
      <br>By default, none of your uploaded data is sent to the selected language learning model. However, if you choose, there is a settings option that allows our Services to randomly select at most five rows from your data to send to the selected Large Language Model (LLM) to enhance coding results.
      <br>We do collect data structure information and the prompts you send while using our Services. If you prefer not to share this information, you can opt out through the settings.
      <br>Additionally, column names and data type information from your uploaded data are sent to the selected language learning model to improve the functionality of the Services. This includes category names of all factor variables within your uploaded data.</p>
      <h2>2. How We Collect This Data</h2>
      <p>To collect data structure information and prompts, we utilize Google Analytics and an SQL database. Google Analytics helps us track and analyze user interactions with our Services, providing insights into how data is structured, and which prompts are most commonly used. This information is collected in real-time and stored securely in our SQLite database. The SQLite database allows us to manage and retrieve this data efficiently, ensuring that it is used to improve the quality and functionality of our Services.</p>
      <h2>3. Why We Collect This Data</h2>
      <p>We may use information from users for the following purposes:<br><br>
      To operate and improve the Services by providing you with more effective customer service.<br>
      To perform research and analysis aimed at improving the Services, or other products and technologies of Orditus LLC.<br>
      To diagnose or fix problems with the Services.<br>
      To communicate with you, especially through the sending of emails, information on technical Service issues, security announcements, information on new Services available, legal notices, response to your requests, or other information that we think might be relevant to you.<br>
      To protect against harm to the rights, property or safety of Orditus LLC, our users, yourself, or the public.<br>
      To enforce any terms of use, including investigations of potential violations of the terms.<br>
      To comply with any law, regulation, legal process, or governmental requests.</p>
      <h2>4. How We Use This Data</h2>
      <p>Orditus LLC will use the information collected from you to pursue legitimate interests such as legal or regulatory compliance, security control, business operations, scientific research, or any other interest reasonably held as legitimate.</p>
      <h2>5. Sharing This Data</h2>
      <p>We will not lease, sell, or rent your personal, usage, or upload information except in the following instances described in this policy.<br></p>
      <h3>5.1 Third Party Service Providers</h3>
      <p>We may occasionally hire other companies to provide limited services on our behalf, such as providing customer support, hosting websites, or performing statistical analysis of our Services. Those companies will be allowed to obtain only the information they need to deliver the relevant service. They will be required to maintain the confidentiality of the information and are prohibited from using it for any other purpose.<br><br>
      Third Party Service Providers:<br>
      Google Analytics<br>
      SQLite<br>
      AWS<br></p>
      <h3>5.2 With your Consent</h3>
      <p>At any time during your use of our Services, or upon explicit request from us, you may consent to the disclosure of your information.<br></p>
      <h3>5.3 For Security and Safety Purposes</h3>
      <p>In the event of any fraud, security threats, or incidents, we reserve the right to disclose your information without your consent for the purposes of maintaining security on our website and for our users and addressing fraud or security issues. We reserve the right to disclose your information without your consent for the purpose of protecting against harm to the rights, property, or safety of Orditus LLC, our users, yourself, or the public.<br></p>
      <h3>5.4 For Legal or Regulatory Purposes</h3>
      <p>We reserve the right to disclose your information without your consent to comply with any applicable law, regulation, legal process, or governmental requests.</p>
      <h2>6. Opting Out of Sharing Your Data</h2>
      <h3>6.1 Access your Information</h3>
      <p>You may be entitled under data protection laws to access and review personal information that Orditus LLC holds related to you. You may access, modify, or delete the information we collected by controlling the content that you share or upload at any time.<br><br>
      Any other request should be addressed to: ge@orditus.com. Such inquiries should be clearly marked as data protection queries, and you should indicate if the request is time sensitive.<br></p>
      <h3>6.2 Data Retention</h3>
      <p>We retain your Information for as long as necessary to deliver our Services, to comply with any applicable legal requirements, to maintain security and prevent incidents, and, in general, to pursue our legitimate interests. If you wish to request the erasure of all your personal information that we process, you may do so by sending a written request to ge@orditus.com.</p>
      <h3>6.3 Opting Out</h3>
      <p>You may opt out of sending information relevant to your data uploads and prompt requests at any time by switching off the checkbox in the RTutor.ai settings.</p>
      <h2>7. Data Security</h2>
      <p>All interactions with our Services occur over HTTPS which ensures that data transmitted between your device and our servers is encrypted and secure.</p>
      <h2>8. Children’s Privacy</h2>
      <p>Orditus.com is neither directed to nor structured to attract children under the age of 13 years. Accordingly, Orditus LLC does not intend to collect personal information from anyone it knows to be under 13 years of age. Orditus LLC will direct potential users under 13 years of age not to use the Services. If the Company learns that personal information of persons less than 13 years of age has been collected without verifiable parental consent, Orditus LLC will take the appropriate steps to delete this information.<br><br>
      To make such a request, please contact Orditus LLC at: ge@orditus.com.</p>
      <h2>9. Contact Us</h2>
      <p style="padding-bottom: 90px;">If you have questions about this policy, please email us at: ge@orditus.com.</p>
  </div>
  ')
}

