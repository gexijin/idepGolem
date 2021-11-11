#' fct_08_bicluster.R This file holds all of the main data analysis functions
#' associated with eighth tab of the iDEP website.
#'
#'
#' @section fct_07_bicluster.R functions:
#'
#'
#' @name fct_08_bicluster.R
NULL

#' BICLUSTERING MAIN FUNCTION
get_biclustering <- function(
  data,
  n_genes,
  biclust_method
) {
  if(n_genes > dim(data)[1]) {
    n_genes <- dim(x)[1]
  }
  if(n_genes < 10) {
    n_genes <- 10
  } 
  if(n_genes > 2000 ) {
    n_genes <- 2000
  } 			
  data <- as.matrix(data[1:n_genes, ]) - apply(data[1:n_genes, ], 1, mean)
			
  if(biclust_method == "BCXmotifs()") {
    data <- biclust::discretize(data)
  }
  
  run_biclust <- paste(
    "res <- biclust::biclust(as.matrix(data), method = ", biclust_method, ")"
  )
  eval(parse(text = run_biclust))

  return(list(
    data = data,
    res = res
  ))	 
}

bicluster_summary_message <- function(
  biclustering,
  select_bicluster
) {
  res <- biclustering$res
	if(res@Number == 0) { 
		return(
      "No cluster found! Perhaps sample size is too small."
    ) 
	}	else  { 
		data <- biclust::bicluster(
      biclustering$data,
      res,
      as.numeric(select_bicluster)
    )[[1]]
    return(
      paste(
        res@Number,
        "clusters found. Cluster ",
        select_bicluster,
        " has",
        dim(data)[1],
        "genes correlated across",
        dim(data)[2],
        "samples."
      )	
    )
	}
				

		
}