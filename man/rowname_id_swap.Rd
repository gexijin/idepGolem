% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_analysis_random.R
\name{rowname_id_swap}
\alias{rowname_id_swap}
\title{Swap rowname IDs for data matrix}
\usage{
rowname_id_swap(data_matrix, all_gene_names, select_gene_id)
}
\arguments{
\item{data_matrix}{Data matrix with ensembl or user gene ID rownames}

\item{all_gene_names}{Data frame of gene names}

\item{select_gene_id}{String designating the desired ID type for rownames
("User_ID", "ensembl_ID", "symbol")}
}
\value{
Data matrix with changed rownames
}
\description{
This function uses the all_gene_names dataframe
to swap the current rownames of a data matrix with
the desired gene ID. For instance, if the rownames
were currently ensembl, this function is able to
switch them back to the original form.
}
