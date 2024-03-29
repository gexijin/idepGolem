% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_05_deg.R
\name{volcano_data}
\alias{volcano_data}
\title{Data processing for volcano and ma plot}
\usage{
volcano_data(
  select_contrast,
  comparisons,
  top_genes,
  limma_p_val,
  limma_fc,
  processed_data,
  contrast_samples,
  all_gene_names,
  select_gene_id
)
}
\arguments{
\item{select_contrast}{Comparison from DEG analysis to filter
for the significant genes}

\item{comparisons}{The comparisons vector from the results list
of the limma_value function}

\item{top_genes}{top_genes list from results list of the
limma_value function}

\item{limma_p_val}{Significant p-value to use to in determining
the expressed genes}

\item{limma_fc}{Minimum fold change value to use in determining
the expressed genes}

\item{processed_data}{Data matrix that has gone through
pre-processing}

\item{contrast_samples}{Samples that are included in the selected
comparison}

\item{all_gene_names}{Gene symbols, etc.}

\item{select_gene_id}{Selected gene id type, symbol, Ensmembl, etc.}
}
\value{
A list containing processed data for volcano plot in
 \code{plot_volcano()} & \code{plot_ma()} and list of differently expressed
 genes for labeling
}
\description{
Data processing for volcano and ma plot
}
