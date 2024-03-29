% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_03_clustering.R
\name{prep_download}
\alias{prep_download}
\title{Prep heatmap data for download}
\usage{
prep_download(heatmap, heatmap_data, cluster_meth = c(1, 2))
}
\arguments{
\item{heatmap}{Heatmap object from the \code{\link{heatmap_main}()} function.}

\item{heatmap_data}{Matrix of heatmap data from the
\code{\link{process_heatmap_data}()} function}

\item{cluster_meth}{Integer designating the clustering method used. 1 for
hierarchical and 2 for kmeans}
}
\value{
A dataframe with the heatmap data and associated clusters.
}
\description{
Prep heatmap data for download or additional analysis by merging gene ids
with the clusters from kmean clustering or hierarchical clustering.
}
