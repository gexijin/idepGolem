% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_analysis_random.R
\name{basic_heat_sub}
\alias{basic_heat_sub}
\title{Heatmap of User brush selection}
\usage{
basic_heat_sub(ht_brush, ht, ht_pos_main, heatmap_data)
}
\arguments{
\item{ht_brush}{Brush information from the User on the main
heatmap}

\item{ht}{Main heatmap to create the sub-heatmap from}

\item{ht_pos_main}{Position information from the main heatmap
to use for the sub-heatmap}

\item{heatmap_data}{Data matrix that is being plotted in the
main heatmap}
}
\value{
A ComplexHeatmap object that will be inputted into the
 draw function in the server, the sub-heatmap data matrix, the
 group color mapping for the annotation, and the groups that
 the columns fall into.
}
\description{
Create a heatmap from the brush selection of the main heatmap.
Used in iDEP to create an interactive heatmap and enable the
User to zoom in on areas they find interesting.
}
