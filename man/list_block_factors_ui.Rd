% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_05_deg.R
\name{list_block_factors_ui}
\alias{list_block_factors_ui}
\title{Get the block factor choices}
\usage{
list_block_factors_ui(sample_info, select_factors_model)
}
\arguments{
\item{sample_info}{Experiment file information for grouping}

\item{select_factors_model}{The selected factors for the model
expression}
}
\value{
This function returns a vector of choices for a batch
 effect or paired samples.
}
\description{
This function uses the sample_info file and the selected
factors for the model to create a selection for the batch
effect. Returns a vector that turns into a checkbox for
the User.
}
