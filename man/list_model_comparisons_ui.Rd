% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_05_deg.R
\name{list_model_comparisons_ui}
\alias{list_model_comparisons_ui}
\title{Create model comparisons choices}
\usage{
list_model_comparisons_ui(sample_info, select_factors_model, processed_data)
}
\arguments{
\item{sample_info}{Experiment file information for grouping}

\item{select_factors_model}{The selected factors for the model
expression}

\item{processed_data}{Data that has been through the pre-processing
function}
}
\value{
Returns a list containing a vector of choices and a
 title for the UI element.
}
\description{
This function uses the sample_info file and the selected
factors to create a list of options for model comparisons.
Changes with the input of select_factors_model. If there
is no selected factor then it defaults to comparisons that
can be created from the processed data.
}
