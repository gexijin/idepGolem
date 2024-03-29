% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_01_load_data.R
\name{input_data}
\alias{input_data}
\title{Load basic data information}
\usage{
input_data(
  expression_file,
  experiment_file,
  go_button,
  demo_data_file,
  demo_metadata_file
)
}
\arguments{
\item{expression_file}{Data file path for the expression file, should be
accessed with \code{expression_file$datapath}}

\item{experiment_file}{Data file path for the experiment file, should be
accessed with \code{experiment_file$datapath}}

\item{go_button}{TRUE/FALSE to load the demo data files}

\item{demo_data_file}{Expression demo data path \code{idep_data$demo_data_file}}

\item{demo_metadata_file}{Experiment demo data path
\code{idep_data$demo_metadata_file}}
}
\value{
This returns a list that contains the expression data
and the sample information. If there is no experiment file it
only returns the expression data.
}
\description{
This function does the immediate loading of the data and
sample info to present in the data preview table and the
sample info table. The data undergoes very basic filtering
and transformation before entering the table.
}
\seealso{
Other load data functions: 
\code{\link{convert_data}()},
\code{\link{gene_info}()},
\code{\link{get_all_gene_names}()}
}
\concept{load data functions}
