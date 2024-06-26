% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_06_pathway.R
\name{pgsea_data}
\alias{pgsea_data}
\title{Pathway analysis with the PGSEA package}
\usage{
pgsea_data(
  processed_data,
  gene_sets,
  my_range,
  pathway_p_val_cutoff,
  n_pathway_show
)
}
\arguments{
\item{processed_data}{Matrix of gene data that has been through
\code{\link{pre_process}()}}

\item{gene_sets}{List of vectors with each vector being the
set of genes that correspond to a particular pathway in
the database. See returned list from \code{\link{read_gene_sets}()}}

\item{my_range}{Vector of the (min_set_size, max_set_size)}

\item{pathway_p_val_cutoff}{Significant p-value to determine
enriched pathways}

\item{n_pathway_show}{Number of significant pathways to show}
}
\value{
A list with a data frame and a numeric value that is used
 in the \code{\link{plot_pgsea}()} to create a heatmap.
}
\description{
Run pathway analysis with the PGSEA package using the results
from the limma_value function.
}
\seealso{
Other pathway functions: 
\code{\link{fgsea_data}()},
\code{\link{find_overlap}()},
\code{\link{gage_data}()},
\code{\link{get_gsva_plot_data}()},
\code{\link{get_pathway_list_data}()},
\code{\link{get_pgsea_plot_all_samples_data}()},
\code{\link{get_pgsea_plot_data}()},
\code{\link{gsva_data}()},
\code{\link{kegg_pathway}()},
\code{\link{pathway_select_data}()},
\code{\link{pgsea_plot_all}()},
\code{\link{plot_gsva}()},
\code{\link{plot_pgsea}()},
\code{\link{reactome_data}()},
\code{\link{remove_pathway_id}()},
\code{\link{remove_pathway_id_second}()}
}
\concept{pathway functions}
