% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_06_pathway.R
\name{get_pgsea_plot_all_samples_data}
\alias{get_pgsea_plot_all_samples_data}
\title{Data from PGSEA plot all samples}
\usage{
get_pgsea_plot_all_samples_data(
  data,
  select_contrast,
  gene_sets,
  my_range,
  pathway_p_val_cutoff,
  n_pathway_show
)
}
\arguments{
\item{data}{Matrix of gene data that has been through
\code{\link{pre_process}()}}

\item{select_contrast}{String designating the comparison from DEG analysis to
filter for the significant genes. See the 'comparison' element from the list
returned from \code{\link{limma_value}()} for options.}

\item{gene_sets}{List of vectors with each vector being the
set of genes that correspond to a particular pathway in
the database. See list returned from \code{\link{read_gene_sets}()}}

\item{my_range}{Vector of the (min_set_size, max_set_size)}

\item{pathway_p_val_cutoff}{Significant p-value to determine
enriched pathways}

\item{n_pathway_show}{Number of pathways to return in final
result}
}
\value{
Data matrix with the rownames the descriptions of pathways
 and the matrix the returned expression calculation from the PGSEA
 package.
}
\description{
Get the data matrix that is plotted in the heatmap created by
the pgsea_plot_all function.
}
\seealso{
Other pathway functions: 
\code{\link{fgsea_data}()},
\code{\link{find_overlap}()},
\code{\link{gage_data}()},
\code{\link{get_gsva_plot_data}()},
\code{\link{get_pathway_list_data}()},
\code{\link{get_pgsea_plot_data}()},
\code{\link{gsva_data}()},
\code{\link{kegg_pathway}()},
\code{\link{pathway_select_data}()},
\code{\link{pgsea_data}()},
\code{\link{pgsea_plot_all}()},
\code{\link{plot_gsva}()},
\code{\link{plot_pgsea}()},
\code{\link{reactome_data}()},
\code{\link{remove_pathway_id}()},
\code{\link{remove_pathway_id_second}()}
}
\concept{pathway functions}
