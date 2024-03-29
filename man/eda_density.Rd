% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_02_pre_process.R
\name{eda_density}
\alias{eda_density}
\title{Density plot for the processed data}
\usage{
eda_density(processed_data, sample_info, plots_color_select)
}
\arguments{
\item{processed_data}{Matrix of gene data that has gone through
\code{\link{pre_process}()}}

\item{sample_info}{Sample_info from the experiment file}

\item{plots_color_select}{Vector of colors for plots}
}
\value{
A density plot as a \code{ggplot} object
}
\description{
This function takes in the processed data and sample info
and creates a density plot for the distribution of sequences
that are mapped to each sample.
}
\seealso{
Other plots: 
\code{\link{PCA_plot}()},
\code{\link{chr_counts_ggplot}()},
\code{\link{chr_normalized_ggplot}()},
\code{\link{cor_plot}()},
\code{\link{draw_sample_tree}()},
\code{\link{eda_boxplot}()},
\code{\link{eda_scatter}()},
\code{\link{gene_counts_ggplot}()},
\code{\link{individual_plots}()},
\code{\link{k_means_elbow}()},
\code{\link{mean_sd_plot}()},
\code{\link{rRNA_counts_ggplot}()},
\code{\link{sd_density}()},
\code{\link{t_SNE_plot}()},
\code{\link{total_counts_ggplot}()}

Other preprocess functions: 
\code{\link{chr_counts_ggplot}()},
\code{\link{chr_normalized_ggplot}()},
\code{\link{eda_boxplot}()},
\code{\link{eda_scatter}()},
\code{\link{gene_counts_ggplot}()},
\code{\link{individual_plots}()},
\code{\link{mean_sd_plot}()},
\code{\link{pre_process}()},
\code{\link{rRNA_counts_ggplot}()},
\code{\link{total_counts_ggplot}()}
}
\concept{plots}
\concept{preprocess functions}
