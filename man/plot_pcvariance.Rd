% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_pcvariance}
\alias{plot_pcvariance}
\title{Plot PCA variance}
\usage{
plot_pcvariance(fsce, n_dims = NULL, expt = "rnaseq")
}
\arguments{
\item{fsce}{An object of class \code{\link{FunctionalSingleCellExperiment}}.}

\item{n_dims}{specify the number of dimensions from "dr" to use for
clustering, defaults to all dimensions}

\item{expt}{Data to use for calculating variable features
(default is \code{rnaseq}). Must be present in \code{names(fsce)}.}
}
\description{
Plots proportion of variance explained by computed prinicpal components
}
\examples{
plot_pcvariance(fsce_small)

}
\seealso{
Other plot fuctions: \code{\link{add_stats}},
  \code{\link{plot_heatmap}}
}
\concept{plot fuctions}
