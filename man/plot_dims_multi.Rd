% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_dims_multi}
\alias{plot_dims_multi}
\title{Plot multiple 2D plots in a grid}
\usage{
plot_dims_multi(df, features, ...)
}
\arguments{
\item{df}{plot data}

\item{features}{list of features}

\item{...}{params to pass to \code{\link[=plot_dims]{plot_dims()}}}
}
\description{
Plot multiple 2D plots in a grid
}
\examples{
plot_dims_multi(
  fsce_tidy,
  features = c("k_cluster", "Uracil_45", "IL7R", "GNLY"),
  x = UMAP1, y = UMAP2, size = 0.5
)

}
