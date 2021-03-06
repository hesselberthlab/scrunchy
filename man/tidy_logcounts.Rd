% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidiers.R
\name{tidy_logcounts}
\alias{tidy_logcounts}
\title{Tidy logcounts data}
\usage{
tidy_logcounts(fsce)
}
\arguments{
\item{fsce}{An object of class \code{\link{FunctionalSingleCellExperiment}}.}
}
\description{
Tidy logcounts data
}
\section{Tidying}{

This is a data tidier for a \code{\link{FunctionalSingleCellExperiment}}. Returns data
from \code{\link{SingleCellExperiment}} in a tidy format, where variables are columns
and observations are rows.

If the \code{\link{FunctionalSingleCellExperiment}} contains more than one
\code{\link{SingleCellExperiment}}, data from each \code{sce} are joined using the \code{cell_id}
variable name, and a new \code{experiment} column contains the name of the source
\code{sce}.
}

\examples{
x <- fsce_small[c("Uracil_45"), , "haircut"]
tidy_logcounts(x)

}
\seealso{
Other tidiers: \code{\link{tidy_all}},
  \code{\link{tidy_coldata}}, \code{\link{tidy_counts}},
  \code{\link{tidy_dims}}
}
\concept{tidiers}
