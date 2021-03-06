% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expt-haircut.R
\name{calc_hairpin_coverage}
\alias{calc_hairpin_coverage}
\title{Calculate coverage across hairpins}
\usage{
calc_hairpin_coverage(fsce, cell_ids = NULL, activities = NULL,
  expt = "haircut", meta = "rnaseq")
}
\arguments{
\item{fsce}{\code{\link{FunctionalSingleCellExperiment}}}

\item{cell_ids}{cell ids to include in coverage calculation or category from colData in rnaseq expt}

\item{activities}{vector of activites to calculate coverage over. default is all present.}

\item{expt}{Data to use for calculating variable features
(default is \code{rnaseq}). Must be present in \code{names(fsce)}.}

\item{meta}{Data to use to find cell_ids categories default is \code{rnaseq}}
}
\description{
Calculate coverage across hairpins
}
\examples{
calc_hairpin_coverage(fsce_small)
}
