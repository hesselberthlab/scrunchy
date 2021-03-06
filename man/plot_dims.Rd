% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_dims}
\alias{plot_dims}
\title{Scatter plot of cells in a two-dimensional embedding.}
\usage{
plot_dims(df, x, y, color = "cell_id", size = 0.1, alpha = 1,
  palette = "cloupe", labels = NULL, label_legend = TRUE,
  label_groups = FALSE)
}
\arguments{
\item{df}{plot data}

\item{x}{variable for x-axis}

\item{y}{variable for y-axis}

\item{color}{variable for point colors (default is black)}

\item{size}{size for \code{\link{geom_point}}}

\item{alpha}{alpha for \code{\link{geom_point}}}

\item{palette}{palette for continuous colors. One of cloupe (the default),
brewer, viridis.}

\item{labels}{labels for groups}

\item{label_legend}{add labels to legend}

\item{label_groups}{add labels to points}
}
\description{
This is the base plot for superimposing annotations like cell types,
cluster assignments, and measured activities.
}
\details{
Embeddings can be calculated by \code{\link[=calc_umap]{calc_umap()}} and \code{\link[=calc_tsne]{calc_tsne()}}, and
retrieved with \code{\link[=tidy_dims]{tidy_dims()}}.
}
\examples{
plot_dims(fsce_tidy, UMAP1, UMAP2, size = 1)

plot_dims(fsce_tidy, UMAP1, UMAP2, IL7R, size = 1)

plot_dims(fsce_tidy, UMAP1, UMAP2, Uracil_45, size = 1)

plot_dims(fsce_tidy, UMAP1, UMAP2, k_cluster, size = 1)

plot_dims(fsce_tidy, UMAP1, UMAP2, k_cluster, labels = LETTERS[1:6])

plot_dims(fsce_tidy, UMAP1, UMAP2, k_cluster,
          labels = LETTERS[1:6], label_groups = TRUE)

}
\seealso{
Other plot functions: \code{\link{plot_activity}},
  \code{\link{plot_hairpin_coverage}}
}
\concept{plot functions}
