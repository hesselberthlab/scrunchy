#' scrunchy: single-cell analysis of functional heterogeneity
#'
#' scrunchy provides tools to analyze and visualize data from
#' single-cell assays for functional heterogeneity.
#'
#' To learn more about scrunchy, start with the vignette:
#' `browseVignettes(package = "scrunchy")`
#'
#' @author Jay Hesselberth <jay.hesselberth@@gmail.com>
#' @author Kent Riemondy <kent.riemondy@@gmail.com>
#' @author Amanda Richer <mandy.l.richer@@gmail.com>
#'
#' @docType package
#' @name scrunchy
#'
#' @seealso Report bugs at [https://github.com/hesselberthlab/scrunchy/issues]
#'
#' @import Matrix
#' @import stringr
#' @import readr
#' @import cowplot
#' @import rlang
#' @import dplyr
#' @import tibble
#' @import fs
#' @import purrr
#' @import ggplot2
#' @import umap
#' @importFrom stats kmeans
#' @importFrom magrittr %>%
#' @importFrom irlba prcomp_irlba
#' @import Rtsne
#' @import SingleCellExperiment
#' @import MultiAssayExperiment
#' @import Seurat
"_PACKAGE"
