#' scrunchy: single-cell reconstruction of functional heterogeneity
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
#' @seealso Report bugs at <https://github.com/hesselberthlab/scrunchy/issues>
#'
#' @importFrom methods as
#' @importFrom stats median var kmeans as.formula
#'
#' @import Matrix
#' @import readr
#' @import dplyr
#' @import tibble
#' @import fs
#' @import purrr
#' @import ggplot2
#' @import umap
#'
#' @importFrom tidyr gather nest crossing
#'
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom glue glue
#' @importFrom cowplot theme_cowplot plot_grid
#' @importFrom rlang parse_quosure is_character is_logical quo_text
#' @importFrom rlang list2 .data :=
#' @importFrom magrittr %>%
#' @importFrom irlba prcomp_irlba
#' @importFrom Rtsne Rtsne
#'
#' @import SingleCellExperiment
#' @import MultiAssayExperiment
#' @importFrom S4Vectors SimpleList
#'
#' @useDynLib scrunchy, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
