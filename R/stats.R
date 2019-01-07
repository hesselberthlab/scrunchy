#' Compare activity differences between groups.
#'
#' Data is assumed to be tidy with a single grouping variable; all other
#' variables are treated as activities to compare. Descriptive columns like
#' `cell_id` should not be included.
#'
#' Applies [`stats::wilcox.test()`] to unique pairs of groups for each measured
#' variable.
#'
#' @param df tidied vdata from a `SingleCellExperiment`
#' @param group variable for generating combinations
#' @param complete If `TRUE`, generate complete group combinations (useful for
#'   e.g. matrix visulatization if p-values). Default is `FALSE`, generating
#'   unique groups combinations.
#'
#' @examples
#' x <- fsce_tidy[c("k_cluster", "Uracil_45", "riboG_44")]
#' stat_activity_grouped(x, group = k_cluster)
#'
#' @importFrom broom tidy
#' @importFrom stats wilcox.test
#'
#' @return tibble sorted by `p.value` of the test.
#'
#' @family statistical tests
#'
#' @references \doi{10.1038/nmeth.4612}
#'
#' @export
stat_activity_grouped <- function(df, group, complete = FALSE) {

  group <- enquo(group)

  ## gather, group, and flatten to vectors
  x <- gather(df, activity, value, -!!group)
  x <- nest(group_by(x, !!group, activity))
  x <- mutate(x, data = flatten(data))

  ## split groups by activity
  groups <- x$activity
  splits <- split(x, groups)

  ## generate ucombinations of groups
  crossed <- purrr::map(splits, cross_groups, complete)

  res <- purrr::map_dfr(crossed, group_stat, group, .id = "activity")

  arrange(res, p.value)
}

# Utilities ---------------------------------------------------------

cross_groups <- function(x, complete) {
  ## set standardized names for crossed data
  names(x) <- c("group", "activity", "data")

  xx <- crossing(x, x)

  if (complete) {
    return(xx)
  }

  xx <- xx[unique_inds(xx, group, group1), ]

  filter(xx, group != group1)
}

group_stat <- function(x, group) {
  res <- mutate(x,
    stat = purrr::map2(data, data1, tidy_wilcoxon),
    p.value = purrr::map_dbl(stat, "p.value")
  )

  select(res, group, group1, p.value)
}

tidy_wilcoxon <- function(x, y) {
  broom::tidy(suppressWarnings(stats::wilcox.test(x, y)))
}

unique_inds <- function(x, ...) {
  inds <- select(x, ...)

  # https://stat.ethz.ch/pipermail/r-help/2011-July/282836.html
  !duplicated(t(apply(inds, 1, sort)))
}
