#' Calculate statistics between pairs of groups
#'
#' Applies [`stats::wilcox.test()`] to all combinations of variables across
#' groups in a tidied data frame.
#'
#' Data is assumed to contain a single grouping variable, and all other
#' variables are treated as activitys to compare. Descriptive columns like
#' `cell_id` should not be included.
#'
#' @param df tidied version of data from a `SingleCellExperiment`
#' @param group variable for generating combinations
#'
#' @examples
#' x <- fsce_tidy[c("k_cluster", "Uracil_45", "riboG_44")]
#' x
#'
#' calc_group_stats(x, group = k_cluster)
#'
#' @importFrom broom tidy
#'
#' @export
calc_group_stats <- function(df, group) {

  group <- enquo(group)

  ## gather, group, and flatten to vectors
  x <- gather(df, activity, value, -!!group)
  x <- nest(group_by(x, !!group, activity))
  x <- mutate(x, data = flatten(data))

  ## split groups by activity
  groups <- x$activity
  splits <- split(x, groups)

  ## generate unique combinations of groups
  crossed <- purrr::map(splits, cross_groups)

  res <- purrr::map_dfr(crossed, group_stat, group, .id = "activity")

  arrange(res, p.value)
}

cross_groups <- function(x, group) {
  ## set standardized names for crossed data
  names(x) <- c("group", "activity", "data")

  xx <- crossing(x, x)
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
  broom::tidy(suppressWarnings(wilcox.test(x, y)))
}

unique_inds <- function(x, ...) {
  inds <- select(x, ...)

  # https://stat.ethz.ch/pipermail/r-help/2011-July/282836.html
  !duplicated(t(apply(inds, 1, sort)))
}
