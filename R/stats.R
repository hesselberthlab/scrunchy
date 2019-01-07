#' Compare activity differences between groups.
#'
#' Data is assumed to be tidy with a single grouping variable; all other
#' variables are treated as activities to compare. Descriptive columns like
#' `cell_id` should not be included.
#'
#' Applies [`stats::wilcox.test()`] to unique pairs of groups for each measured
#' variable.
#'
#' @param tbl data from a `SingleCellExperiment`
#' @param group variable for generating combinations
#' @param complete If `TRUE`, generate complete group combinations (useful for
#'   e.g., matrix visulatization of p-values). Default is `FALSE`, generating
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
stat_activity_grouped <- function(rbl, group, complete = FALSE) {
  group <- enquo(group)

  ## gather, group, and flatten to vectors
  x <- gather(df, activity, value, -!!group)
  x <- nest(group_by(x, !!group, activity))
  x <- mutate(x, data = flatten(data))

  ## split groups by activity
  groups <- x$activity
  splits <- split(x, groups)

  ## generate combinations of groups
  crossed <- purrr::map(splits, cross_groups, complete)

  res <- purrr::map_dfr(crossed, group_stat, group, .id = "activity")

  arrange(res, p.value)
}

#' Analysis of variance of activities across groups
#'
#' @param df tidied data from a `SingleCellExperiment`
#' @param group variable for generating combinations
#' @param tidy tidy the results
#'
#' @examples
#' x <- fsce_tidy[c("k_cluster", "Uracil_45", "riboG_44")]
#' x$k_cluster <- as.factor(x$k_cluster)
#'
#' # default is list of aov models
#' stat_anova_grouped(x, k_cluster)
#'
#' ## tidy results
#' stat_anova_grouped(x, k_cluster, tidy = TRUE)
#'
#' @export
stat_anova_grouped <- function(tbl, group, tidy = FALSE) {
  group <- enquo(group)

  tbl <- gather(tbl, activity, value, -!!group)
  groups <- tbl$activity
  tbl_split <- split(tbl, groups)

  res <- purrr::map(
    tbl_split,
    anova_fun,
    as.formula(paste("value ~", quo_name(group)))
  )

  if (tidy) {
    return(purrr::map_dfr(res, broom::tidy, .id = "activity"))
  }

  res
}

#' Post-hoc analysis of ANOVA results
#'
#' @inheritParams stat_anova_grouped
#'
#' @examples
#' x <- fsce_tidy[c("k_cluster", "Uracil_45", "riboG_44")]
#' x$k_cluster <- as.factor(x$k_cluster)
#'
#' res <- stat_anova_grouped(x, k_cluster)
#' stat_anova_tukey(res, tidy = TRUE)
#'
#' @export
stat_anova_tukey <- function(tbl, group, tidy = FALSE) {
  res <- purrr::map(tbl, tukey_fun, group)

  if (tidy) {
    return(purrr::map_dfr(res, broom::tidy, .id = "activity"))
  }

  res
}

# Utilities ---------------------------------------------------------

#' @noRd
#' @importFrom stats aov
anova_fun <- function(x, fmla) {
  stats::aov(fmla, data = x)
}

#' @noRd
#' @importFrom multcomp glht mcp
tukey_fun <- function(x, group) {
  multcomp::glht(
    x, linfct = multcomp::mcp(k_cluster = "Tukey")
  )
}

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
