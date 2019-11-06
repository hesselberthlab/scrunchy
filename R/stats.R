#' Compare activity differences between groups.
#'
#' Data is assumed to be tidy with a single grouping variable; all other
#' variables are treated as activities to compare. Descriptive columns like
#' `cell_id` should not be included.
#'
#' Applies [`stats::wilcox.test()`] to unique pairs of groups for each measured
#' variable. Adjusted p-values are calculated using [`qvalue::qvalue()`]. A
#' ratio of activity in each group is reported; `Inf` values indicate that
#' no activities were measured for the second group.
#'
#' @param tbl data from a `SingleCellExperiment`
#' @param group variable for generating combinations
#' @param complete If `TRUE`, generate complete group combinations (useful for
#'   e.g., matrix visulatization of p-values). Default is `FALSE`, generating
#'   unique groups combinations.
#' @param ... additional parameters to pass to [`calc_qvalues()`].
#'
#' @examples
#' x <- fsce_tidy[c("k_cluster", "Uracil_45", "riboG_44")]
#' stat_activity_grouped(x, group = k_cluster)
#'
#' @importFrom broom tidy
#' @importFrom stats wilcox.test
#'
#' @return tibble, sorted by `q.value`
#' - `activity`
#' - `group`
#' - `group1`
#' - `ratio` (ratio of the median signals in `group` over `group1`)
#' - `p.value`
#' - `q.value`
#'
#' @family statistical tests
#'
#' @references \doi{10.1038/nmeth.4612}
#'
#' @export
stat_activity_grouped <- function(tbl, group, complete = FALSE, ...) {
  group <- enquo(group)

  activity <- sym("activity")
  value <- sym("value")
  data <- sym("data")
  q.value <- sym("q.value")

  ## gather, group, and flatten to vectors
  x <- gather(tbl, !!activity, !!value, -!!group)
  x <- nest(group_by(x, !!group, !!activity))

  x <- mutate(x, !!data := flatten(!!data))

  ## split groups by activity
  groups <- x$activity
  splits <- split(x, groups)

  ## generate combinations of groups
  crossed <- purrr::map(splits, cross_groups, complete)

  res <- purrr::map_dfr(crossed, group_stat, group, .id = "activity")

  res <- calc_qvalues(res, id = "p.value", ...)

  arrange(res, !!q.value)
}

#' Analysis of variance of activities across groups
#'
#' @param tbl tidied data from a `SingleCellExperiment`
#' @param group variable for generating combinations
#'
#' @examples
#' x <- fsce_tidy[c("k_cluster", "Uracil_45", "riboG_44")]
#' x$k_cluster <- as.factor(x$k_cluster)
#'
#' ## default is list of aov models
#' stat_anova_grouped(x, k_cluster)
#'
#' @export
stat_anova_grouped <- function(tbl, group = NULL) {
  group <- enquo(group)

  activity <- sym("activity")
  value <- sym("value")

  if (is.null(group)) {
    stop("must specify a `group` for the ANOVA", call. = FALSE)
  }

  tbl <- gather(tbl, !!activity, !!value, -!!group)
  groups <- tbl$activity
  tbl_split <- split(tbl, groups)

  purrr::map(
    tbl_split,
    anova_fun,
    stats::as.formula(paste("value ~", quo_name(group)))
  )
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
#'
#' ## first result
#' stat_anova_tukey(res, k_cluster)[[1]]
#'
#' @export
stat_anova_tukey <- function(tbl, group = NULL) {
  group <- enquo(group)

  if (is.null(group)) {
    stop("must specify a `group` for tukey contrasts", call. = FALSE)
  }

  purrr::map(tbl, tukey_fun, group)
}

#' Tidy grouped statistics
#'
#' @param x list of named statistical results, e.g. from [`stat_anova_grouped`].
#' @param id name of id variable in the output (default: "activity")
#'
#' @examples
#' x <- fsce_tidy[c("k_cluster", "Uracil_45")]
#' x$k_cluster <- as.factor(x$k_cluster)
#'
#' res <- stat_anova_grouped(x, k_cluster)
#' tidy_stats_grouped(res)
#'
#' @export
tidy_stats_grouped <- function(x, id = "activity") {
  purrr::map_dfr(x, broom::tidy, .id = id)
}

# Utilities ---------------------------------------------------------

#' Calculate q-values
#'
#' @param x data frame
#' @param id name of variable with p-values
#' @param bh If `TRUE`, use the Benjamini-Hochberg adjustment
#'
#' @importFrom qvalue qvalue
#'
#' @export
calc_qvalues <- function(x, id = "p.value", bh = TRUE) {
  pvs <- x[[id]]

  if (bh) {
    ## pi0 = 1 uses the BH procecure
    qvs <- qvalue::qvalue(pvs, pi0 = 1)
  } else {
    qvs <- qvalue::qvalue(pvs)
  }

  mutate(x, q.value = qvs$qvalues)
}

#' @noRd
#' @importFrom stats aov
anova_fun <- function(x, fmla) {
  stats::aov(fmla, data = x)
}

#' @noRd
#' @importFrom multcomp glht mcp
tukey_fun <- function(x, group) {
  multcomp::glht(
    x, linfct = do.call(
      multcomp::mcp, rlang::list2(!!group := "Tukey")
    )
  )
}

cross_groups <- function(x, complete) {
  group <- sym("group")
  group1 <- sym("group1")

  y <- x

  ## set standardized names for crossed data
  names(x) <- c("group", "activity", "data")
  names(y) <- c("group1", "activity1", "data1")

  ## need to have unique column names for second tibble
  xx <- crossing(x, y)

  if (complete) {
    return(xx)
  }

  xx <- xx[unique_inds(xx, group, group1), ]

  filter(xx, !!group != !!group1)
}

group_stat <- function(x, group) {
  data <- sym("data")
  data1 <- sym("data1")

  res <- mutate(x,
    stat = purrr::map2(!!data, !!data1, tidy_wilcoxon),
    p.value = purrr::map_dbl(stat, "p.value"),
    ratio = purrr::map2_dbl(!!data, !!data1, fold_change)
  )

  cols <- syms(c("group", "group1", "ratio", "p.value"))
  select(res, !!!cols)
}

fold_change <- function(x, y) {
  median(x) / median(y)
}

tidy_wilcoxon <- function(x, y) {
  broom::tidy(suppressWarnings(stats::wilcox.test(x, y)))
}

unique_inds <- function(x, ...) {
  cols <- list2(...)
  inds <- select(x, !!!cols)

  # https://stat.ethz.ch/pipermail/r-help/2011-July/282836.html
  !duplicated(t(apply(inds, 1, sort)))
}
