#' Calculate coverage across hairpins
#'
#' @param fsce [`FunctionalSingleCellExperiment`]
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param cell_ids cell ids to include in coverage calculation or category from colData in rnaseq expt
#' @param meta Data to use to find cell_ids caategories default is `rnaseq`
#'
#' @examples
#' calc_hairpin_coverage(fsce_small)
#' @export
calc_hairpin_coverage <- function(fsce,
                                  cell_ids = NULL,
                                  activities = NULL,
                                  expt = "haircut",
                                  meta = "rnaseq"
                                  ) {
  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce"), call. = FALSE)
  }


  if (is.null(cell_ids)) {
    category <- "all_cells"
    cell_ids <- colnames(counts(fsce[[expt]]))
  }

  metadata <- colData(fsce[[meta]])

  if(length(cell_ids) == 1){
    if (cell_ids %in% colnames(metadata)) {
    category <- metadata[,cell_ids]
    cell_ids <- colnames(counts(fsce[[expt]]))
    }
  }

  metadata <- metadata[cell_ids, , drop = FALSE]

  if(is.null(activities)){
    activities <- rownames(fsce[[expt]])
  }

  fsce_to_plot <- fsce[[expt]][activities, , drop = FALSE]
  counts <- counts(fsce_to_plot)

  cell_ids <- split(rownames(metadata), category)

  ## compute per sample per adduct position counts
  res <- purrr::map_dfr(
    cell_ids,
    function(ids) {
      hairpin_info <- as.data.frame(rowData(fsce_to_plot))
      hairpin_info$count <- log1p(rowMeans(counts[, ids, drop = FALSE]))
      hairpin_info
    },
    .id = "cell_id"
  )

  res$position <- as.numeric(res$position)
  as_tibble(res)
}

#' Plot coverage across hairpins
#'
#' @inheritParams calc_hairpin_coverage
#' @param color variable to use for coloring lines (defaults to "cell_id")
#' @param use_points if TRUE add points to plot
#'
#' @examples
#' plot_hairpin_coverage(fsce_small) + ggplot2::facet_wrap(~hairpin)
#'
#' plot_hairpin_coverage(fsce_small, cell_ids = "k_cluster") + ggplot2::facet_wrap(~hairpin)
#'
#' @family plot functions

#' @export
plot_hairpin_coverage <- function(fsce,
                                  cell_ids = NULL,
                                  activities = NULL,
                                  expt = "haircut",
                                  meta = "rnaseq",
                                  color = "cell_id",
                                  use_points = FALSE) {

  res <- calc_hairpin_coverage(fsce, cell_ids, activities, expt, meta)

  color <- enquo(color)
  colorN <- quo_name(color)

  p <- ggplot(res, aes(x = position, y = count)) +
    geom_line(aes_string(color = colorN), alpha = 0.7, size = 0.8) +
    cowplot::theme_cowplot() +
    scale_color_brewer(palette = "Set1") +
    labs(
      x = "Position",
      y = "Counts"
    )

  if (use_points){
    p <- p + geom_point(aes_string(color = colorN))
  }
  p
}
