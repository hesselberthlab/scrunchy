#' Calculate coverage across hairpins
#'
#' @param fsce [`FunctionalSingleCellExperiment`]
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param cell_ids cell ids to include in coverage calculation
#'
#' @examples
#' calc_hairpin_coverage(fsce_small)
#' @export
calc_hairpin_coverage <- function(fsce,
                                  expt = "haircut",
                                  cell_ids = NULL) {
  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce"), call. = FALSE)
  }

  if (is.null(cell_ids)) {
    cell_ids <- colnames(counts(fsce[[expt]]))
  }

  metadata <- colData(fsce[[expt]])[cell_ids, , drop = FALSE]
  counts <- counts(fsce[[expt]])

  cell_ids <- split(rownames(metadata), metadata$cell_id)

  ## compute per sample per adduct position counts
  res <- purrr::map_dfr(
    cell_ids,
    function(ids) {
      hairpin_info <- as.data.frame(rowData(fsce[[expt]]))
      hairpin_info$count <- rowSums(counts[, ids, drop = FALSE])
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
#' @param color variable to use for coloring lines (defaults to "hairpin")
#'
#' @examples
#' plot_hairpin_coverage(fsce_small) + ggplot2::facet_wrap(~hairpin)
#' @export
plot_hairpin_coverage <- function(fsce,
                                  expt = "haircut",
                                  cell_ids = NULL,
                                  color = "hairpin") {
  res <- calc_hairpin_coverage(fsce, expt, cell_ids)

  color <- enquo(color)

  ggplot(res, aes(x = position, y = count)) +
    geom_line(aes(color = !!color), alpha = 0.7, size = 0.8) +
    cowplot::theme_cowplot() +
    scale_color_brewer(palette = "Set1") +
    labs(
      x = "Position",
      y = "Counts"
    )
}
