# Plots -------------------------------------------------------------

#' Scatter plot of cells in a two-dimensional embedding.
#'
#' This is the base plot for superimposing annotations like cell types,
#' cluster assignments, and measured activities.
#'
#' Embeddings can be calculated by [`calc_umap()`] and [`calc_tsne()`], and
#' retrieved with [`tidy_dims()`].
#'
#' @param df plot data
#' @param x variable for x-axis
#' @param y variable for y-axis
#' @param color variable for point colors (default is black)
#' @param size size for [`geom_point`]
#' @param alpha alpha for [`geom_point`]
#' @param palette palette for continuous colors. One of cloupe (the default),
#'   brewer, viridis.
#' @param title plot title
#'
#' @examples
#' plot_dims(fsce_tidy, UMAP1, UMAP2, size = 1)
#'
#' plot_dims(fsce_tidy, UMAP1, UMAP2, Uracil_45, size = 1)
#'
#' plot_dims(fsce_tidy, UMAP1, UMAP2, k_cluster, size = 1)
#'
#' @export
plot_dims <- function(df, x, y, color = "cell_id",
                      size = 0.1, alpha = 1,
                      palette = "cloupe", title = NULL) {
  x <- enquo(x)
  y <- enquo(y)
  color <- enquo(color)

  p <- ggplot(df, aes(x = !!x, y = !!y)) +
    geom_point(
      aes(color = !!color),
      size = size, alpha = alpha
    )

  ## theme default
  p <- p + cowplot::theme_cowplot()

  llim <- legend_limits(df, color)

  if (is_discrete(pull(df, !!color))) {

    ## legend aesthetics
    p <- p + guides(
      colour = guide_legend(
        override.aes = list(size = 4)
      )
    ) + theme(legend.title = element_blank())

    ## color aesthetics
    p <- p + scale_color_manual(
      values = discrete_palette_default
    )

  } else {

    if (palette == "cloupe") {
      p <- p + scale_color_gradientn(
        colors = loupe_palette,
        limits = llim
      )
    } else if (palette == "viridis") {
      p <- p + scale_color_viridis_c(
        option = "inferno",
        limits = llim,
        direction = -1
      )
    } else if (palette == "brewer") {
      p <- p + scale_color_distiller(
        palette = "Reds",
        limits = llim,
        direction = 1
      )
    }

  } # discrete?

  if (!is.null(title)) {
    p <- p + labs(title = title)
  }

  p
}

#' Plot activities per cluster
#'
#' Generates a beeswarm plot of activity across specified groups
#'
#' @param data data to plot
#' @param activity activity variable
#' @param group grouping variable
#' @param legend display legend
#'
#' @examples
#' cowplot::plot_grid(
#'   plotlist = list(
#'     plot_dims(fsce_tidy, UMAP1, UMAP2, k_cluster, title = "mRNA expression"),
#'     plot_activity(fsce_tidy, Uracil_45, k_cluster, legend = FALSE),
#'     plot_activity(fsce_tidy, riboG_44, k_cluster, legend = FALSE)
#'   )
#' )
#'
#' @export
plot_activity <- function(data, activity, group = NULL, legend = TRUE) {

  activity <- enquo(activity)
  group <- enquo(group)

  p <- ggplot(data, aes(x = !! activity, y = !! group, color = !! group)) +
    ggbeeswarm::geom_quasirandom(size = 0.5, groupOnX = FALSE) +
    scale_color_OkabeIto() +
    cowplot::theme_cowplot() +
    labs(
      x = "Activity",
      y = "Group",
      title = glue(
        "Activity: {activity}",
        activity = rlang::quo_text(activity)
      )
    )

  if (!legend) {
    p <- p + theme(legend.position = "none")
  }

  p
}

#' Heatmap of signals
#'
#' Plots `logcounts` or `counts` from an experiment for specified rows.
#'
#' @import ComplexHeatmap
#'
#' @param mtx Matrix of `logcounts` or `counts`
#' @param rows names of rows to select for heatmap
#' @param ... params for [`ComplexHeatmap::Heatmap`]
#'
#' @examples
#' rows <- paste("Uracil", 10:50, sep = "_")
#' mtx <- logcounts(fsce_small[["haircut"]])
#'
#' plot_heatmap(mtx, rows, name = "Uracil")
#'
#' @export
plot_heatmap <- function(mtx, rows = NULL, ...) {
  if (!is.null(rows)) {
    mtx <- mtx[rows, ]
  }

  # traspose and strip rownames (cell ids)
  mtx <- t(mtx)
  rownames(mtx) <- NULL

  ComplexHeatmap::Heatmap(mtx, cluster_columns = FALSE, ...)
}

# Palettes ----------------------------------------------------------

loupe_palette <- rev(scales::brewer_pal(palette = "RdGy")(11)[c(1:5, 7)])

#' @noRd
#' @include reexport-colorblindr.R
discrete_palette_default <- c(
  palette_OkabeIto_black,
  scales::brewer_pal(palette = "Paired")(12),
  scales::brewer_pal(palette = "Set1")(9),
  scales::brewer_pal(palette = "Set2")(8),
  scales::brewer_pal(palette = "Dark2")(8)
)

# Utilities ---------------------------------------------------------

legend_limits <- function(x, var) {
  if (is_discrete(pull(x, !!var))) {
    c(NA, NA)
  }
  c(0, max(pull(x, !!var)))
}

is_discrete <- function(x) {
  is_character(x) | is_logical(x) | is.factor(x)
}
