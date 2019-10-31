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
#' @param labels labels for groups
#' @param label_legend add labels to legend
#' @param label_groups add labels to points
#'
#' @examples
#' plot_dims(fsce_tidy, UMAP1, UMAP2, size = 1)
#'
#' plot_dims(fsce_tidy, UMAP1, UMAP2, IL7R, size = 1)
#'
#' plot_dims(fsce_tidy, UMAP1, UMAP2, Uracil_45, size = 1)
#'
#' plot_dims(fsce_tidy, UMAP1, UMAP2, k_cluster, size = 1)
#'
#' plot_dims(fsce_tidy, UMAP1, UMAP2, k_cluster, labels = LETTERS[1:6])
#'
#' plot_dims(fsce_tidy, UMAP1, UMAP2, k_cluster,
#'           labels = LETTERS[1:6], label_groups = TRUE)
#'
#' @family plot functions
#'
#' @importFrom forcats fct_count
#'
#' @export
plot_dims <- function(df, x, y, color = "cell_id",
                      size = 0.1, alpha = 1,
                      palette = "cloupe",
                      labels = NULL,
                      label_legend = TRUE,
                      label_groups = FALSE) {
  x <- enquo(x)
  y <- enquo(y)
  color <- enquo(color)

  p <- ggplot(df, aes(x = !!x, y = !!y)) +
    geom_point(
      aes(color = !!color),
      size = size, alpha = alpha
    )

  ## theme default
  p <- p + cowplot::theme_minimal_grid(line_size = 0.2)

  if (is_discrete(pull(df, !!color))) {

    ## get labels
    n_col <- n_colors(df, color, labels)

    ## legend aesthetics
    p <- p + guides(
      colour = guide_legend(
        override.aes = list(size = 4)
      )
    ) + theme(legend.title = element_blank())

    ## color aesthetics
    if (label_legend && !is.null(labels)) {
      lbls <- labels
    } else {
      lbls <- n_col$f
    }

    p <- p + scale_color_manual(
      values = discrete_palette_default,
      labels = lbls
    )

  } else {

    llim <- legend_limits(df, color)

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

  if (label_groups) {
    p <- add_group_labels(p, x, y, color, labels)
  }

  p
}

#' Plot multiple 2D plots in a grid
#'
#' @param df plot data
#' @param features list of features
#' @param ... params to pass to [`plot_dims()`]
#'
#' @examples
#' plot_dims_multi(
#'   fsce_tidy,
#'   features = c("k_cluster", "Uracil_45", "IL7R", "GNLY"),
#'   x = UMAP1, y = UMAP2, size = 0.5
#' )
#'
#' @export
plot_dims_multi <- function(df, features, ...) {
  plts <- list()

  for (i in seq_along(features)) {
    feat <- features[i]
    plts[[i]] <- plot_dims(df, color = !!sym(feat), ...)
  }

  cowplot::plot_grid(plotlist = plts)
}

#' Plot activities per cluster
#'
#' Generates a beeswarm plot of activity across specified groups
#'
#' @param data data to plot
#' @param activity activity variable
#' @param group grouping variable
#' @param labels legend labels
#' @param vertical should activity be on the y axis. Default is FALSE
#' @param stats dataframe to use to add p values on plot. Stats can be a dataframe output from [`stat_activity_grouped()`].
#'     Output must include columsn with names group, group1, q.value and/or p.value. If sepcificied, `vertical` is TRUE
#' @param ... params for [`add_stats()`]
#'
#'
#' @examples
#' plot_activity(fsce_tidy, Uracil_45, k_cluster)
#'
#' plot_activity(fsce_tidy, riboG_44, k_cluster, labels = LETTERS[1:6])
#'
#' plot_activity(fsce_tidy, Uracil_45, k_cluster, vertical = TRUE)
#'
#' x <- fsce_tidy[c("k_cluster", "Uracil_45")]
#' stats <- stat_activity_grouped(x, group = k_cluster)
#' stats <- subset(stats, q.value < 0.01)
#' stats[,"y_loc"] = seq(max(x$Uracil_45), max(x$Uracil_45) + 3,  length.out = length(stats$group))
#' plot_activity(fsce_tidy, Uracil_45, k_cluster, stats = stats)
#'
#' @family plot functions
#'
#' @export
plot_activity <- function(data, activity, group = NULL, labels = NULL, vertical = FALSE,
                          stats = NULL, ...) {
  x <- enquo(activity)
  y <- enquo(group)
  group <- enquo(group)
  groupOnX = FALSE
  x_lab = "Activity"
  y_lab = "Group"

  n_col <- n_colors(data, group, labels)

  if(!is.null(stats)){
    vertical <- TRUE
  }

  if(vertical){
    x <- enquo(group)
    y <- enquo(activity)
    groupOnX = TRUE
    x_lab = "Group"
    y_lab = "Activity"
  }


  p <- ggplot(data, aes(x = !!x, y = !!y, color = !!group)) +
    ggbeeswarm::geom_quasirandom(size = 0.5, groupOnX = groupOnX) +
    scale_color_OkabeIto(use_black = TRUE, labels = labels %||% n_col$f) +
    cowplot::theme_cowplot() +
    labs(x = x_lab, y = y_lab)

  if(!is.null(stats)){
    p <- add_stats(p, stats, ...)
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
#' @param columns names of columns to select for heatmap
#' @param ... params for [`ComplexHeatmap::Heatmap`]
#'
#' @examples
#' mtx <- SingleCellExperiment::logcounts(fsce_small[["haircut"]])
#' rows <- paste("Uracil", 1:61, sep = "_")
#'
#' plot_heatmap(mtx, rows, name = "Uracil")
#'
#' @family plot fuctions
#'
#' @export
plot_heatmap <- function(mtx, rows = NULL, columns = NULL, ...) {
  if (!is.null(rows)) {
    mtx <- mtx[rows, ]
  }

  if (!is.null(columns)) {
    mtx <- mtx[ ,columns]
  }

  if (class(mtx) %in% c("dgCMatrix")) {
    mtx <- as.matrix(mtx)
  }

  # traspose and strip rownames (cell ids)
  mtx <- t(mtx)
  rownames(mtx) <- NULL

  ComplexHeatmap::Heatmap(mtx, cluster_columns = FALSE, ...)
}

#' Plot PCA variance
#'
#' Plots proportion of variance explained by computed prinicpal components
#'
#' @param fsce An object of class [`FunctionalSingleCellExperiment`].
#' @param n_dims specify the number of dimensions from "dr" to use for
#'   clustering, defaults to all dimensions
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @examples
#' plot_variance(fsce_small)
#'
#' @family plot fuctions
#'
#' @export
plot_pcvariance <- function(fsce, n_dims = NULL, expt = "rnaseq") {

  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce "), call. = FALSE)
  }

  if (!"PCA" %in% names(reducedDims(fsce[[expt]]))) {
    stop("PCA values not found in expt", call. = FALSE)
  }

  var_df <- pcvariance_tbl(fsce[[expt]])

  if(!is.null(n_dims)){
    var_df <- var_df[1:n_dims, ]
  }

  ggplot(var_df, aes(PCs, `Variance Explained`)) +
    geom_point() +
    cowplot::theme_cowplot()
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

centroids <- function(df, x, y, group = NULL) {
  if (!is.null(group)) {
    df <- group_by(df, !!group)
  }
  summarize(df, x = median(!!x), y = median(!!y))
}

#' @importFrom ggrepel geom_label_repel
add_group_labels <- function(p, x, y, group, labels) {
  cents <- centroids(p$data, x, y, group)
  p + geom_label_repel(data = cents, aes(x = x, y = y, label = labels))
}

n_colors <- function(x, color, labels) {
  n_col <- fct_count(pull(x, !!color))

  if (!is.null(labels) && (length(labels) != nrow(n_col))) {
    stop(glue("`labels` ({nl}) must match factors in `{color}` ({nc})",
              color = rlang::quo_text(color),
              nl = length(labels),
              nc = nrow(n_col)), call. = FALSE)
  }

  n_col
}

#' Add stats to plot activities per cluster
#'
#' Adds stat comparison bars to a beeswarm plot of activity across specified groups
#'
#' @importFrom ggsignif geom_signif
#' @param p plot to add stats to
#' @param df dataframe output from [`stat_activity_grouped()`].
#'     Output must include columnw with names group, group1, q.value and/or p.value.
#' @param val values to add to plot. Default is q.value from  [`stat_activity_grouped()`] output
#' @param y_loc location of comparisson bars on graph. Default it `y_loc` column of `df`. Can also be numeric values.
#' @param xmin start location of comparison bar. Default is `group` column of `df`.
#' @param xmax stop location of comparison bar. Default is `group1` column of `df`.
#'
#' @family plot fuctions
#'
#' @export

add_stats <- function(p, df,
                      val = q.value,
                      y_loc = y_loc,
                      xmin = group,
                      xmax = group1) {

  val <- enquo(val)
  y_loc <- enquo(y_loc)
  xmin <- enquo(xmin)
  xmax = enquo(xmax)

  p + geom_signif(data=df,
                  aes(xmin = !!xmin,
                      xmax = !!xmax,
                      annotations = signif(!!val, digits = 3),
                      y_position = !!y_loc
                  ),
                  color='black',
                  manual = TRUE)
}
