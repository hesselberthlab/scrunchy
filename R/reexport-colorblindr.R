#' Color palette proposed by Okabe and Ito
#'
#' These functions are copied from [colorblindr](https://github.com/clauswilke/colorblindr).
#'
#' Two color palettes taken from the article "Color Universal Design" by Okabe and Ito, http://jfly.iam.u-tokyo.ac.jp/color/.
#' The variant `palette_OkabeIto` contains a gray color, while `palette_OkabeIto_black` contains black instead.
#' @export
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

#' @rdname palette_OkabeIto
#' @export
palette_OkabeIto_black <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

#' @rdname scale_OkabeIto
#' @param aesthetics aesthetic mappings
#' @export
#' @usage NULL
scale_colour_OkabeIto <- function(aesthetics = "colour", ...) {
  scale_OkabeIto(aesthetics, ...)
}

#' @rdname scale_OkabeIto
#' @export
scale_color_OkabeIto <- scale_colour_OkabeIto

#' @rdname scale_OkabeIto
#' @inheritParams scale_colour_OkabeIto
#' @export
scale_fill_OkabeIto <- function(aesthetics = "fill", ...) {
  scale_OkabeIto(aesthetics, ...)
}

#' Okabe-Ito color scale
#'
#' This is a color-blind friendly, qualitative scale with eight different
#' colors. See [palette_OkabeIto] for details.
#'
#' This code is copied and modified from colorblindr to remove the `darken`
#' param, which is only available in unreleased colorspace 0.4.1.
#'
#' @param use_black If `TRUE`, scale includes black, otherwise includes gray.
#' @param order Numeric vector listing the order in which the colors should be used. Default is 1:8.
#' @param alpha Alpha transparency level of the color. Default is no transparency.
#' @param ... common discrete scale parameters: `name`, `breaks`, `labels`, `na.value`, `limits`, `guide`, and `aesthetics`.
#'  See [discrete_scale] for more details.
#'
#' @examples
#' library(ggplot2)
#' ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
#'   geom_point() + scale_color_OkabeIto()
#' ggplot(iris, aes(Sepal.Length, fill = Species)) +
#'   geom_density(alpha = 0.7) + scale_fill_OkabeIto(order = c(1, 3, 5))
#'
#' @export
scale_OkabeIto <- function(aesthetics, use_black = FALSE, order = 1:8, alpha = NA, ...) {
  if (use_black) {
    values <- palette_OkabeIto_black[order]
  } else {
    values <- palette_OkabeIto[order]
  }

  n <- length(values)
  alpha <- rep_len(alpha, n)

  ai <- !is.na(alpha)
  if (sum(ai) > 0) { # at least one color needs alpha
    values[ai] <- scales::alpha(values[ai], alpha[ai])
  }

  pal <- function(n) {
    if (n > length(values)) {
      warning("Insufficient values in manual scale. ", n, " needed but only ",
        length(values), " provided.",
        call. = FALSE
      )
    }
    values
  }
  ggplot2::discrete_scale(aesthetics, "manual", pal, ...)
}
