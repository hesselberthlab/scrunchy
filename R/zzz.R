
#' @noRd
#' @importFrom reticulate import
.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  leidenalg_py <<- reticulate::import("leidenalg", delay_load = TRUE)
  igraph_py <<- reticulate::import("igraph", delay_load = TRUE)
}
