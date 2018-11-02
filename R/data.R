#' Generate Seurat object from barnyard data
#'
#' Load data from barnyard experiment peformed with
#' UNG and RNASEH2C Hap1 cells and A:U and rG hairpin
#' substrates.
#'
#' @importFrom magrittr %>%
#' @import Seurat
#'
#' @export
barnyard_data <- function() {

  mrna_tab <- load_csv(scrunchy_data("mrna.csv.gz"))
  haircut_tab <- load_csv(scrunchy_data("haircut.csv.gz"))

  sce <- CreateSeuratObject(raw.data = mrna_tab) %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE, y.cutoff = 0.5) %>%
    ScaleData(display.progress = FALSE) %>%
    RunPCA(pcs.print = 0) %>%
    FindClusters(dims.use = 1:6, print.output = FALSE) %>%
    RunTSNE(dims.use = 1:6)

  sce <- SetAssayData(
    sce, assay.type = 'CITE',
    slot = "raw.data",
    new.data = haircut_tab
  ) %>%
    NormalizeData(assay.type = "CITE", normalization.method = "genesCLR") %>%
    ScaleData(assay.type = "CITE", display.progress = FALSE)

  sce
}

scrunchy_data <- function(x) {
  system.file("extdata", x, package = "scrunchy", mustWork = TRUE)
}

load_csv <- function(x) {
  read.csv(x, sep = ",", header = TRUE, row.names = 1)
}
