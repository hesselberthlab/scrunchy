# process mRNA and repair product data through Seurat pipeline

library(Seurat)
library(tidyverse)

mrna_tab <- read.csv("data-raw/mrna.csv.gz", sep = ',', header = TRUE, row.names = 1)
haircut_tab <- read.csv("data-raw/haircut.csv.gz", sep = ',', header = TRUE, row.names = 1)

sce <- CreateSeuratObject(raw.data = mrna_tab) %>%
  NormalizeData() %>%
  FindVariableGenes(do.plot = FALSE, y.cutoff = 0.5) %>%
  ScaleData(display.progress = FALSE) %>%
  RunPCA(pcs.print = 0) %>%
  FindClusters(dims.use = 1:6, print.output = FALSE) %>%
  RunTSNE(dims.use = 1:6)

TSNEPlot(sce, do.label = TRUE, pt.size = 0.5)

sce <- SetAssayData(
  sce, assay.type = 'CITE',
  slot = "raw.data",
  new.data = haircut_tab
) %>%
  NormalizeData(assay.type = "CITE", normalization.method = "genesCLR") %>%
  ScaleData(assay.type = "CITE", display.progress = FALSE)

FeaturePlot(
  sce,
  features.plot = c("Uracil_45", "riboG_44","ENSG00000076248", "ENSG00000172922"),
  cols.use = c("lightgrey", "blue")
)
