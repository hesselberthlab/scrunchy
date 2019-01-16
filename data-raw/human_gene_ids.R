library(scrunchy)
library(tidyverse)

tab <- read_tsv(
  "~/devel/scrunchy/inst/extdata/features.tsv.gz",
  col_names = c("id", "symbol", "type")
)

human_gene_ids <- select(tab, id, symbol) %>%
  dplyr::rename(ensembl_id = id, gene_symbol = symbol)

## remove duplicated gene_symbols
dups <- duplicated(human_gene_ids$gene_symbol)
human_gene_ids <- human_gene_ids[!dups, ]

usethis::use_data(human_gene_ids, compress = "xz", overwrite = TRUE)
