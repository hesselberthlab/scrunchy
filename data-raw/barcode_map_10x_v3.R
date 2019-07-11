library(tidyverse)
library(usethis)

# load 3M-february-2018.txt.gz barcodes from:
# cellranger/3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/translation

barcode_map_10x_v3 <- readr::read_tsv(
  "data-raw/3M-february-2018.txt.gz",
  col_names = c("id1", "id2")
)

usethis::use_data(barcode_map_10x_v3, compress = "xz")
