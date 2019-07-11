library(tidyverse)
library(usethis)
library(vroom)

# load 3M-february-2018.txt.gz barcodes from:
# cellranger/3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/translation

# info on the file at:
# https://github.com/10XGenomics/cellranger/commit/79d5089d1a7842a95f6ed117b0249beaa3030f71

barcode_map_10x_v3 <- vroom::vroom(
  "data-raw/3M-february-2018.txt.gz",
  col_names = c("src", "dst")
)

usethis::use_data(barcode_map_10x_v3, compress = "xz", overwrite = TRUE)
