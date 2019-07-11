library(tidyverse)
library(usethis)
library(vroom)

# load 3M-february-2018.txt.gz barcodes from:
# cellranger/3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/translation

barcode_map_10x_v3 <- vroom::vroom(
  "data-raw/3M-february-2018.txt.gz",
  col_names = c("src", "dst")
)

usethis::use_data(barcode_map_10x_v3, compress = "xz", overwrite = TRUE)

