library(scrunchy)
library(tidyverse)
library(usethis)

features_hairpin <- c("Normal_45", "Uracil_45",
                      "riboG_44", "Abasic_46")
features_mrna <- c("IL7R", "CD14", "LYZ", "MS4A1",
              "CD8A", "FCGR3A", "MS4A7",
              "GNLY", "NKG7", "FCER1A",
              "CST3", "PPBP")


fsce_tidy <- purrr::reduce(
  list(
    tidy_dims(fsce_small) %>%
      select(cell_id, starts_with("UMAP"), -experiment),
    tidy_coldata(fsce_small),
    tidy_logcounts(fsce_small[features_hairpin, , ]) %>%
      select(-experiment),
    tidy_logcounts(fsce_small[features_mrna, , ]) %>%
      select(-experiment)
  ),
  left_join,
  by = "cell_id"
)

usethis::use_data(fsce_tidy, compress = "xz", overwrite = TRUE)
