library(scrunchy)
library(usethis)

hc_features <- calc_var_features(fsce_small, expt = "haircut", n = 10)
mrna_features <- c(
 "IL7R",
 "CD14",
 "LYZ",
 "MS4A1",
 "CD8A",
 "FCGR3A",
 "MS4A7",
 "GNLY",
 "NKG7",
 "FCER1A",
 "CST3",
 "PPBP"
)

features <- c(hc_features, mrna_features)

fsce_tidy <- purrr::reduce(
  list(
    tidy_dims(fsce_small) %>%
      select(cell_id, starts_with("UMAP"), -experiment),
    tidy_coldata(fsce_small),
    tidy_logcounts(fsce_small[features, , ]) %>%
      select(-experiment)
  ),
  left_join,
  by = "cell_id"
)

usethis::use_data(fsce_tidy, compress = "xz", overwrite = TRUE)
