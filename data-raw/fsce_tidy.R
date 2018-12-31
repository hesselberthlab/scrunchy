library(scrunchy)
library(usethis)

features <- calc_var_features(fsce_small, expt = "haircut", n = 10)

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
