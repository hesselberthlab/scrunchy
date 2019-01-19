library(scrunchy)
library(usethis)

seed <- 42

fsce_small <- create_fsce(
  list(
    rnaseq = create_sce_rnaseq(scrunchy_data("mrna")),
    haircut = create_sce_haircut(scrunchy_data("haircut"))
  )
)

# pre-calculate PCA, UMAP, t-SNE, and k-means clusters
var_genes <- calc_var_features(fsce_small, "rnaseq", n = 2000)

fsce_small <- calc_pca(fsce_small, n_pcs = 20, genes = var_genes, seed = seed)
fsce_small <- calc_umap(fsce_small, n_dims = 6, seed = seed)
fsce_small <- calc_tsne(fsce_small, n_dims = 6, seed = seed)

fsce_small <- cluster_kmeans(fsce_small, k = 6, seed = seed)
fsce_small <- cluster_leiden(fsce_small, seed = seed)

fsce_small <- calc_cell_cycle(fsce_small)

usethis::use_data(fsce_small, compress = "xz", overwrite = TRUE)
