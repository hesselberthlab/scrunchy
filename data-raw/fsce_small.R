library(scrunchy)
library(usethis)

fsce_small <- create_fsce(
  list(
    rnaseq = create_sce_rnaseq(scrunchy_data("mrna/")),
    haircut = create_sce_haircut(scrunchy_data("haircut/"))
  )
)

# pre-calculate PCA, UMAP, t-SNE, and k-means clusters
var_genes <- calc_var_features(fsce_small, "rnaseq", n = 5000)

fsce_small <- calc_pca(fsce_small, n_pcs = 20, genes = var_genes)
fsce_small <- calc_umap(fsce_small, n_dims = 6)
fsce_small <- calc_tsne(fsce_small, n_dims = 6)
fsce_small <- calc_kmeans(fsce, k = 6)

usethis::use_data(fsce_small, compress = "xz", overwrite = TRUE)