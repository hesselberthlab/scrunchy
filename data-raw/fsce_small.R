library(scrunchy)
library(usethis)

fsce_small <- create_fsce(
  list(
    rnaseq = create_sce_rnaseq(scrunchy_data("mrna/")),
    haircut = create_sce_haircut(scrunchy_data("haircut/"))
  )
)

usethis::use_data(fsce_small, compress = "xz")
