
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scrunchy <img src="man/figures/logo.png" align="right" alt="" width="120" />

[![Travis build
status](https://travis-ci.org/hesselberthlab/scrunchy.svg?branch=master)](https://travis-ci.org/hesselberthlab/scrunchy)
[![Codecov test
coverage](https://codecov.io/gh/hesselberthlab/scrunchy/branch/master/graph/badge.svg)](https://codecov.io/gh/hesselberthlab/scrunchy?branch=master)
[![Lifecycle
badge](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/scrunchy)](https://cran.r-project.org/package=scrunchy)

scrunchy provides analysis tools for the **s**ingle-**c**ell
**r**econstruction of f**unc**tional **h**eterogeneit**y**.

New methods to study heterogeneity at cellular resolution in complex
tissues are rapidly transforming human biology. These methods measure
differences in gene expression, chromatin accessibility, and protein
levels across thousands of cells to understand developmental
trajectories of tissues, tumors, and whole organisms. However, their
reliance on measurements of steady-state abundance of DNA, RNA, and
protein limits our ability to extract dynamic information from single
cells.

To propel the study of heterogeneity among single cells, we are
developing functional assays as a new modality for single-cell
experiments. Instead of measuring the molecular abundance of DNA, RNA,
or protein in single cells and predicting functional states, our key
innovation is to directly measure enzymatic activities in single cells
by analyzing the conversion of substrates to products by single-cell
extracts in a high-throughput DNA sequencing experiment.

### Functional heterogeneity of DNA repair in immune cells

Our first functional method simultaneously measures the activity of DNA
repair enzymes and the abundance of mRNAs from thousands of single
cells. We measure DNA repair activities by encapsulating synthetic DNA
oligonucleotides with defined lesions with single cells. Cellular DNA
repair enzymes recognize and catalyze incision of these substrates,
which we subsequently capture in a [modified library construction
protocol](https://www.protocols.io/view/single-cell-analysis-of-functional-heterogeneity-i-uhyet7w).

An example data set in scrunchy contains a subset of data from an
experiment in which we simultaneously measure mRNA expression and
specific types of DNA repair activities in thousands of single cells
(human primary mononuclear blood cells).

The following plot shows cells classified by mRNA expression for each
cell type embedded in a two-dimensional UMAP projection. This data set
is too small to robustly classify cell types, but the full data set
(containing \~4,000 cells) indicates that clusters 1, 2, 3, and 5 are
T-cell subtypes.

<img src="man/figures/README-mrna_plot-1.png" width="100%" />

The following plots show DNA repair activities measured for each cluster
above. The repair substrates included an unmodified DNA substrate
(**A**), an A:U base-pair (**B**), a ribonucleotide (rG-D, **C**), and
an abasic site (**D**). Activity is measured as a normalized count of
incisions observed at expected repair positions.

These data show that whereas little repair activity is associated with
the unmodified DNA (**A**), uracil base excision (**B**), ribonucleotide
excision (**C**), and abasic site processing (**D**) activities can be
measured in single cells. Moreover, differences in activities between
clusters may reflect inherent differences in the levels of these DNA
repair pathways.

<img src="man/figures/README-activity_plot-1.png" width="100%" />

### Installation

scrunchy is under active development. You can install the R package from
github:

``` r
# install.packages("remotes")
remotes::install_github("hesselberthlab/scrunchy")
```
