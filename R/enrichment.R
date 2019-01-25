#' Identify haircut signals enriched over empty droplet distribution
#'
#' @param fsce An object of class [`FunctionalSingleCellExperiment`]
#' @param unfiltered_mat matrix with haircut umi counts for all barcodes, or path to
#' unfiltered matrix able to be read by [`scrunchy::read_matrix()`]
#' @param umi_cutoff minimum umi counts used to define background. Droplets
#' with fewer than this number of UMIs are considered backgroud, default = 50
#' @param count_cutoff minimum UMI counts at each site required for testing, default = 10
#' @param qval_cutoff minimum q.value to consider as enriched,
#' default = 0.001
#' @param enrich_cutoff minimum log2 increase to consider as enriched,
#' default = 2
#' @param per_cluster Default = FALSE, calculate enrichment per cluster
#' @param cluster_type if per_cluster than clustering column to group by,
#' one of "k_cluster" or "leiden_cluster", defaults to "leiden_cluster".
#'
#' @return `fsce` with enrichment values added to haircut `rowData``
#'
#' @importFrom SummarizedExperiment rowData<-
#'
#' @export
calc_signal_enrichment <- function(fsce,
                                   unfiltered_mat,
                                   umi_cutoff = 50,
                                   count_cutoff = 10,
                                   qval_cutoff = 0.001,
                                   enrich_cutoff = 2,
                                   per_cluster = FALSE,
                                   cluster_type = "leiden_cluster") {

  unfiltered_mat <- apply_if_path(unfiltered_mat, read_matrix)

  # get background distribution
  bc_ranks <- rank_barcodes(unfiltered_mat)
  bg_barcodes <- filter(bc_ranks, counts < umi_cutoff)
  bg_barcodes <- bg_barcodes[["cell_id"]]

  bg_umi_counts <- Matrix::rowSums(unfiltered_mat[, bg_barcodes])
  bg_umi_counts <- tibble(feature_id = names(bg_umi_counts),
                          bg_umi_count = bg_umi_counts)

  bg_umi_counts <- mutate(bg_umi_counts,
                          bg_total_count = sum(bg_umi_count))

  # get counts per site from haircut data including all colData
  hc_counts <- tidy_counts_full(fsce, "haircut")

  if (per_cluster) {
    # group by adduct and cluster
    hc_counts <- group_by(hc_counts,
                          !!sym(cluster_type),
                          !!sym("feature_id"),
                          !!sym("hairpin"))
  } else {
    # group just by each adduct type
    hc_counts <- group_by(hc_counts,
                          !!sym("feature_id"),
                          !!sym("hairpin"))
  }

  # calc total UMIs per group
  hc_counts <- summarize(hc_counts,
                         cell_counts = sum(counts))
  hc_counts <- ungroup(hc_counts)

  hc_counts <- mutate(hc_counts,
                      total_cell_counts = sum(cell_counts))

  # combine for stats test
  hc_counts <- left_join(hc_counts,
                         bg_umi_counts,
                         by = "feature_id")

  # run hypergeometric
  hc_site_stats <- mutate(hc_counts,
                          stats = purrr::pmap(
                            list(
                              cell_counts,
                              bg_umi_count,
                              bg_total_count - bg_umi_count,
                              total_cell_counts
                            ),
                            tidy_hyper
                          ))

  hc_site_stats <- tidyr::unnest(hc_site_stats)
  hc_site_stats <- calc_qvalues(hc_site_stats)

  hc_site_stats <- mutate(
    hc_site_stats,
    enriched_over_background = cell_counts > count_cutoff &
      q.value < qval_cutoff &
      enrichment > enrich_cutoff
  )

  if(per_cluster) {
    # report clusters with enriched signal as csv field
    hc_site_stats <- group_by(hc_site_stats,
                              !!sym("feature_id"))

    # only report clusters that are enriched in csv field
    hc_site_stats <- summarize(hc_site_stats,
                               enriched_over_background = any(enriched_over_background),
                               enriched_clusters = paste0(.data[[cluster_type]][enriched_over_background],
                                                          collapse = ","))
    hc_site_stats <- ungroup(hc_site_stats)
    hc_site_stats <- mutate(hc_site_stats,
                            enriched_clusters = ifelse(enriched_clusters == "",
                            NA,
                            enriched_clusters))

    keep_cols <- c("feature_id", "enriched_over_background", "enriched_clusters")
    hc_site_stats <- select(hc_site_stats, one_of(keep_cols))
  } else {
    hc_site_stats <- select(hc_site_stats,
                            one_of(c("feature_id", "enriched_over_background")))
  }

  # reorder to match rowData feature_id
  hc_site_stats <- left_join(tidy_rowdata(fsce, "haircut")[, "feature_id"],
            hc_site_stats,
            by = "feature_id")

  hc_site_stats <- select(hc_site_stats, -!!sym("feature_id"))

  rowData(fsce[["haircut"]]) <- cbind(rowData(fsce[["haircut"]]),
                                      hc_site_stats)

  fsce
}

# Utilities ---------------------------------------------------------

#' Rank barcodes by UMI count
#'
#' @param mat matrix with umi counts per cell
#' @noRd
rank_barcodes <- function(mat) {
  umi_counts <- Matrix::colSums(mat)
  bc_ranks <- tibble(cell_id = names(umi_counts),
                     counts = umi_counts)
  bc_ranks <- arrange(bc_ranks, desc(counts))
  bc_ranks <- mutate(bc_ranks, bc_rank = dplyr::row_number())
  bc_ranks
}

#' Tidyer for phyper
#'
#' @param success_in_sample # of successes in sample
#' @param success_in_population # of success in population
#' @param fail_in_population # of failures in population
#' @param sample_size sample size
#' @noRd
tidy_hyper <- function(success_in_sample,
                       success_in_population,
                       fail_in_population,
                       sample_size) {
  pval <- phyper(
    success_in_sample  - 1,
    success_in_population,
    fail_in_population,
    sample_size,
    lower.tail = FALSE
  )

  sample_prop <- success_in_sample / sample_size
  bg_prop <-  success_in_population /
    (success_in_population + fail_in_population)
  enrich <- log2(sample_prop / bg_prop )

  tibble(p.value = pval,
         enrichment = enrich)
}

#' If input is a path, then apply function
#'
#' @param input # of failures in population
#' @param .f function to apply
#' @param ... additional arguments to .f
#' @param verbose message which function is applied
#' @noRd
apply_if_path <- function(input, .f, ..., verbose = TRUE) {

  if (class(input) != "character") {
    return(input)
  }

  fxn_name <- deparse(substitute(.f))
  if (fs::dir_exists(input)) {

      if(verbose){
        message(glue::glue("applying {fxn_name} to {input}"))
      }

      res <- .f(input, ...)
      return(res)

  } else {
      stop(
        glue::glue(
          "{input} must be path or an object generated by {fxn_name}"
        ),
        call. = FALSE
      )
  }

}


