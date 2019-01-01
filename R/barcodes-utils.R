
#' Read molecule file into memory
#'
#' @param molecule_file molecule.tsv.gz flatfile produced by scrunchy pipeline
#' @param return_ids convert barcodes and features to names rather than indexes
#'   in barcodes.tsv.gz and features.tsv.gz (default = TRUE)
#'
#' @importFrom fs path_join path_dir file_exists
#' @importFrom readr read_tsv
#' @export
read_molecules <- function(molecule_file,
                           return_ids = TRUE) {

  if(!fs::file_exists(molecule_file)){
    stop(bcs_fn, " does not exist", call. = FALSE)
  }

  dat <- readr::read_tsv(molecule_file,
                         col_types = c('iici'),
                         col_names = c("barcode",
                                       "feature",
                                       "umi",
                                       "count"))


  if(return_ids){
    base_path <- fs::path_dir(molecule_file)
    bcs_fn <- fs::path_join(c(base_path, "barcodes.tsv"))
    features_fn <- fs::path_join(c(base_path, "features.tsv"))

    if(!fs::file_exists(bcs_fn)){
      stop(bcs_fn, " does not exist", call. = FALSE)
    }

    if(!fs::file_exists(features_fn)){
      stop(features_fn, " does not exist", call. = FALSE)
    }

    bcs <- readr::read_tsv(bcs_fn, col_names = "barcode_seq", col_types = "c")
    features <- readr::read_tsv(features_fn, col_names = "feature_name", col_types = "c")

    bcs[["idx"]] <- seq_len(nrow(bcs))
    features[["idx"]] <- seq_len(nrow(features))

    dat <- left_join(dat, bcs, by = c("barcode" = "idx"))
    dat <- left_join(dat, features, by = c("feature" = "idx"))

    dat <- select(dat,
                  barcode = barcode_seq,
                  feature = feature_name,
                  umi, count)
  }

  dat <- molecule_df(dat)
  dat
}

#' Test if the object is a molecule_df.
#'
#' @param x An object
#' @return `TRUE` if the object inherits from the [molecule_df()] class.
#' @export
is.molecule_df <- function(x) {
  "molecule_df" %in% class(x)
}

#' data.frame that contains molecule info
#'
#' @param x An object
#' @return `TRUE` if the object inherits from the [molecule_df()] class.
#' @export
molecule_df <- function(x) {
  expect_names <- c("barcode", "feature", "umi", "count")

  missing <- setdiff(expect_names, names(x))
  if (length(missing) != 0) {
    stop(sprintf(
      "expected %d required names, missing: %s",
      length(expected),
      paste0(missing, collapse = ", ")
    ))
  }

  class(x) <- union("molecule_df", class(x))
  x
}

#' Filter barcodes in molecule file and write to disk
#'
#' Useful for only keeping cell associated molecules for example.
#'
#' @param molecule_file path molecule.tsv.gz flatfile produced by scrunchy
#'   pipeline
#' @param output_file output file name, will be compressed by default.
#' @param bcs_to_use character vector of barcode sequences to use keep in output
#'
#' @importFrom readr write_tsv write_lines
#' @importFrom fs path
#' @export
filter_molecules <- function(molecule_file,
                             output_file,
                             bcs_to_keep) {

  dat <- read_molecules(molecule_file)

  dat <- filter(dat, barcode %in% bcs_to_keep)

  # convert features and barcode sequences back to indexes
  base_path <- fs::path_dir(molecule_file)
  out_path <- fs::path_dir(output_file)

  bcs_fn <- fs::path_join(c(base_path, "barcodes.tsv.gz"))
  features_fn <- fs::path_join(c(base_path, "features.tsv.gz"))

  bcs <- readr::read_tsv(bcs_fn, col_names = "barcode_seq", col_types = "c")
  features <- readr::read_tsv(features_fn, col_names = "feature_name", col_types = "c")

  bcs_subset <- filter(bcs, barcode_seq %in% dat$barcode)
  feature_subset <- filter(features, feature_name %in% dat$feature)

  # order based on cell barcode  and feature order in dat
  bcs_subset <- left_join(tibble(barcode_seq = unique(dat$barcode)),
                          bcs_subset,
                          by = "barcode_seq")

  # order based on cell barcode order in dat
  feature_subset <- left_join(tibble(feature_name = unique(dat$feature)),
                              feature_subset,
                          by = "feature_name")

  # make an index to store in molecules output
  bcs_subset[["bc_idx"]] <- seq_len(nrow(bcs_subset))
  feature_subset[["f_idx"]] <- seq(nrow(feature_subset))

  dat <- left_join(dat, bcs_subset, by = c("barcode" = "barcode_seq"))
  dat <- left_join(dat, feature_subset, by = c("feature" = "feature_name"))

  subset_dat <- select(dat,
                       bc_idx,
                       f_idx,
                       umi,
                       count)

  if(!fs::dir_exists(out_path)){
    fs::dir_create(out_path)
  }

  output_fn <- ifelse(endsWith(output_file, ".gz"),
                      output_file,
                      paste0(output_file, ".gz"))

  readr::write_tsv(subset_dat, output_fn, col_names = FALSE)


  # subset barcodes and features
  bc_dat <- unique(select(dat, barcode, bc_idx))
  f_dat <- unique(select(dat, feature, f_idx))

  bc_dat <- arrange(bc_dat, bc_idx)
  f_dat <- arrange(f_dat, f_idx)

  # write out barcodes and features to same directory as molecules
  readr::write_lines(bc_dat[["barcode"]],
                     fs::path_join(c(out_path, "barcodes.tsv.gz")))

  readr::write_lines(f_dat[["feature"]],
                     fs::path_join(c(out_path, "features.tsv.gz")))
}

#' Plot sequencing saturation
#'
#' Reads per cell are downsampled and the sequencing saturation is computed.
#' Sequencing saturation is defined as 1 - UMIs / total reads.
#'
#' @param molecules path molecule.tsv.gz flatfile produced by scrunchy pipeline
#'   or molecule_df data.frame
#' @param bcs_to_use character vector of barcode sequences to use for computing
#'   sequencing saturation. Defaults to use all barcodes.
#' @param proportions vector of proportions between 0 and 1 to downsample reads
#'   (defaults to `seq(0, 1, by = 0.1)`)
#' @param return_data if TRUE return plotting data instead of plotting
#'   distribution
#'
#' @importFrom tidyr gather
#' @importFrom scales comma
#' @export
plot_saturation <- function(molecules,
                             bcs_to_use = NULL,
                             proportions = seq(0, 1, by = 0.1),
                             return_data = FALSE) {

  # read from file unless is a molecule data.frame
  if(is.character(molecules)){
    message("reading in molecule file")
    molecules <- read_molecules(molecules)
  } else if (!is.molecule_df(molecules)){
    stop("molecules must be a molecule_df or a path to a molecules file", call. = FALSE)
  }

  dat <- select(molecules, -feature)

  if(!is.null(bcs_to_use)){
    dat <- filter(dat, barcode %in% bcs_to_use)
  }

  message("downsampling reads")

  # build list of rbinom function calls
  fxn_args <- paste0("~ rbinom(length(.), ., ", proportions, ")")

  lambdas <- map(fxn_args, as.formula)
  names(lambdas) <- paste0("downsampled_", proportions)

  # downsample reads per cell
  dat <- group_by(dat, barcode)
  umis_ds <- mutate_at(dat,
                       .vars = "count",
                       .funs = lambdas)

  message("calculating sequence saturation")
  # tidy data
  umis_ds <- ungroup(umis_ds)
  umis_ds <- select(umis_ds, -count)
  umis_ds <- tidyr::gather(umis_ds, proportion,
                           count, -barcode, -umi)

  # compute summaries
  umis_ds <- mutate(umis_ds, proportion = factor(proportion))
  umis_ds <- filter(umis_ds, count > 0)
  umis_ds <- group_by(umis_ds, proportion, barcode)
  umis_ds <- mutate(umis_ds, n_per_cell_reads = sum(count))
  umis_ds <- group_by(umis_ds, proportion)

  umis_ds_summary <- summarize(umis_ds,
                               mean_reads_per_cell = mean(n_per_cell_reads),
                               n_umi = n(),
                               n_reads = sum(count),
                               seq_saturation = 1 - (n_umi / n_reads))

  # handle NaNs produced by downsampling reads to 0
  umis_ds_summary[is.na(umis_ds_summary)] <- 0

  if(return_data) {
    return(umis_ds_summary)
  }

  p <- ggplot(umis_ds_summary, aes(mean_reads_per_cell,
                              seq_saturation)) +
    geom_line() +
    scale_x_continuous(labels = scales::comma) +
    cowplot::theme_cowplot() +
    labs(x = "Mean reads per cell",
         y = "Sequencing saturation")
  p

}

#' Plot UMI count distribution for barcodes
#'
#' @param molecules path molecule.tsv.gz flatfile produced by scrunchy pipeline
#'   or molecule_df data.frame
#' @param log_y plot with y axis as log
#' @param return_data if TRUE return plotting data instead of plotting
#'   distribution
#' @importFrom cowplot theme_cowplot
#' @importFrom scales comma
#' @export
plot_barcodes <- function(molecules,
                          log_y = TRUE,
                          return_data = FALSE) {

  # read from file unless is a molecule data.frame
  if(is.character(molecules)){
    message("reading in molecule file")
    molecules <- read_molecules(molecules)
  } else if (!is.molecule_df(molecules)){
    stop("molecules must be a molecule_df or a path to a molecules file", call. = FALSE)
  }

  dat <- group_by(molecules, barcode)
  umis <- summarize(dat, n_umi = n())
  umis <- ungroup(umis)
  umis <- arrange(umis, desc(n_umi))
  umis <- mutate(umis, bc_rank = row_number())

  if (return_data) {
    return(umis)
  }

  p <- ggplot(umis, aes(bc_rank, n_umi)) +
    geom_line() +
    scale_x_log10(labels = scales::comma) +
    cowplot::theme_cowplot() +
    labs(x = "Barcodes",
         y = "UMI counts")

  if(log_y){
    p <- p + scale_y_log10(labels = scales::comma)
  }

  p
}
