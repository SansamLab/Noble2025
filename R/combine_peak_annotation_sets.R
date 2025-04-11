
#' Combine multiple labeled peak annotation sets into a single data.frame
#'
#' Calls `annotate_and_label_peaks()` for each set and binds the results
#' into a single tidy data frame.
#'
#' @param peak_sets Named list of GRanges objects representing different peak sets.
#' @param txdb A TxDb object for gene annotation.
#' @param tss_region A numeric vector of length 2 for TSS region. Default: c(-3000, 3000)
#'
#' @return A data.frame containing annotation summaries for all peak sets.
#' @import dplyr
#' @export
combine_peak_annotation_sets <- function(peak_sets, txdb, tss_region = c(-3000, 3000)) {
  annotation_dfs <- lapply(names(peak_sets), function(set_name) {
    annotate_and_label_peaks(
      peaks = peak_sets[[set_name]],
      set_name = set_name,
      txdb = txdb,
      tss_region = tss_region
    )
  })

  dplyr::bind_rows(annotation_dfs)
}
