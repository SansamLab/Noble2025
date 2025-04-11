
#' Annotate and label a set of genomic peaks
#'
#' Uses ChIPseeker::annotatePeak to annotate peaks and adds a column
#' indicating the peak set.
#'
#' @param peaks A GRanges object representing peaks to annotate.
#' @param set_name Character string labeling the source of the peaks.
#' @param txdb A TxDb object (e.g. from GenomicFeatures::makeTxDbFromGFF).
#' @param tss_region A numeric vector of length 2 specifying the TSS region.
#' Default: c(-3000, 3000)
#'
#' @return A data.frame of annotation statistics with a new column "Set"
#' @import ChIPseeker
#' @import dplyr
annotate_and_label_peaks <- function(peaks, set_name, txdb, tss_region = c(-3000, 3000)) {
  ChIPseeker::annotatePeak(peaks, tssRegion = tss_region, TxDb = txdb, annoDb = "org.Hs.eg.db")@annoStat %>%
    dplyr::mutate(Set = set_name)
}
