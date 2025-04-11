#' Define a Unified Reference Peak Set
#'
#' This function merges multiple peak sets into a single unified peak set,
#' reducing overlapping peaks into non-overlapping regions.
#'
#' @param peaks_gRanges_list A list of `GRanges` objects, each representing a set of peaks.
#'
#' @return A `GRanges` object representing the unified peak set.
#'
#' @import GenomicRanges
#'
#' @seealso \code{\link{make_peak_overlap_df}}, \code{\link{make_euler_plot_showing_peak_overlaps}}
#' @export
#' @examples
#' peaks_list <- list(
#'   GRanges(seqnames="chr1", IRanges(start=c(100, 300), width=50)),
#'   GRanges(seqnames="chr1", IRanges(start=c(120, 350), width=50))
#' )
#' unified_peaks <- define_unified_reference_peaks(peaks_list)
#' print(unified_peaks)
define_unified_reference_peaks <- function(peaks_gRanges_list){
  peaks_gRanges_list %>% GRangesList %>% unlist %>% reduce
}
