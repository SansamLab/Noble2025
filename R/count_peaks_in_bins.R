#' Count Peaks Overlapping Genomic Bins
#'
#' This function counts the number of peaks that overlap with each genomic bin.
#'
#' @param peaks_gRanges A `GRanges` object containing peak locations.
#' @param genomicRange_bins A `GRanges` object representing genomic bins.
#'
#' @return A numeric vector of peak counts per bin.
#' @importFrom GenomicRanges countOverlaps
#' @export
#' @family noble_peak_density_heatmap
#' @examples
#' bins <- create_bins_across_genomic_range("chr1", 100000, 200000, 5000)
#' peaks <- GRanges("chr1", IRanges(start = c(100100, 150000), width = 200))
#' count_peaks_in_bins(peaks, bins)
count_peaks_in_bins <- function(peaks_gRanges, genomicRange_bins){
  peak_counts_in_bins <- countOverlaps(genomicRange_bins, peaks_gRanges)
  return(peak_counts_in_bins)
}
