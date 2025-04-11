#' Create Genomic Bins Across a Range
#'
#' This function generates evenly spaced bins of a given size across a specified genomic range.
#'
#' @param range_chromosome Character. The chromosome name (e.g., "chr1").
#' @param range_start Numeric. The start position of the genomic range.
#' @param range_end Numeric. The end position of the genomic range.
#' @param binSize Numeric. The width of each bin.
#'
#' @return A `GRanges` object containing genomic bins.
#' @family noble_peak_density_heatmap
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @export
#' @examples
#' create_bins_across_genomic_range("chr1", 100000, 200000, 5000)
create_bins_across_genomic_range <- function(range_chromosome,
                                             range_start,
                                             range_end,
                                             binSize){
  starts <- seq(range_start, range_end - binSize + 1, by = binSize)
  genomic_range_bins <- GRanges(seqnames = range_chromosome,
                                ranges = IRanges(start = starts, width = binSize))
  return(genomic_range_bins)
}
