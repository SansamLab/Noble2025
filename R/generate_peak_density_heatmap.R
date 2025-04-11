#' Generate a Peak Density Heatmap
#'
#' This function creates a heatmap representing peak density across genomic bins.
#'
#' @param peaks_gRanges A `GRanges` object containing peak locations.
#' @param range_chromosome Character. The chromosome name (e.g., "chr1").
#' @param range_start Numeric. The start position of the genomic range.
#' @param range_end Numeric. The end position of the genomic range.
#' @param binSize Numeric. The width of each bin.
#' @param colors Character vector of colors corresponding to discrete peak count values.
#'   Defaults to a `viridis` color scale based on the number of breakpoints.
#' @param breaks Numeric vector defining the breakpoints for color assignment.
#'   Defaults to `c(0:4)`, meaning counts are capped at 4.
#'
#' @return A `ggplot` object representing the heatmap.
#'
#' @details This function performs the following steps:
#'   1. Creates genomic bins across the specified range using \code{\link{create_bins_across_genomic_range}}.
#'   2. Counts peaks overlapping each bin using \code{\link{count_peaks_in_bins}}.
#'   3. Converts the count vector into a data frame using \code{\link{convert_vector_to_dataframe_for_plot}}.
#'   4. Generates a heatmap with \code{\link{make_heatmap_of_peak_counts_in_bins}}.
#'
#' @family noble_peak_density_heatmap
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @import ggplot2
#' @importFrom viridisLite viridis
#' @export
#' @seealso \code{\link{create_bins_across_genomic_range}}, \code{\link{count_peaks_in_bins}},
#'   \code{\link{convert_vector_to_dataframe_for_plot}}, \code{\link{make_heatmap_of_peak_counts_in_bins}}
#'
#' @examples
#' # Example GRanges object with peaks
#' peaks <- GRanges("chr1", IRanges(start = c(100100, 150000, 180000), width = 200))
#'
#' # Generate a heatmap for the given genomic range
#' generate_peak_density_heatmap(peaks, "chr1", 100000, 200000, binSize = 5000,
#'                               colors = c("white", "red", "green", "blue"),
#'                               breaks = c(0, 1, 2, 3))
generate_peak_density_heatmap <- function(peaks_gRanges,
                                          range_chromosome,
                                          range_start,
                                          range_end,
                                          binSize,
                                          colors = viridis(length(breaks)),
                                          breaks = c(0:4)) {
  # Generate genomic bins
  genomic_range_bins <- create_bins_across_genomic_range(range_chromosome,
                                                         range_start,
                                                         range_end, binSize)

  # Count peaks in bins
  peak_counts_in_bins <- count_peaks_in_bins(peaks_gRanges,
                                             genomic_range_bins)

  # Convert counts to a data frame for plotting
  peak_counts_in_bins_df <- convert_vector_to_dataframe_for_plot(peak_counts_in_bins)

  # Generate the heatmap
  heat_map_of_peak_counts <- make_heatmap_of_peak_counts_in_bins(peak_counts_in_bins_df,
                                                                 colors,
                                                                 breaks)
  return(heat_map_of_peak_counts)
}
