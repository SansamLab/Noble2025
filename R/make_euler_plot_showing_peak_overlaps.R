#' Generate an Euler Plot Showing Peak Overlaps from Multiple BED Files
#'
#' This function reads multiple BED files, defines a unified reference peak set,
#' determines peak overlaps, and creates an Euler plot.
#'
#' @param peaks_bed_paths_list A named list of file paths to BED files, where each represents a peak set.
#' @param colors_vector A vector of colors to use for the Euler plot.
#'
#' @return A `ggplot` object representing the Euler plot of peak overlaps.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom eulerr euler
#' @export
#' @seealso \code{\link{define_unified_reference_peaks}}, \code{\link{make_peak_overlap_df}}, \code{\link{make_euler_plot_of_overlaps_with_reference_peaks}}
#'
#' @examples
#' peaks_bed_files <- list(
#'   "sample1" = "path/to/sample1.bed",
#'   "sample2" = "path/to/sample2.bed"
#' )
#' colors <- c("#1b9e77", "#d95f02")
#' euler_plot <- make_euler_plot_showing_peak_overlaps(peaks_bed_files, colors)
#' print(euler_plot)
make_euler_plot_showing_peak_overlaps <- function(peaks_bed_paths_list, colors_vector){
  peaks_gRanges_list <- lapply(peaks_bed_paths_list, readBed)
  names(peaks_gRanges_list) <- names(peaks_bed_paths_list)

  reference_peaks <- define_unified_reference_peaks(peaks_gRanges_list)
  overlap_df <- make_peak_overlap_df(reference_peaks, peaks_gRanges_list)

  euler_plot <- make_euler_plot_of_overlaps_with_reference_peaks(overlap_df, colors_vector)
  return(euler_plot)
}

