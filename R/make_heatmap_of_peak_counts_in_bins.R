#' Generate a Heatmap of Peak Counts in Genomic Bins
#'
#' This function creates a heatmap showing peak density across genomic bins.
#'
#' @param peak_counts_in_bins_df A data frame containing peak counts per bin.
#'   It must have two columns: `counts` (factor or numeric) and `x` (bin index).
#' @param heatmap_fill_colors A vector of colors corresponding to unique peak count values.
#' @param heatmap_breaks A numeric vector defining the breakpoints for color assignment.
#'
#' @return A `ggplot` object representing the heatmap.
#' @import ggplot2
#' @importFrom scales viridis_pal
#' @export
#' @seealso \code{\link{create_bins_across_genomic_range}}, \code{\link{count_peaks_in_bins}}
#' @family noble_peak_density_heatmap
#' @examples
#' bins_df <- data.frame(counts = c(0, 1, 2, 5, 3), x = 1:5)
#' make_heatmap_of_peak_counts_in_bins(bins_df, c("white", "red", "green", "blue"), c(0, 1, 2, 5))
make_heatmap_of_peak_counts_in_bins <- function(peak_counts_in_bins_df,
                                                heatmap_fill_colors,
                                                heatmap_breaks){
  cap_value <- max(heatmap_breaks)

  # Cap the counts if they exceed the highest break value
  if (!is.null(cap_value)) {
    peak_counts_in_bins_df$counts[peak_counts_in_bins_df$counts > cap_value] <- cap_value
  }

  # Convert counts to a factor to ensure discrete color mapping
  peak_counts_in_bins_df$counts <- factor(peak_counts_in_bins_df$counts,
                                          levels = heatmap_breaks)

  # Validate color and breaks match the number of levels
  if (length(heatmap_fill_colors) != length(levels(peak_counts_in_bins_df$counts))) {
    stop("Error: The number of colors must match the number of unique count values.")
  } else {
    message("Levels match the number of colors.")
  }

  heat_map_of_peak_counts <- ggplot(peak_counts_in_bins_df,
                                    aes(x = x, y = 1, fill = counts)) +
    geom_tile() +
    scale_fill_manual(values = setNames(heatmap_fill_colors, heatmap_breaks)) +
    theme_void() +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme(
          # legend.position = "none",
          plot.margin = margin(0, 0, 0, 0, "pt"))

  return(heat_map_of_peak_counts)
}
