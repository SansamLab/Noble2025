#' Generate an Euler Plot of Peak Overlaps
#'
#' This function creates an Euler plot based on a data frame of peak overlaps.
#'
#' @param peak_overlap_df A data frame where each row represents a peak in the reference set,
#'        and each column corresponds to a query peak set, with `TRUE`/`FALSE` values indicating overlaps.
#' @param plot_colors A vector of colors to use for the Euler plot.
#'
#' @return A `ggplot` object representing the Euler plot.
#'
#' @import eulerr
#' @export
#' @seealso \code{\link{make_peak_overlap_df}}, \code{\link{make_euler_plot_showing_peak_overlaps}}
#'
#' @examples
#' overlap_df <- data.frame(Sample1 = c(TRUE, FALSE, TRUE), Sample2 = c(FALSE, TRUE, TRUE))
#' colors <- c("#1b9e77", "#d95f02")
#' euler_plot <- make_euler_plot_of_overlaps_with_reference_peaks(overlap_df, colors)
#' print(euler_plot)
make_euler_plot_of_overlaps_with_reference_peaks <- function(peak_overlap_df, plot_colors){
  fit <- eulerr::euler(peak_overlap_df)
  eulerPlot <- plot(fit,
                    fills = plot_colors,
                    labels = list(fontsize = 8),
                    quantities = list(fontsize = 6))
  return(eulerPlot)
}
