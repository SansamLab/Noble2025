
#' Generate a Two-Sample Euler Plot from GRanges Objects
#'
#' This function creates a reduced union of two GRanges objects, computes overlaps,
#' and plots a two-set Euler diagram using specified fill colors.
#'
#' @param gr1 GRanges. First genomic range set.
#' @param gr2 GRanges. Second genomic range set.
#' @param labels Character vector of length 2. Names for the sets (default: c("Set1", "Set2")).
#' @param fills Character vector of length 2. Fill colors for the Euler plot (default: c("grey", "white")).
#'
#' @return A ggplot2 object showing the Euler diagram.
#' @examples
#' plot <- plot_two_sample_euler(kumagai_peaks_gr, noble_peaks_gr, labels = c("Kumagai", "Noble"))
#' print(plot)
#'
#' @import GenomicRanges
#' @import IRanges
#' @import eulerr
#' @import magrittr
#'
#' @export
plot_two_sample_euler <- function(gr1, gr2, labels = c("Set1", "Set2"), fills = c("grey", "white")) {
  # Ensure GRanges inputs
  stopifnot(inherits(gr1, "GRanges"), inherits(gr2, "GRanges"))

  # Create union of input GRanges and reduce overlaps
  reduced_peaks_gr <- reduce(c(gr1, gr2))

  # Create binary overlap matrix
  attendance_test_df <- data.frame(
    reduced_peaks_gr %over% gr1,
    reduced_peaks_gr %over% gr2
  ) %>% set_names(.,labels)

  # Fit Euler diagram
  euler_fit <- eulerr::euler(attendance_test_df)

  # Plot Euler diagram
  euler_plot <- plot(
    euler_fit,
    fills = fills,
    labels = list(fontsize = 8),
    quantities = list(fontsize = 6)
  )

  return(euler_plot)
}
