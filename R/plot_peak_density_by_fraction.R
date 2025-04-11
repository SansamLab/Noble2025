
#' Plot Peak Density by Replication Timing Fraction
#'
#' Generates a line plot showing peak density across S phase replication timing fractions.
#'
#' @param df A data frame returned from \code{calculate_peak_density_by_fraction()}.
#' @param title Optional title for the plot. Default is \code{NULL}.
#'
#' @return A \code{ggplot} object.
#'
#' @import ggplot2
#' @export
#'
plot_peak_density_by_fraction <- function(df, title = NULL) {
  ggplot(df, aes(x = S.fraction.numeric, y = peakDensity)) +
    geom_line() +
    geom_point(size = 1) +
    theme_bw() +
    sansam_theme() +
    xlab("Rep. timing") +
    ylab("Peaks per 50 kb") +
    ggtitle(title)
}
