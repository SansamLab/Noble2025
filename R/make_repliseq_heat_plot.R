#' Create a Repli-seq heatmap for a genomic region
#'
#' This function generates a heatmap of Repli-seq signal intensity across a specified
#' genomic region. It extracts Repli-seq signal from a preloaded GRanges object
#' (assumed to be named \code{HCT_Repliseq_gr}) and plots the signal across S-phase fractions.
#'
#' The heatmap has genomic coordinates on the x-axis and S-phase fraction bins (S1–S16) on the y-axis.
#'
#' @param chromosomeToPlot Chromosome name (e.g., "chr1").
#' @param startToPlot Start coordinate (integer).
#' @param endToPlot End coordinate (integer).
#'
#' @return A \code{ggplot2} object representing the Repli-seq heatmap.
#'
#' @details The function expects a \code{GRanges} object named \code{HCT_Repliseq_gr}
#' to be available, where metadata columns S1 to S16 contain replication signal across S-phase.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom scales oob_squish
#' @export
make_repliseq_heat_plot <- function(chromosomeToPlot, startToPlot, endToPlot) {
  # Load Repli-seq data (must be available as HCT_Repliseq_gr)
  RepSeqData <- data(HCT_Repliseq_gr)

  # Subset for the desired chromosome and coordinates, then transform into a plotting-ready dataframe
  RepSeqData_subset <- RepSeqData[
    seqnames(RepSeqData) == chromosomeToPlot &
      start(RepSeqData) >= startToPlot &
      end(RepSeqData) <= endToPlot
  ] %>%
    as.data.frame() %>%
    {
      # Compute midpoint for each genomic bin
      .$x <- round(rowMeans(data.frame(.$start, .$end)))
      .
    } %>%
    # Keep only midpoint and S-phase fraction signal columns
    .[c("x", paste0("S", 1:16))] %>%
    # Convert from wide to long format
    melt(id.vars = "x") %>%
    {
      # Rename columns to generic plotting names
      names(.) <- c("X", "Y", "Z")
      .
    } %>%
    {
      # Ensure y-axis shows S16 at top and S1 at bottom
      .$Y <- factor(.$Y, levels = paste0("S", 16:1))
      .
    }

  # Generate the heatmap plot
  htmp <- ggplot(RepSeqData_subset, aes(X, Y, fill = Z)) +
    geom_tile() +
    scale_fill_gradient(
      limits = c(5, 20),           # Clamp fill scale to [5, 20]
      low = "white", high = "black",
      oob = scales::oob_squish     # Prevent color scale warnings
    ) +
    theme_void() +                 # Remove axes, ticks, etc.
    theme(
      legend.position = "none",   # Suppress color legend
      plot.margin = margin(0, 0, 0, 0)
    ) +
    scale_x_continuous(expand = c(0, 0))

  return(htmp)
}
