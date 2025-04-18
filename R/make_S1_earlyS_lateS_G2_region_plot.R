#' Create a multi-track genome browser-style plot for a genomic region
#'
#' This function generates a genome browser-style figure showing log2 signal-to-background
#' ratios across a specified genomic region, along with a Repli-seq heatmap.
#'
#' Each track represents a different cell cycle stage, using input and treatment BAM files
#' and optional peak annotations. The plot includes a subregion marker line and custom y-axis limits.
#'
#' @param chromosomeToPlot Chromosome name as a character string (e.g., "chr1").
#' @param startToPlot Start position of the genomic region as a numeric value.
#' @param endToPlot End position of the genomic region as a numeric value.
#' @param windowSize Size of sliding window for signal smoothing (default: 50000).
#' @param stepSize Step size between windows (default: 5000).
#' @param YLimits Numeric vector of length 2 specifying y-axis limits for log2(Signal/Input).
#' @param Bams_Beds_and_Colors_List A named list of lists, one per condition (e.g., "G1", "EarlyS").
#'   Each inner list must include:
#'   \describe{
#'     \item{txBamFile}{Path to treatment BAM file.}
#'     \item{inBamFile}{Path to input BAM file.}
#'     \item{peaksBedFile}{Path to peaks BED file (unused in this function but included for completeness).}
#'     \item{track_color}{Color for the plotted track.}
#'   }
#' @param subregion_marker_color Color of horizontal line marker to indicate subregion of interest (default: "red").
#'
#' @return A ggplot2 object composed using cowplot::plot_grid, combining genome tracks and a heatmap.
#'
#' @import ggplot2
#' @import cowplot
#' @export
make_S1_earlyS_lateS_G2_region_plot <- function(
    chromosomeToPlot, startToPlot, endToPlot,
    windowSize = 50000, stepSize = 5000,
    YLimits = c(-1.7, 1.4),
    Bams_Beds_and_Colors_List = list(
      "G1" = list(
        txBamFile = "G1_tx.bam",
        inBamFile = "G1_in.bam",
        peaksBedFile = "peaks.bed",
        track_color = G1_Color
      ),
      "EarlyS" = list(
        txBamFile = "EarlyS_tx.bam",
        inBamFile = "EarlyS_in.bam",
        peaksBedFile = "peaks.bed",
        track_color = EarlyS_Color
      ),
      "LateS" = list(
        txBamFile = "LateS_tx.bam",
        inBamFile = "LateS_in.bam",
        peaksBedFile = "peaks.bed",
        track_color = LateS_Color
      ),
      "G2" = list(
        txBamFile = "G2_tx.bam",
        inBamFile = "G2_in.bam",
        peaksBedFile = "peaks.bed",
        track_color = G2_Color
      )
    ),
    subregion_marker_color = "red"
) {
  chromosomeToPlot <- as.character(chromosomeToPlot)
  startToPlot <- as.numeric(startToPlot)
  endToPlot <- as.numeric(endToPlot)

  Lg2FCPlots <- lapply(Bams_Beds_and_Colors_List, function(Bam_Bed_Color) {
    makeLog2FcGenomePlot(
      chrom = chromosomeToPlot,
      trackStart = startToPlot,
      trackEnd = endToPlot,
      windowSize = windowSize,
      stepSize = stepSize,
      txBamFile = Bam_Bed_Color[[1]],
      inBamFile = Bam_Bed_Color[[2]],
      autoCalculateYLimits = FALSE,
      yLimits = YLimits,
      fillColor = Bam_Bed_Color[[4]]
    )
  })

  htmp <- makeRepliSeqHeatPlot(
    chromosomeToPlot,
    startToPlot,
    endToPlot
  )

  theme_genome_track <- theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      axis.text.y = element_text(color = "black", size = 6),
      axis.title.y = element_text(color = "black", size = 8, angle = 90, vjust = 1),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.1, linetype = "solid")
    )

  Lg2FCPlots <- lapply(Lg2FCPlots, function(plt) {
    plt + theme_genome_track + ylab("Log2\nSig/Bkg")
  })

  heights <- c(rep(1, length(Bams_Beds_and_Colors_List)), 1)

  Tracks <- plot_grid(
    Lg2FCPlots[[1]] + geom_hline(yintercept = YLimits[2], color = subregion_marker_color, linewidth = 2),
    Lg2FCPlots[[2]],
    Lg2FCPlots[[3]],
    Lg2FCPlots[[4]],
    htmp,
    align = "v",
    ncol = 1,
    rel_heights = heights
  )

  return(Tracks)
}
