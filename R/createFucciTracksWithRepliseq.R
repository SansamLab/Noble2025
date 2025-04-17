#' Create Genome Track Plot for MTBP Binding Across Cell Cycle Phases with Repli-seq and Peak Density
#'
#' This function generates a genome browser-style figure for a specified genomic region,
#' visualizing MTBP log2 fold change signals and peak densities from G1, Early S,
#' Late S, and G2 sorted samples. Repli-seq signal and ideogram tracks are included,
#' along with shaded highlight regions and optional filtering of blacklisted genomic regions.
#'
#' @param range_to_plot Genomic region to visualize in UCSC-style format (e.g. "chr1:65748366-67420714").
#' @param highlight_regions_starts Numeric vector of start positions for highlighted regions.
#' @param highlight_regions_widths Numeric vector of widths for each highlight region.
#' @param highlight_regions_colors Character vector of colors for highlight region outlines.
#' @param mtbp_paths A named list with one entry for each phase ("G1", "EarlyS", "LateS", "G2").
#'        Each entry must be a list with keys: `tx` (treatment BAM), `input` (control BAM), and `peaks` (BED file).
#' @param lg2fc_ylimits Numeric vector of two values defining Y-axis limits for log2 fold change tracks.
#' @param repliseq_gr_rds_path Path to an RDS file containing a GRanges object for Repli-seq signal.
#' @param ideogram_csv_path Path to a CSV file defining ideogram bands (columns: chrom, chromStart, chromEnd, name, gieStain).
#' @param blacklist_bed_path Path to BED file of genomic blacklist regions to exclude from log2FC calculation.
#' @param peak_heatplot_bin_size Bin size (in bp) to use when converting peak annotations to histogram-style density plots.
#' @param phase_colors Named vector of colors for each phase.
#' @param track_sizes Named numeric vector specifying the relative heights of tracks in the plot.
#' @param figure_width Width of output figure in inches.
#' @param figure_height Height of output figure in inches.
#' @param axis_track_scale The distance in base pairs between tick marks on the genome axis.
#'
#' @return A captured grid graphics object (using `grid.grabExpr`) representing the final genome track plot.
#'
#' @import Gviz
#' @import grid
#' @import magrittr
#' @import GenomicRanges
#' @export
createMtbpPhasesTracksWithRepliseq <- function(
    range_to_plot = "chr1:65748366-67420714", highlight_regions_starts = c(
      65748366,
      90349536, 100202878
    ), highlight_regions_widths = c(
      1672350,
      8656743, 1369815
    ), highlight_regions_colors = c("blue", "orange", "purple"),
    mtbp_paths, lg2fc_ylimits = c(-2, 1.5), repliseq_gr_rds_path, ideogram_csv_path,
    blacklist_bed_path, peak_heatplot_bin_size = 50000,
    phase_colors = c(
      G1 = "#114b5F", EarlyS = "#1a936f", LateS = "#88d498",
      G2 = "#cc2936"
    ), track_sizes = c(
      itrack = 0.4, axisTrack = 0.8,
      MTBP_G1_log2fc_dtrack = 1, MTBP_G1_peaks_track = 0.5,
      MTBP_EarlyS_log2fc_dtrack = 1, MTBP_EarlyS_peaks_track = 0.5,
      MTBP_LateS_log2fc_dtrack = 1, MTBP_LateS_peaks_track = 0.5,
      MTBP_G2_log2fc_dtrack = 1, MTBP_G2_peaks_track = 0.5,
      dTrack_Repliseq = 1
    ), figure_width = 16, figure_height = 16,
    axis_track_scale = 10000000) {
  setCustomGvizScheme()
  range_split <- strsplit(range_to_plot, ":|-")[[1]]
  chromosome <- range_split[1]
  start_pos <- as.numeric(range_split[2])
  end_pos <- as.numeric(range_split[3])
  itrack <- IdeogramTrack(genome = "hg38", bands = read.csv(ideogram_csv_path,
    col.names = c(
      "chrom", "chromStart", "chromEnd", "name",
      "gieStain"
    )
  ), ucscChromosomeNames = FALSE, showId = FALSE)
  axisTrack <- GenomeAxisTrack(
    fontsize = 8, labelPos = "below", from = start_pos, to = end_pos, scale = axis_track_scale,
    col.axis = "black", fontcolor = "black", col = "black", showTitle = TRUE
  )
  dTrack_Repliseq <- readRDS(repliseq_gr_rds_path) %>% DataTrack(
    name = "RepliSeq",
    type = "heatmap", gradient = c(
      "white", "#A1A0A1", "#4B4A4B",
      "black"
    ), ylim = c(5, 20), showTitle = FALSE, genome = "hg38"
  )
  all_tracks <- list(itrack = itrack)
  for (phase in names(phase_colors)) {
    paths <- mtbp_paths[[phase]]
    color <- phase_colors[[phase]]
    message("Creating tracks for phase: ", phase)
    log2fc_track <- make1SampleLog2FcCoverageDataTrack(
      chromosome = chromosome, blacklist_bed_path,
      start = start_pos, end = end_pos, windowSize = 50000,
      stepSize = 5000, txBamFile = paths$tx, inBamFile = paths$input,
      trackName = paste("MTBP", phase, "log2FC"), HistogramColor = color,
      yLimits = lg2fc_ylimits
    )

    genomic_range_bins <- create_bins_across_genomic_range(
      chromosome,
      start_pos,
      end_pos,
      peak_heatplot_bin_size
    )
    peak_counts_in_bins <- count_peaks_in_bins(
      importBED(paths$peaks, chromosomesToImport = chromosome), genomic_range_bins
    )
    genomic_range_bins$peak_counts <- peak_counts_in_bins
    peaks_track <- DataTrack(
      genomic_range_bins,
      name = paste("MTBP", phase, "peaks"), fill = color, col.histogram = "transparent",
      type = "hist", showTitle = FALSE, genome = "hg38", ylim = c(0, 8)
    )
    all_tracks[[paste0("MTBP_", phase, "_log2fc_dtrack")]] <- log2fc_track
    all_tracks[[paste0("MTBP_", phase, "_peaks_track")]] <- peaks_track
  }
  all_tracks[["dTrack_Repliseq"]] <- dTrack_Repliseq
  all_tracks[["axisTrack"]] <- axisTrack
  all_tracks <- all_tracks[!vapply(all_tracks, is.null, logical(1))]
  ht <- HighlightTrack(
    trackList = all_tracks[names(track_sizes)],
    start = highlight_regions_starts, width = highlight_regions_widths,
    chromosome = chromosome, fill = rep("transparent", length(highlight_regions_starts)),
    col = highlight_regions_colors, inBackground = FALSE, lwd = 0.5
  )
  plotTracks(ht,
    chromosome = chromosome, from = start_pos,
    to = end_pos, sizes = as.numeric(track_sizes),
    showTitle = FALSE
  ) %>% grid.grabExpr(
    width = figure_width,
    height = figure_height
  )
}
