
#' Create Genome Track Plot for TRESLIN, MTBP, and Replication Features
#'
#' Generates a genome browser-style figure for a specified genomic region,
#' visualizing log2 fold change signals from TRESLIN and MTBP ChIP-seq,
#' peak annotations, and origin mapping (Repli-seq, SNS-seq, Ini-seq).
#' Allows for customizable track inclusion and sizing via a named vector input.
#'
#' @param range_to_plot UCSC-style genomic region to visualize (e.g. "chr1:40000000-120000000").
#' @param treslin_tx_bam_1_path Path to TRESLIN ChIP replicate 1 BAM file.
#' @param treslin_in_bam_1_path Path to input control for TRESLIN replicate 1.
#' @param treslin_tx_bam_2_path Path to TRESLIN ChIP replicate 2 BAM file.
#' @param treslin_in_bam_2_path Path to input control for TRESLIN replicate 2.
#' @param treslin_color Color for TRESLIN signal and peaks. Default: "#006494"
#' @param treslin_log2fc_ylimits Y-axis limits for TRESLIN log2FC signal. Default: c(-0.8, 0.8)
#' @param treslin_log2fc_step_size Step size for sliding window of TRESLIN signal. Default: 5000
#' @param treslin_log2fc_window_size Window size for summing reads (TRESLIN). Default: 25000
#' @param mtbp_tx_bam_1_path Path to MTBP ChIP replicate 1 BAM file.
#' @param mtbp_in_bam_1_path Path to input control for MTBP replicate 1.
#' @param mtbp_tx_bam_2_path Path to MTBP ChIP replicate 2 BAM file.
#' @param mtbp_in_bam_2_path Path to input control for MTBP replicate 2.
#' @param mtbp_color Color for MTBP signal and peaks. Default: "#55C1FF"
#' @param mtbp_log2fc_ylimits Y-axis limits for MTBP log2FC signal. Default: c(-0.8, 0.8)
#' @param mtbp_log2fc_step_size Step size for sliding window of MTBP signal. Default: 5000
#' @param mtbp_log2fc_window_size Window size for summing reads (MTBP). Default: 25000
#' @param treslin_peaks_bed_path Path to BED file with TRESLIN peaks.
#' @param mtbp_peaks_bed_path Path to BED file with MTBP peaks.
#' @param snsseq_bed_path Path to BED file for SNS-seq origins.
#' @param snsseq_track_color Color for SNS-seq annotation track. Default: "#C1666B"
#' @param iniseq_bed_path Path to BED file for Ini-seq origins.
#' @param iniseq_track_color Color for Ini-seq annotation track. Default: "#D4B483"
#' @param repliseq_gr_rds_path Path to RDS file containing Repli-seq GRanges.
#' @param ideogram_csv_path Path to CSV file with ideogram data (columns: chrom, chromStart, chromEnd, name, gieStain).
#' @param highlight_regions_starts Numeric vector of start coordinates for highlight regions.
#' @param highlight_regions_widths Numeric vector of widths for highlight regions.
#' @param figure_width Width of output figure in inches. Default: 16
#' @param figure_height Height of output figure in inches. Default: 16
#' @param genome_build Genome build identifier (e.g. "hg38"). Default: "hg38"
#' @param tracks_to_plot_with_sizes Named numeric vector. Names must match track objects created inside the function. Values control their heights. Tracks appear in the specified order.
#'
#' @return A grid graphics object representing the genome track plot.
#'
#' @import Gviz
#' @import grid
#' @import magrittr
#' @import utils
#' @export

createAsynchronousTracksWithRepliseq <- function(
    range_to_plot = "chr1:40000000-120000000",
    treslin_tx_bam_1_path,
    treslin_in_bam_1_path,
    treslin_tx_bam_2_path,
    treslin_in_bam_2_path,
    treslin_color = "#006494",
    treslin_log2fc_ylimits = c(-0.8, 0.8),
    treslin_log2fc_step_size = 5000,
    treslin_log2fc_window_size = 25000,
    mtbp_tx_bam_1_path,
    mtbp_in_bam_1_path,
    mtbp_tx_bam_2_path,
    mtbp_in_bam_2_path,
    mtbp_color = "#55C1FF",
    mtbp_log2fc_ylimits = c(-0.8, 0.8),
    mtbp_log2fc_step_size = 5000,
    mtbp_log2fc_window_size = 25000,
    treslin_peaks_bed_path,
    mtbp_peaks_bed_path,
    snsseq_bed_path,
    snsseq_track_color = "#C1666B",
    iniseq_bed_path,
    iniseq_track_color = "#D4B483",
    repliseq_gr_rds_path,
    ideogram_csv_path,
    highlight_regions_starts,
    highlight_regions_widths,
    figure_width = 16,
    figure_height = 16,
    genome_build = "hg38",
    tracks_to_plot_with_sizes = c(
      "itrack" = 0.4,
      "axisTrack" = 0.8,
      "MTBP_Lg2FC_dTrack" = 1,
      "MTBP_Unified_Peaks_Track" = 0.5,
      "TRESLIN_Lg2FC_dTrack" = 1,
      "TRESLIN_Unified_Peaks_Track" = 0.5,
      "IniSeqTrack" = 0.5,
      "SnsSeq_Track" = 0.5,
      "dTrack_Repliseq" = 1
    )
){
  setCustomGvizScheme()
  track_sizes <- as.numeric(tracks_to_plot_with_sizes)

  message("Parsing range: ", range_to_plot)
  range_to_plot_split <- strsplit(range_to_plot, ":|-") %>%
    unlist() %>%
    magrittr::set_names(c("chromosome", "start", "end"))

  chromosome <- range_to_plot_split[["chromosome"]]
  start_pos <- as.numeric(range_to_plot_split[["start"]])
  end_pos <- as.numeric(range_to_plot_split[["end"]])

  message("Chromosome: ", chromosome)
  message("Start: ", start_pos)
  message("End: ", end_pos)

  message("Loading ideogram from: ", ideogram_csv_path)
  ideogram_df <- read.csv(ideogram_csv_path)
  names(ideogram_df) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
  itrack <- IdeogramTrack(genome = genome_build, bands = ideogram_df, ucscChromosomeNames = FALSE, showId = FALSE)

  axisTrack <- GenomeAxisTrack(fontsize = 8, distFromAxis = 9)

  message("Loading Repli-seq from: ", repliseq_gr_rds_path)
  dTrack_Repliseq <- readRDS(repliseq_gr_rds_path) %>%
    DataTrack(
      name = "RepliSeq", type = "heatmap",
      gradient = c("white", "#A1A0A1", "#4B4A4B", "black"),
      ylim = c(5, 20), showTitle = FALSE, genome = genome_build
    )

  message("Importing SNS-seq BED from: ", snsseq_bed_path)
  SnsSeq_Track <- AnnotationTrack(
    range = importBED(snsseq_bed_path, chromosomesToImport = chromosome),
    name = "SNSSeq", col = snsseq_track_color,
    fill = snsseq_track_color, col.line = snsseq_track_color,
    showTitle = FALSE, genome = genome_build
  )

  message("Importing Ini-seq BED from: ", iniseq_bed_path)
  IniSeqTrack <- AnnotationTrack(
    range = importBED(iniseq_bed_path, chromosomesToImport = chromosome),
    name = "IniSeq", col = iniseq_track_color,
    fill = iniseq_track_color, col.line = iniseq_track_color,
    showTitle = FALSE, genome = genome_build
  )

  message("Importing MTBP peaks from: ", mtbp_peaks_bed_path)
  MTBP_Unified_Peaks_Track <- AnnotationTrack(
    range = importBED(mtbp_peaks_bed_path, chromosomesToImport = chromosome),
    name = "MTBP peaks", col = mtbp_color,
    fill = mtbp_color, col.line = mtbp_color,
    showTitle = FALSE, genome = genome_build
  )

  message("Importing TRESLIN peaks from: ", treslin_peaks_bed_path)
  TRESLIN_Unified_Peaks_Track <- AnnotationTrack(
    range = importBED(treslin_peaks_bed_path, chromosomesToImport = chromosome),
    name = "TRESLIN peaks", col = treslin_color,
    fill = treslin_color, col.line = treslin_color,
    showTitle = FALSE, genome = genome_build
  )

  message("Creating TRESLIN signal track...")
  TRESLIN_Lg2FC_dTrack <- make2SampleLog2FcCoverageDataTrack(
    chromosome = chromosome,
    start = start_pos,
    end = end_pos,
    windowSize = treslin_log2fc_window_size,
    stepSize = treslin_log2fc_step_size,
    txBamFile1 = treslin_tx_bam_1_path,
    inBamFile1 = treslin_in_bam_1_path,
    txBamFile2 = treslin_tx_bam_2_path,
    inBamFile2 = treslin_in_bam_2_path,
    trackName = "TRESLIN log2FC",
    HistogramColor = treslin_color,
    yLimits = treslin_log2fc_ylimits
  )

  message("Creating MTBP signal track...")
  MTBP_Lg2FC_dTrack <- make2SampleLog2FcCoverageDataTrack(
    chromosome = chromosome,
    start = start_pos,
    end = end_pos,
    windowSize = mtbp_log2fc_window_size,
    stepSize = mtbp_log2fc_step_size,
    txBamFile1 = mtbp_tx_bam_1_path,
    inBamFile1 = mtbp_in_bam_1_path,
    txBamFile2 = mtbp_tx_bam_2_path,
    inBamFile2 = mtbp_in_bam_2_path,
    trackName = "MTBP log2FC",
    HistogramColor = mtbp_color,
    yLimits = mtbp_log2fc_ylimits
  )

  all_tracks_list <- list(
    "itrack" = itrack, "axisTrack" = axisTrack,
    "MTBP_Lg2FC_dTrack" = MTBP_Lg2FC_dTrack,
    "MTBP_Unified_Peaks_Track" = MTBP_Unified_Peaks_Track,
    "TRESLIN_Lg2FC_dTrack" = TRESLIN_Lg2FC_dTrack,
    "TRESLIN_Unified_Peaks_Track" = TRESLIN_Unified_Peaks_Track,
    "IniSeqTrack" = IniSeqTrack, "SnsSeq_Track" = SnsSeq_Track,
    "dTrack_Repliseq" = dTrack_Repliseq
  )

  tracks_to_plot <- all_tracks_list[names(tracks_to_plot_with_sizes)]

  message("Creating highlight track and rendering figure...")
  ht <- HighlightTrack(
    trackList = tracks_to_plot,
    start = highlight_regions_starts,
    width = highlight_regions_widths,
    chromosome = chromosome,
    fill = rep("transparent", length(highlight_regions_starts)),
    col = c("#ED6A5A"),
    inBackground = FALSE, lwd = 2
  )

  track_plot <- plotTracks(
    ht,
    chromosome = chromosome,
    from = start_pos,
    to = end_pos,
    sizes = track_sizes,
    showTitle = FALSE
  ) %>% grid.grabExpr(width = figure_width, height = figure_height)

  message("Track plot complete.")
  return(track_plot)
}
