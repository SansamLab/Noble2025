#' Create Genome Track Plot for MTBP Binding Across Cell Cycle Phases
#'
#' Generates a genome browser-style figure for a specified genomic region,
#' visualizing MTBP log2 fold change signals and peak calls from G1, Early S,
#' Late S, and G2 sorted samples. Repli-seq data is included for comparison.
#'
#' @param range_to_plot UCSC-style genomic region to visualize (e.g. "chr1:65748366-67420714").
#' @param highlight_regions_starts Numeric vector of start coordinates for highlight regions.
#' @param highlight_regions_widths Numeric vector of widths for highlight regions.
#' @param mtbp_paths Named list of BAM and peak BED file paths for MTBP ChIP in each phase.
#'        Each element should be a list with elements `tx`, `input`, `peaks`.
#'        Must be named with: "G1", "EarlyS", "LateS", "G2".
#' @param repliseq_gr_rds_path Path to RDS file containing Repli-seq GRanges.
#' @param ideogram_csv_path Path to CSV file with ideogram data.
#' @param phase_colors Named vector of colors for each phase.
#' @param track_sizes Named numeric vector indicating relative heights of each track.
#' @param figure_width Width of output figure in inches. Default: 16.
#' @param figure_height Height of output figure in inches. Default: 16.
#'
#' @return A grid graphics object representing the genome track plot.
#' @import Gviz
#' @import grid
#' @import magrittr
#' @import GenomicRanges
#' @export

createMtbpPhasesTracksWithRepliseq <- function(
    range_to_plot = "chr1:65748366-67420714",
    highlight_regions_starts = c(65748366, 90349536, 100202878),
    highlight_regions_widths = c(1672350, 8656743, 1369815),
    mtbp_paths,
    repliseq_gr_rds_path,
    ideogram_csv_path,
    phase_colors = c(G1 = "#114b5F", EarlyS = "#1a936f", LateS = "#88d498", G2 = "#cc2936"),
    track_sizes = c(
      itrack = 0.4, axisTrack = 0.8,
      MTBP_G1_log2fc_dtrack = 1, MTBP_G1_peaks_track = 0.5,
      MTBP_EarlyS_log2fc_dtrack = 1, MTBP_EarlyS_peaks_track = 0.5,
      MTBP_LateS_log2fc_dtrack = 1, MTBP_LateS_peaks_track = 0.5,
      MTBP_G2_log2fc_dtrack = 1, MTBP_G2_peaks_track = 0.5,
      dTrack_Repliseq = 1
    ),
    figure_width = 16,
    figure_height = 16
) {
  setCustomGvizScheme()
  range_split <- strsplit(range_to_plot, ":|-")[[1]]
  chromosome <- range_split[1]
  start_pos <- as.numeric(range_split[2])
  end_pos <- as.numeric(range_split[3])

  itrack <- IdeogramTrack(
    genome = "hg38",
    bands = read.csv(ideogram_csv_path, col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain")),
    ucscChromosomeNames = FALSE, showId = FALSE
  )

  axisTrack <- GenomeAxisTrack(fontsize = 8, distFromAxis = 9)

  dTrack_Repliseq <- readRDS(repliseq_gr_rds_path) %>%
    DataTrack(
      name = "RepliSeq", type = "heatmap",
      gradient = c("white", "#A1A0A1", "#4B4A4B", "black"),
      ylim = c(5, 20), showTitle = FALSE, genome = "hg38"
    )

  all_tracks <- list(itrack = itrack, axisTrack = axisTrack)

  for (phase in names(phase_colors)) {
    paths <- mtbp_paths[[phase]]
    color <- phase_colors[[phase]]

    message("Creating tracks for phase: ", phase)

    # Log2FC track creation with error handling
    log2fc_track <- tryCatch({
      make1SampleLog2FcCoverageDataTrack(
        chromosome = chromosome,
        start = start_pos,
        end = end_pos,
        windowSize = 25000,
        stepSize = 5000,
        txBamFile = paths$tx,
        inBamFile = paths$input,
        trackName = paste("MTBP", phase, "log2FC"),
        HistogramColor = color,
        yLimits = c(-0.8, 0.8)
      )
    }, error = function(e) {
      warning(sprintf("Error generating log2FC track for %s: %s", phase, e$message))
      NULL
    })

    # Peak annotation track creation with error handling
    peaks_track <- tryCatch({
      AnnotationTrack(
        range = importBED(paths$peaks, chromosomesToImport = chromosome),
        name = paste("MTBP", phase, "peaks"),
        col = color, fill = color, col.line = color,
        showTitle = FALSE, genome = "hg38"
      )
    }, error = function(e) {
      warning(sprintf("Error generating peak track for %s: %s", phase, e$message))
      NULL
    })

    all_tracks[[paste0("MTBP_", phase, "_log2fc_dtrack")]] <- log2fc_track
    all_tracks[[paste0("MTBP_", phase, "_peaks_track")]] <- peaks_track
  }

  all_tracks[["dTrack_Repliseq"]] <- dTrack_Repliseq

  # Filter out NULL tracks that failed
  all_tracks <- all_tracks[!vapply(all_tracks, is.null, logical(1))]

  ht <- HighlightTrack(
    trackList = all_tracks[names(track_sizes)],
    start = highlight_regions_starts,
    width = highlight_regions_widths,
    chromosome = chromosome,
    fill = rep("transparent", length(highlight_regions_starts)),
    col = c("#ED6A5A"),
    inBackground = FALSE, lwd = 2
  )

  plotTracks(
    ht,
    chromosome = chromosome,
    from = start_pos,
    to = end_pos,
    sizes = as.numeric(track_sizes[names(all_tracks)]),
    showTitle = FALSE
  ) %>% grid.grabExpr(width = figure_width, height = figure_height)
}
