#' Create MTBP Phase Genome Tracks with Repli-seq Overlay
#'
#' This function generates genome browser-style tracks showing MTBP log2 fold change signals
#' for each cell cycle phase (G1, Early S, Late S, G2), peak density heatmaps, and Repli-seq data.
#' Ideograms and highlighted genomic regions are overlaid for context.
#'
#' @param range_to_plot UCSC-style genomic coordinates (e.g., \"chr1:65748366-67420714\").
#' @param highlight_regions List containing 'ranges' (vector of UCSC strings) and 'colors'.
#' @param log2fc_track_parameters Named list specifying parameters for each cell cycle phase.
#' @param auto_calculate_y_limits Logical indicating whether Y-axis limits are auto-calculated.
#' @param colors_for_heat_map_function Function mapping signals to colors (e.g., colorRamp2).
#' @param breaks_for_heat_map Numeric vector defining breaks for the heatmap.
#' @param blacklist_bed_path Path to a blacklist BED file (currently unused).
#' @param additional_ggplot_theme Additional ggplot theme adjustments.
#'
#' @importFrom ggbio plotIdeogram
#' @importFrom cowplot plot_grid
#' @importFrom GenomicRanges seqnames start end
#' @importFrom tidyr pivot_longer
#' @importFrom circlize colorRamp2
#' @importFrom ggplot2 theme element_text element_line geom_vline aes scale_color_identity theme_void ylab
#' @return A list with ideogram plot and genome tracks as cowplot objects.
#' @export
createMtbpPhasesTracksWithRepliseq2 <- function(
    range_to_plot = "chr1:65748366-67420714",
    highlight_regions = list(
      ranges = c("chr1:65748366-67420714", "chr1:90349536-99016989", "chr1:100202878-101572693"),
      colors = c("blue", "orange", "purple")
    ),
    log2fc_track_parameters = list(
      G1 = list(
        txBam = "G1_tx.bam", inBam = "G1_in.bam", color = "#114b5F",
        window_bp = 50000, step_bp = 5000, lg2fc_ylimits = c(-2, 1.5), peaks_bed = "peaks.bed"
      ),
      EarlyS = list(
        txBam = "EarlyS_tx.bam", inBam = "EarlyS_in.bam", color = "#1a936f",
        window_bp = 50000, step_bp = 5000, lg2fc_ylimits = c(-2, 1.5), peaks_bed = "peaks.bed"
      ),
      LateS = list(
        txBam = "LateS_tx.bam", inBam = "LateS_in.bam", color = "#88d498",
        window_bp = 50000, step_bp = 5000, lg2fc_ylimits = c(-2, 1.5), peaks_bed = "peaks.bed"
      ),
      G2 = list(
        txBam = "G2_tx.bam", inBam = "G2_in.bam", color = "#cc2936",
        window_bp = 50000, step_bp = 5000, lg2fc_ylimits = c(-2, 1.5), peaks_bed = "peaks.bed"
      )
    ),
    auto_calculate_y_limits = FALSE,
    colors_for_heat_map_function = circlize::colorRamp2(c(0, 5), c("transparent", "#3f007d")),
    breaks_for_heat_map = c(0, 1, 2, 3, 4, 5),
    blacklist_bed_path = NULL,
    additional_ggplot_theme = theme(
      axis.text.y = element_text(color = "black", size = 6, hjust = 1),
      axis.title.y = element_text(color = "black", size = 8, angle = 90, vjust = 1),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.1, linetype = "solid")
    )) {
  range_granges <- parseUCSCtoGRanges(range_to_plot)

  ideogram_plot <- ggbio::plotIdeogram(
    hg38_ideogram_granges,
    which = range_granges,
    fill = "transparent",
    zoom.offset = 4,
    main = ""
  )
  ideogram_plot@subchr <- NULL

  overlay_plots <- lapply(names(log2fc_track_parameters), function(stage) {
    params <- log2fc_track_parameters[[stage]]
    makeLog2FcGenomePlot(
      chrom = as.character(seqnames(range_granges)),
      trackStart = start(range_granges),
      trackEnd = end(range_granges),
      windowSize = params$window_bp,
      stepSize = params$step_bp,
      txBamFile = params$txBam,
      inBamFile = params$inBam,
      autoCalculateYLimits = auto_calculate_y_limits,
      yLimits = params$lg2fc_ylimits,
      fillColor = params$color,
      blacklist_granges = hg38_blacklist_v2_granges
    )
  })
  names(overlay_plots) <- names(log2fc_track_parameters)

  peakGRs <- lapply(log2fc_track_parameters, function(prms) importBED(prms$peaks_bed))
  peakGRs <- lapply(peakGRs, function(gr) {
    start(gr) <- start(gr) + 1
    gr
  })

  peak_density_heat_maps <- lapply(peakGRs, function(gr) {
    generate_peak_density_heatmap(
      peaks_gRanges = gr,
      range_chromosome = as.character(seqnames(range_granges)),
      range_start = start(range_granges),
      range_end = end(range_granges),
      binSize = 100000,
      colors = colors_for_heat_map_function(breaks_for_heat_map),
      breaks = breaks_for_heat_map
    ) + theme(legend.position = "none")
  }) %>% setNames(names(log2fc_track_parameters))

  regions_df <- highlight_regions$ranges %>%
    lapply(function(strng) strsplit(gsub(",", "", strng), "[:-]")[[1]]) %>%
    do.call(rbind, .) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    setNames(c("chrom", "start", "end")) %>%
    mutate(across(c(start, end), as.numeric),
           colors = highlight_regions$colors
    ) %>%
    pivot_longer(cols = c(start, end), values_to = "x")

  subregion_layer <- list(
    geom_vline(data = regions_df, aes(xintercept = x, color = colors), inherit.aes = FALSE),
    scale_color_identity()
  )

  genome_tracks <- plot_grid(
    overlay_plots[["G1"]] + theme_void() + additional_ggplot_theme + ylab("Log2\nSig/Bkg") + subregion_layer,
    peak_density_heat_maps[["G1"]],
    overlay_plots[["EarlyS"]] + theme_void() + additional_ggplot_theme + ylab("Log2\nSig/Bkg") + subregion_layer,
    peak_density_heat_maps[["EarlyS"]],
    overlay_plots[["LateS"]] + theme_void() + additional_ggplot_theme + ylab("Log2\nSig/Bkg") + subregion_layer,
    peak_density_heat_maps[["LateS"]],
    overlay_plots[["G2"]] + theme_void() + additional_ggplot_theme + ylab("Log2\nSig/Bkg") + subregion_layer,
    peak_density_heat_maps[["G2"]],
    NULL,
    make_repliseq_heat_plot(as.character(seqnames(range_granges)), start(range_granges), end(range_granges)),
    align = "v", ncol = 1, rel_heights = c(rep(c(1, 0.1), 4), 0.5, 1)
  )

  return(list(ideogram = ideogram_plot@ggplot, genome_tracks = genome_tracks))
}
