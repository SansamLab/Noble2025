# Load required libraries
library(magrittr)
library(GenomicRanges)
library(parallel)
library(bamsignals)
library(Rsamtools)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)
library(viridisLite)
library(ChIPseeker)
library(eulerr)
library(UpSetR)
library(corrplot)
library(zoo)


#' Create a Data Frame Indicating Peak Overlaps
#'
#' This function takes a set of subject peaks and a list of query peak sets,
#' then determines whether each query peak set overlaps with the subject peaks.
#'
#' @param subjectPeaks A `GRanges` object representing the reference peak set.
#' @param queryPeaksList A named list of `GRanges` objects, each representing a query peak set.
#'
#' @return A data frame where each row corresponds to a peak in `subjectPeaks`,
#'         and each column corresponds to a query peak set. Values are `TRUE` if the
#'         subject peak overlaps the query peak set and `FALSE` otherwise.
#'
#' @importFrom GenomicRanges `%over%`
#'
#' @seealso \code{\link{define_unified_reference_peaks}}, \code{\link{make_euler_plot_of_overlaps_with_reference_peaks}}
#'
#' @examples
#' subjectPeaks <- GRanges(seqnames="chr1", IRanges(start=c(100,200), width=50))
#' queryPeaksList <- list(
#'   Sample1 = GRanges(seqnames="chr1", IRanges(start=c(120, 250), width=50)),
#'   Sample2 = GRanges(seqnames="chr1", IRanges(start=c(180, 300), width=50))
#' )
#' overlap_df <- make_peak_overlap_df(subjectPeaks, queryPeaksList)
#' print(overlap_df)
make_peak_overlap_df <- function(subjectPeaks, queryPeaksList){
  logical_overlap_vector_list <- lapply(queryPeaksList, function(query_gr){
    subjectPeaks %over% query_gr
  })
  logical_overlap_df <- do.call(cbind, logical_overlap_vector_list)
  logical_overlap_df <- as.data.frame(logical_overlap_df)
  colnames(logical_overlap_df) <- names(queryPeaksList)
  return(logical_overlap_df)
}

#' Define a Unified Reference Peak Set
#'
#' This function merges multiple peak sets into a single unified peak set,
#' reducing overlapping peaks into non-overlapping regions.
#'
#' @param peaks_gRanges_list A list of `GRanges` objects, each representing a set of peaks.
#'
#' @return A `GRanges` object representing the unified peak set.
#'
#' @import GenomicRanges
#'
#' @seealso \code{\link{make_peak_overlap_df}}, \code{\link{make_euler_plot_showing_peak_overlaps}}
#'
#' @examples
#' peaks_list <- list(
#'   GRanges(seqnames="chr1", IRanges(start=c(100, 300), width=50)),
#'   GRanges(seqnames="chr1", IRanges(start=c(120, 350), width=50))
#' )
#' unified_peaks <- define_unified_reference_peaks(peaks_list)
#' print(unified_peaks)
define_unified_reference_peaks <- function(peaks_gRanges_list){
  peaks_gRanges_list %>% GRangesList %>% unlist %>% reduce
}

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
#'
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

#' Generate an Euler Plot Showing Peak Overlaps from Multiple BED Files
#'
#' This function reads multiple BED files, defines a unified reference peak set,
#' determines peak overlaps, and creates an Euler plot.
#'
#' @param peaks_bed_paths_list A named list of file paths to BED files, where each represents a peak set.
#' @param colors_vector A vector of colors to use for the Euler plot.
#'
#' @return A `ggplot` object representing the Euler plot of peak overlaps.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom eulerr euler
#'
#' @seealso \code{\link{define_unified_reference_peaks}}, \code{\link{make_peak_overlap_df}}, \code{\link{make_euler_plot_of_overlaps_with_reference_peaks}}
#'
#' @examples
#' peaks_bed_files <- list(
#'   "sample1" = "path/to/sample1.bed",
#'   "sample2" = "path/to/sample2.bed"
#' )
#' colors <- c("#1b9e77", "#d95f02")
#' euler_plot <- make_euler_plot_showing_peak_overlaps(peaks_bed_files, colors)
#' print(euler_plot)
make_euler_plot_showing_peak_overlaps <- function(peaks_bed_paths_list, colors_vector){
  peaks_gRanges_list <- lapply(peaks_bed_paths_list, readBed)
  names(peaks_gRanges_list) <- names(peaks_bed_paths_list)

  reference_peaks <- define_unified_reference_peaks(peaks_gRanges_list)
  overlap_df <- make_peak_overlap_df(reference_peaks, peaks_gRanges_list)

  euler_plot <- make_euler_plot_of_overlaps_with_reference_peaks(overlap_df, colors_vector)
  return(euler_plot)
}


#' Annotate Genomic Features Overlapping Peaks
#'
#' This function annotates genomic regions overlapping with peaks using ChIPseeker.
#'
#' @param peaks_gRanges A `GRanges` object containing genomic peak locations.
#' @param annotationDb Character. The annotation database for gene annotations. Defaults to `"org.Hs.eg.db"` for human genes.
#' @param TxDb A `TxDb` object specifying transcript annotations. Defaults to `txdb`.
#' @param SetName Character. The name of the peak set (used for labeling in downstream analyses). Defaults to `"All Peaks"`.
#'
#' @return A data frame summarizing the genomic features overlapping the given peaks.
#'   The data frame includes:
#'   - `Feature`: Genomic feature (e.g., "Promoter", "Exon", "Intron").
#'   - `Frequency`: The number of peaks overlapping each feature.
#'   - `Set`: The label for the peak set.
#'
#' @import ChIPseeker
#' @importFrom GenomicFeatures TxDb
#' @importFrom AnnotationDbi select
#'
#' @seealso \code{\link{make_stacked_barchart_of_genomic_features_overlapping_peaks}}
#'
#' @examples
#' library(GenomicRanges)
#' peaks <- GRanges("chr1", IRanges(start = c(100100, 150000, 180000), width = 200))
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' df <- make_df_of_genomic_features_overlapping_peaks(peaks, TxDb = txdb)
#' head(df)
make_df_of_genomic_features_overlapping_peaks <- function(
  peaks_gRanges, annotationDb = "org.Hs.eg.db", TxDb = txdb, SetName = "All Peaks"
) {
  # Annotate peaks with genomic features
  ChiPSeeker_annotatePeaks_output <- annotatePeak(
    peaks_gRanges,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = annotationDb
  )

  # Extract annotation statistics
  df_of_genomic_features_overlapping_peaks <- ChiPSeeker_annotatePeaks_output@annoStat
  df_of_genomic_features_overlapping_peaks$Set <- SetName

  return(df_of_genomic_features_overlapping_peaks)
}

#' Create a Stacked Bar Chart of Genomic Features Overlapping Peaks
#'
#' This function generates a stacked bar chart representing the frequency of genomic
#' features overlapping with peaks.
#'
#' @param df_of_genomic_features_overlapping_peaks A data frame containing the frequency of genomic features overlapping with peaks.
#'   Expected columns:
#'   - `Feature`: Genomic feature (e.g., "Promoter", "Exon", "Intron").
#'   - `Frequency`: The number of peaks overlapping each feature.
#'   - `Set`: The label for the peak set.
#' @param feature_colors A named vector of colors for each genomic feature.
#'   The default color scheme follows a standardized palette for promoters, UTRs, exons, introns, and intergenic regions.
#'
#' @return A `ggplot` object representing a stacked bar chart of genomic feature frequencies.
#'
#' @import ggplot2
#' @importFrom grid unit
#'
#' @seealso \code{\link{make_df_of_genomic_features_overlapping_peaks}}
#'
#' @examples
#' df <- data.frame(
#'   Feature = c("Promoter (<=1kb)", "Exon", "Intron", "Distal Intergenic"),
#'   Frequency = c(100, 200, 300, 150),
#'   Set = c("Peaks Set 1", "Peaks Set 1", "Peaks Set 1", "Peaks Set 1")
#' )
#' make_stacked_barchart_of_genomic_features_overlapping_peaks(df)
make_stacked_barchart_of_genomic_features_overlapping_peaks <- function(
  df_of_genomic_features_overlapping_peaks,
  feature_colors = c(
    "Promoter (<=1kb)" = "#9b2226",
    "Promoter (1-2kb)" = "#ae2012",
    "Promoter (2-3kb)" = "#bb3e03",
    "5' UTR" = "#ca6702",
    "3' UTR" = "#ee9b00",
    "1st Exon" = "#e9d8a6",
    "Other Exon" = "#F7F1DE",
    "1st Intron" = "#94d2bd",
    "Other Intron" = "#0a9396",
    "Downstream (<=300)" = "#005f73",
    "Distal Intergenic" = "#001219"
  )
) {
  # Ensure the data frame has the required columns
  if (!all(c("Feature", "Frequency", "Set") %in% colnames(df_of_genomic_features_overlapping_peaks))) {
    stop("Error: Input data frame must contain 'Feature', 'Frequency', and 'Set' columns.")
  }

  # Generate stacked bar chart
  FeaturePlot <- ggplot(df_of_genomic_features_overlapping_peaks, aes(fill = Feature, y = Frequency, x = Set)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = feature_colors) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      text = element_text(size = 8),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 8),
      legend.key.size = unit(2, 'mm'),
      panel.background = element_rect(fill = 'transparent'),
      plot.background = element_rect(fill = 'transparent', color = NA)
    ) +
    coord_flip()

  return(FeaturePlot)
}



#' Create Genomic Bins Across a Range
#'
#' This function generates evenly spaced bins of a given size across a specified genomic range.
#'
#' @param range_chromosome Character. The chromosome name (e.g., "chr1").
#' @param range_start Numeric. The start position of the genomic range.
#' @param range_end Numeric. The end position of the genomic range.
#' @param binSize Numeric. The width of each bin.
#'
#' @return A `GRanges` object containing genomic bins.
#' @family noble_peak_density_heatmap
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @examples
#' create_bins_across_genomic_range("chr1", 100000, 200000, 5000)
create_bins_across_genomic_range <- function(range_chromosome,
                                             range_start,
                                             range_end,
                                             binSize){
  starts <- seq(range_start, range_end - binSize + 1, by = binSize)
  genomic_range_bins <- GRanges(seqnames = range_chromosome,
                                ranges = IRanges(start = starts, width = binSize))
  return(genomic_range_bins)
}

#' Count Peaks Overlapping Genomic Bins
#'
#' This function counts the number of peaks that overlap with each genomic bin.
#'
#' @param peaks_gRanges A `GRanges` object containing peak locations.
#' @param genomicRange_bins A `GRanges` object representing genomic bins.
#'
#' @return A numeric vector of peak counts per bin.
#' @importFrom GenomicRanges countOverlaps
#' @family noble_peak_density_heatmap
#' @examples
#' bins <- create_bins_across_genomic_range("chr1", 100000, 200000, 5000)
#' peaks <- GRanges("chr1", IRanges(start = c(100100, 150000), width = 200))
#' count_peaks_in_bins(peaks, bins)
count_peaks_in_bins <- function(peaks_gRanges, genomicRange_bins){
  peak_counts_in_bins <- countOverlaps(genomicRange_bins, peaks_gRanges)
  return(peak_counts_in_bins)
}

#' Convert a Numeric Vector to a Data Frame for Plotting
#'
#' This function converts a numeric vector (e.g., peak counts) into a data frame
#' that can be used for visualization.
#'
#' @param count_vector Numeric vector representing counts (e.g., peak counts per bin).
#'
#' @return A data frame with two columns:
#'   - `counts`: The peak count per bin.
#'   - `x`: The bin index (1-based).
#' @family noble_peak_density_heatmap
#' @examples
#' convert_vector_to_dataframe_for_plot(c(1, 3, 5, 2, 0))
convert_vector_to_dataframe_for_plot <- function(count_vector){
  count_data_frame <- data.frame(counts = count_vector, x = 1:length(count_vector))
  return(count_data_frame)
}

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

#' Generate a Peak Density Heatmap
#'
#' This function creates a heatmap representing peak density across genomic bins.
#'
#' @param peaks_gRanges A `GRanges` object containing peak locations.
#' @param range_chromosome Character. The chromosome name (e.g., "chr1").
#' @param range_start Numeric. The start position of the genomic range.
#' @param range_end Numeric. The end position of the genomic range.
#' @param binSize Numeric. The width of each bin.
#' @param colors Character vector of colors corresponding to discrete peak count values.
#'   Defaults to a `viridis` color scale based on the number of breakpoints.
#' @param breaks Numeric vector defining the breakpoints for color assignment.
#'   Defaults to `c(0:4)`, meaning counts are capped at 4.
#'
#' @return A `ggplot` object representing the heatmap.
#'
#' @details This function performs the following steps:
#'   1. Creates genomic bins across the specified range using \code{\link{create_bins_across_genomic_range}}.
#'   2. Counts peaks overlapping each bin using \code{\link{count_peaks_in_bins}}.
#'   3. Converts the count vector into a data frame using \code{\link{convert_vector_to_dataframe_for_plot}}.
#'   4. Generates a heatmap with \code{\link{make_heatmap_of_peak_counts_in_bins}}.
#'
#' @family noble_peak_density_heatmap
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @import ggplot2
#' @importFrom viridisLite viridis
#'
#' @seealso \code{\link{create_bins_across_genomic_range}}, \code{\link{count_peaks_in_bins}},
#'   \code{\link{convert_vector_to_dataframe_for_plot}}, \code{\link{make_heatmap_of_peak_counts_in_bins}}
#'
#' @examples
#' # Example GRanges object with peaks
#' peaks <- GRanges("chr1", IRanges(start = c(100100, 150000, 180000), width = 200))
#'
#' # Generate a heatmap for the given genomic range
#' generate_peak_density_heatmap(peaks, "chr1", 100000, 200000, binSize = 5000,
#'                               colors = c("white", "red", "green", "blue"),
#'                               breaks = c(0, 1, 2, 3))
generate_peak_density_heatmap <- function(peaks_gRanges,
                                          range_chromosome,
                                          range_start,
                                          range_end,
                                          binSize,
                                          colors = viridis(length(breaks)),
                                          breaks = c(0:4)) {
  # Generate genomic bins
  genomic_range_bins <- create_bins_across_genomic_range(range_chromosome,
                                                         range_start,
                                                         range_end, binSize)

  # Count peaks in bins
  peak_counts_in_bins <- count_peaks_in_bins(peaks_gRanges,
                                             genomic_range_bins)

  # Convert counts to a data frame for plotting
  peak_counts_in_bins_df <- convert_vector_to_dataframe_for_plot(peak_counts_in_bins)

  # Generate the heatmap
  heat_map_of_peak_counts <- make_heatmap_of_peak_counts_in_bins(peak_counts_in_bins_df,
                                                                 colors,
                                                                 breaks)
  return(heat_map_of_peak_counts)
}


#' Generate a Matrix of Smoothed Read Counts Across Genomic Regions
#'
#' This function creates a count matrix representing read density across a set of genomic ranges,
#' smoothed using a rolling sum within each chromosome. It is primarily used to generate signal
#' profiles (e.g., for metaplots or heatmaps) from BAM files using midpoint-aligned paired-end reads.
#'
#' @param bamPath Character. Path to the input BAM file.
#' @param Ranges A \code{GRanges} object representing the genomic regions of interest.
#' @param plottingRange Integer. Total width to which each region should be resized (default: 4001).
#' @param windowSize Integer. Width of the smoothing window (default: 150).
#' @param stepSize Integer. Step size for signal binning (default: 5).
#' @param coreCount Integer. Number of parallel processes to use (default: 8).
#' @param chromosomes Character vector of chromosome names to include (default: \code{allowedChromosomes}).
#'
#' @return A numeric matrix where rows represent bins within each genomic region and columns represent regions.
#'
#' @details Each region in \code{Ranges} is resized to \code{plottingRange} bp centered on its midpoint.
#' Reads are binned using \code{stepSize}, and the signal is smoothed using a rolling window
#' of width \code{windowSize / stepSize}.
#'
#' @importFrom GenomicRanges resize seqnames
#' @importFrom parallel mclapply
#' @importFrom Rsamtools idxstatsBam
#' @importFrom zoo rollsumr
#' @importFrom bamsignals bamProfile alignSignals
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gr <- GRanges("chr1", IRanges(start = c(100000, 200000), width = 4001))
#' mat <- makeCountMatrix("mydata.bam", gr)
#' }
makeCountMatrix <- function(
    bamPath = "../2024Dec30_siTICRRonMTBP_CutRun/ReplicatePeakAnalyzer/results/mergedDownSampledBams/controlMTBP_in.bam",
    Ranges = CombinedRanges,
    plottingRange = 4001,
    windowSize = 150,
    stepSize = 5,
    coreCount = 8,
    chromosomes = allowedChromosomes
){
  Ranges <- resize(Ranges,width=plottingRange,fix="center")
  countMatrix <- mclapply(paste0("chr", 1:22), function(chrm){
    chrom_ranges <- Ranges[seqnames(Ranges) == chrm]
    result <- bamProfile(bamPath, chrom_ranges, binsize = stepSize, paired.end = "midpoint", verbose = TRUE) %>%
      alignSignals(.) %>%
      as.matrix(.) %>%
      rollsumr(., k = windowSize / stepSize, align = "center") + 0.0001
    return(result)
  },
  mc.cores = coreCount) %>%
  do.call(cbind, .)
  return(countMatrix)
}

#' Compute Total Mapped Reads for Canonical Chromosomes
#'
#' This function calculates the total number of mapped reads in a BAM file,
#' restricted to canonical chromosomes (chr1–chr22).
#'
#' @param BamFile Character. Path to a BAM file.
#'
#' @return A single numeric value representing the total number of reads mapped
#'   to chromosomes chr1 through chr22.
#'
#' @details This function uses \code{idxstatsBam()} to obtain alignment statistics
#' for all reference sequences in the BAM file, then filters to autosomes
#' and returns the sum of mapped reads.
#'
#' @importFrom Rsamtools idxstatsBam
#'
#' @export
#'
#' @examples
#' \dontrun{
#' total_reads <- getTotalReadCount("example.bam")
#' }
getTotalReadCount <- function(BamFile){
  idxstatsBam(BamFile) %>%
      {
        .[.$seqnames %in% paste0("chr", c(1:22)), ]
      } %>%
      {
        sum(as.numeric(.[, 3]))
      }
}

#' Create a plot from BED files
#'
#' This function reads a list of BED files, calculates overlaps, and creates a plot
#' based on the specified plot choice (UpSet or Euler). It can be used to visualize
#' genomic overlap patterns between the BED files.
#'
#' @param BedFilenames A character vector containing paths to the BED files.
#' @param plotChoice A character string specifying the type of plot to create ("UpSet" or "Euler").
#' @param sampleLabelsDF A data frame with file labels and corresponding file paths.
#' @param colrs A list of colors for the plot elements.
#' @param alfa Alpha transparency value for plot elements.
#' @param fontScale Scaling factor for text size in the plot.
#' @param minFracOverlap Minimum fraction of overlap required for inclusion in the plot.
#'
#' @return NULL (The function generates and displays the plot without returning a value.)
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce findOverlaps width pintersect seqnames start end queryHits subjectHits
#' @import UpSetR euler plot
#' @export
#'
#' @examples
#' BedFilenames <- c("path/to/bedfile1.bed", "path/to/bedfile2.bed")
#' plotChoice <- "UpSet"
#' sampleLabelsDF <- data.frame(files = BedFilenames, labels = c("Label1", "Label2"))
#' colrs <- list(colors = c("#1b9e77", "#d95f02"))
#' alfa <- 0.7
#' fontScale <- 8L
#' minFracOverlap <- 0L
#' makePlotFromBeds(BedFilenames, plotChoice, sampleLabelsDF, colrs, alfa, fontScale, minFracOverlap)
#'
makePlotFromBeds <- function(BedFilenames, plotChoice, sampleLabelsDF, colrs, alfa, fontScale, minFracOverlap) {
  # Read and convert BedFilenames to GenomicRanges
  Beds.gr <- lapply(BedFilenames, function(file) {
    bed.df <- read.table(file)
    GenomicRanges::makeGRangesFromDataFrame(bed.df, seqnames.field = "V1", start.field = "V2", end.field = "V3")
  })

  # Reduce to get AllBeds.gr
  AllBeds.gr <- GenomicRanges::reduce(do.call("c", Beds.gr))

  # Function to get overlaps
  GetOverlapsWithAll <- function(bed.gr) {
    hits <- GenomicRanges::findOverlaps(AllBeds.gr, bed.gr)
    sigHits <- hits[GenomicRanges::width(GenomicRanges::pintersect(AllBeds.gr[queryHits(hits)], bed.gr[subjectHits(hits)])) / GenomicRanges::width(AllBeds.gr[queryHits(hits)]) > minFracOverlap]
    GRovrlps <- AllBeds.gr[queryHits(sigHits)]
    paste(GenomicRanges::seqnames(GRovrlps), GenomicRanges::start(GRovrlps), GenomicRanges::end(GRovrlps), sep = "_")
  }

  # Get overlaps for all BedFilenames
  lst <- lapply(Beds.gr, GetOverlapsWithAll)

  # Match sampleLabelsDF and remove file extensions
  nmes <- gsub("\\.[^.]*$", "", sampleLabelsDF$labels[match(sampleLabelsDF$files, basename(BedFilenames))])

  # Set names for lst
  names(lst) <- nmes

  # Create upsetList and plot
  upsetList <- UpSetR::fromList(lst)
  if (plotChoice == "UpSet") {
    UpSetR::upset(upsetList, nsets = length(lst), order.by = "freq", text.scale = fontScale)
  } else {
    plot(eulerr::euler(upsetList), quantities = TRUE, fills = colrs$colors, alpha = alfa)
  }
}



#' Read and convert a BED file to a GenomicRanges object
#'
#' This function reads a BED file and converts it into a GenomicRanges object
#' using the GenomicRanges package. The BED file should have columns named
#' "seqnames", "start", and "end" to represent genomic ranges.
#' @import GenomicRanges
#' @param bedFile Path to the BED file to be read and converted.
#'
#' @return A GenomicRanges object representing the data in the BED file.
#'
#' @examples
#' bedFile <- "path/to/your/bedfile.bed"
#' gr <- readBed(bedFile)
readBed <- function(bedFile){
  bedData <- read.table(bedFile)
  if(!is.numeric(bedData[1,2])){
    bedData <- bedData[-1,]
  }
  bedData <- set_names(bedData, c("seqnames", "start", "end"))
  return(makeGRangesFromDataFrame(bedData))
}

#' Create a unified peak set from multiple peak files and write it to a file
#'
#' This function takes a list of peak files and creates a unified peak set
#' by combining the second set of peaks with non-overlapping peaks from the first set.
#' The resulting peak set is then written to a specified output file in BED format.
#'
#' @param pkFiles A list of paths to the peak files to be combined.
#' @param outputPath Path to the output file where the unified peak set will be saved.
#'
#' @return NULL (The function writes the result to a file and doesn't return a value.)
#'
#' @examples
#' pkFiles <- c("path/to/peakfile1.bed", "path/to/peakfile2.bed")
#' outputPath <- "path/to/output/unified_peaks"
#' createAndWriteUnifiedPeakSet(pkFiles, outputPath)
#'
createAndWriteUnifiedPeakSet <- function(pkFiles, smtFiles,outputPath){
  # Read and convert peak files to GenomicRanges objects
  peaks <- lapply(pkFiles, readBed)
  summits <- lapply(smtFiles,readBed)

  # Combine the second set of peaks with non-overlapping peaks from the first set
  unifiedPeaks <- peaks[[2]][!peaks[[2]] %over% peaks[[1]]] %>% c(., peaks[[1]])
  unifiedSummits <- summits[[2]][!peaks[[2]] %over% peaks[[1]]] %>% c(., summits[[1]])

  # Create the output directory if it doesn't exist
  if (!file.exists(dirname(outputPath))) {
    dir.create(dirname(outputPath), recursive = TRUE)
  }

  # Write the unified peak set to the specified output file
  write.table(unifiedPeaks %>%
                as.data.frame %>%
                .[,1:3],
              paste0(outputPath,".bed"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  write.table(unifiedSummits %>%
                as.data.frame %>%
                .[,1:3],
              paste0(outputPath,"_summits.bed"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  mcols(unifiedPeaks)$summits <- start(unifiedSummits)
  saveRDS(unifiedPeaks,
          paste0(outputPath,".rds"))
}

#' Lift over genomic coordinates from hg19 to hg38
#'
#' This function takes a BED file in hg19 coordinates and performs a liftOver
#' to convert the coordinates to hg38 using the provided chain file.
#'
#' @param hg19BedInput Path to the input BED file in hg19 coordinates.
#' @param hg38BedOutput Path to the output BED file in hg38 coordinates.
#'
#' @details The function downloads the hg19 to hg38 liftover chain file from
#' UCSC if it's not already present in the working directory.
#'
#' @return None (output is saved to hg38BedOutput file).
#'
#' @import rtracklayer liftOver
#'
#' @examples
#' hg19ToHg38Liftover("input_hg19.bed", "output_hg38.bed")
#'
#' @export
hg19ToHg38Liftover <- function(hg19BedInput,hg38BedOutput){
  # Define the URL of the chain file
  chain_url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

  # Define the destination file path
  dest_file <- "hg19ToHg38.over.chain.gz"

  # Download the chain file
  if (!file.exists("hg19ToHg38.over.chain.gz")) {
    download.file(chain_url, destfile = dest_file, method = "auto")
    gunzip("hg19ToHg38.over.chain.gz")
  }

  bed_file <- hg19BedInput
  bed_df <- read.table(bed_file,header=F,sep="\t",stringsAsFactors = F)
  bed_gr <- GRanges(bed_df$V1,IRanges(start = bed_df$V2, end = bed_df$V3)) %>%
    {mcols(.)<-bed_df[,-c(1,2,3)];.}

  # Load the hg19 to hg38 liftover chain filename
  hg19_to_hg38_chain <- import.chain("hg19ToHg38.over.chain")

  # Convert hg19 genomic ranges object to hg38 using the chain file
  hg38_gr <- liftOver(bed_gr, hg19_to_hg38_chain) %>% unlist

  export(hg38_gr, hg38BedOutput, format = "bed")
}


#' Create a Complex Heatmap of Domain Log2 Fold Change
#'
#' Generates a Complex Heatmap of log2 fold change (Lg2FC) values for genomic domains.
#'
#' @param domainsBed Path to the BED file containing genomic domains.
#' @param windowSize The size of the sliding window used for calculating Lg2FC.
#' @param stepSize The step size for sliding the window.
#' @param txBamFile1 Path to the BAM file for treatment sample 1.
#' @param inBamFile1 Path to the BAM file for control sample 1.
#' @param txBamFile2 Path to the BAM file for treatment sample 2.
#' @param inBamFile2 Path to the BAM file for control sample 2.
#' @param span The span used for generating Lg2FC values.
#'
#' @return A ComplexHeatmap object representing the log2 fold change heatmap.
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom rtracklayer readBed idxstatsBam
#' @importFrom GenomicRanges split resize width order rev GRanges
#' @importFrom Rsamtools bamProfile
#' @importFrom BiocParallel mclapply
#' @importFrom zoo rollapply
#'
#' @examples
#' heatmap <- makeDomainHeatPlot(domainsBed = "path/to/domains.bed",
#'                               windowSize = 100000,
#'                               stepSize = 5000,
#'                               txBamFile1 = "path/to/txBam1.bam",
#'                               inBamFile1 = "path/to/inBam1.bam",
#'                               txBamFile2 = "path/to/txBam2.bam",
#'                               inBamFile2 = "path/to/inBam2.bam",
#'                               span = 3000000)
#'
#' @export
makeDomainHeatPlot <- function(domainsBed="../GenerateTimingWindowsFromHCTRepliseq/TreslinDepletedLateRepDomains.bed",
                               windowSize=100000,
                               stepSize=5000,
                               txBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                               inBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                               txBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_1000.bam",
                               inBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_1000.bam",
                               span=3000000){
  ranges_gr <- readBed(domainsBed) %>% .[order(width(.))] %>%
    .[width(.)>200000] %>%
    .[width(.)<2000000]
  #windws_lst <- slidingWindows(gr,windowSize,step=stepSize)
  makeCPMS <- function(BamFile=txBamFile1) {
    gr_lst <- split(ranges_gr,f=seqnames(ranges_gr)) %>%
      .[lapply(.,length)>0]
    counts_lst <- mclapply(gr_lst,function(gr){
      gr2 <- resize(gr,width=span,fix="center")
      df <- bamProfile(BamFile, gr2, paired.end = "midpoint", mapqual = 15, binsize = stepSize) %>%
        {.@signals} %>%
        do.call(rbind,.) %>%
        as.data.frame
      row.names(df) <- as.character(gr)
      return(df)
    })
    names(counts_lst) <- NULL
    counts_df <- do.call(rbind,counts_lst)
    counts_rolled <- apply(counts_df,1,function(cnts){
      zoo::rollapply(cnts, width = windowSize / stepSize, sum, align = "left")
    })
    sortingVector <- colnames(counts_rolled) %>% gsub(".*\\.","",.) %>%
      GRanges %>% width %>% order %>% rev
    counts_rolled_sorted <- t(counts_rolled)[sortingVector,]
    totalReadCount <- idxstatsBam(BamFile) %>%
      {
        .[.$seqnames %in% paste0("chr", c(1:22)), ]
      } %>%
      {
        sum(as.numeric(.[, 3]))
      }
    cpms_rolled <- counts_rolled_sorted / totalReadCount * 1000000

    return(cpms_rolled)
  }
  TX_CPMS1 <- makeCPMS(txBamFile1)
  IN_CPMS1 <- makeCPMS(inBamFile1)
  TX_CPMS2 <- makeCPMS(txBamFile2)
  IN_CPMS2 <- makeCPMS(inBamFile2)

  TX_CPMS <- Reduce("+",list(TX_CPMS1,TX_CPMS2))/2
  IN_CPMS <- Reduce("+",list(IN_CPMS1,IN_CPMS2))/2
  Lg2FC <- log2(TX_CPMS/(IN_CPMS))

  col_fun = colorRamp2(c(-0.5, 0, 0.5), c("green", "white", "red"))

  hmp <- ComplexHeatmap::Heatmap(Lg2FC,col=col_fun,cluster_columns = F,cluster_rows = F)

  return(hmp)
}

#' Convert EduSeq Data Frame to GRanges Object
#'
#' This function transforms an EduSeq-format data frame into a \code{GRanges} object.
#' It calculates genomic coordinates based on the `bin` column and filters to include
#' only standard human chromosomes (chr1–chr22).
#'
#' @param Eduseqdata_df A data frame containing at least a column named \code{bin}.
#'   Additional columns are preserved in the resulting GRanges object.
#'
#' @return A \code{GRanges} object with chromosome, start, end, and any extra metadata columns.
#'
#' @details The start coordinate is calculated as \code{(bin + 1) * 10000}, and the end
#' coordinate is set to \code{start + 1}. Only chromosomes chr1 through chr22 are retained.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame seqnames seqlevels
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(bin = 0:4, seqnames = "chr1")
#' gr <- makeGrForEduSeq(df)
#' }
makeGrForEduSeq <- function(Eduseqdata_df) {
  # Compute the start positions for the GRanges object.
  # The start position is calculated as (bin + 1) * 10000.
  Eduseqdata_df$start <- (Eduseqdata_df$bin + 1) * 10000

  # Compute the end positions for the GRanges object.
  # The end position is defined as start + 1.
  Eduseqdata_df$end <- Eduseqdata_df$start + 1

  # Create a GRanges object from the dataframe.
  # - `keep.extra.columns = T` ensures additional columns in the dataframe are preserved.
  # - `ignore.strand = T` specifies that strand information is not required.
  standard_chromosomes <- paste0("chr", c(1:22))
  gr <- makeGRangesFromDataFrame(Eduseqdata_df,
                                 keep.extra.columns = TRUE,
                                 ignore.strand = TRUE) %>%
    .[seqnames(.) %in% standard_chromosomes]
  seqlevels(gr) <- standard_chromosomes
  gr
}


#' Compute CPM Tracks and y-axis Limits for Two Samples
#'
#' This function computes CPM (counts per million) values for two treatment and two control BAM files
#' across a specified genomic region, using sliding windows. It returns GRanges objects for each condition
#' along with quantile-based y-axis limits useful for visualization.
#'
#' @param chrom Character. Chromosome to analyze (e.g., "chr3").
#' @param trackStart Integer. Start coordinate of the genomic region.
#' @param trackEnd Integer. End coordinate of the genomic region.
#' @param windowSize Integer. Size of the rolling window for smoothing (default: 25000).
#' @param stepSize Integer. Step size for sliding the window (default: 5000).
#' @param txBamFile1 Character. Path to the first treatment BAM file.
#' @param inBamFile1 Character. Path to the first control/input BAM file.
#' @param txBamFile2 Character. Path to the second treatment BAM file.
#' @param inBamFile2 Character. Path to the second control/input BAM file.
#'
#' @return A list with three components:
#' \itemize{
#'   \item \code{tx_cpms_gr}: A \code{GRanges} object with CPM values for treatment samples.
#'   \item \code{in_cpms_gr}: A \code{GRanges} object with CPM values for control samples.
#'   \item \code{yLimits}: A numeric vector of y-axis limits based on CPM distribution quantiles.
#' }
#'
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom Rsamtools idxstatsBam
#' @importFrom bamsignals bamProfile
#' @importFrom zoo rollapply
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- make2SampleCpmGRanges(chrom = "chr3", trackStart = 1e6, trackEnd = 20e6)
#' tx <- result$tx_cpms_gr
#' in <- result$in_cpms_gr
#' yLims <- result$yLimits
#' }
make2SampleCpmGRanges <- function(chrom="chr3",
                                  trackStart=5000000,
                                  trackEnd=55000000,
                                  windowSize=25000,
                                  stepSize=5000,
                                  txBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                  inBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                  txBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_1000.bam",
                                  inBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_1000.bam"
){
  gr <- GRanges(seqnames=chrom,ranges = IRanges(start=trackStart,end=trackEnd))
  windws <- slidingWindows(gr,windowSize,step=stepSize) %>% .[[1]]
  makeCPMS <- function(BamFile=txBamFile1) {
    counts <- bamProfile(BamFile, gr, paired.end = "midpoint", mapqual = 15, binsize = stepSize)
    counts_rolled <- zoo::rollapply(counts@signals[[1]], width = windowSize / stepSize, sum, align = "left")
    totalReadCount <- idxstatsBam(BamFile) %>%
      {
        .[.$seqnames %in% paste0("chr", c(1:22)), ]
      } %>%
      {
        sum(as.numeric(.[, 3]))
      }
    cpms_rolled <- counts_rolled / totalReadCount * 1000000
    cpms_rolled_gr <- windws[1:length(cpms_rolled)] %>%
      {
        .$cpms <- cpms_rolled
        .
      }
    return(cpms_rolled_gr)
  }
  combine2Cpms <- function(BamFile1,BamFile2){
    CPMS_gr_lst <- lapply(c(BamFile1,BamFile2),makeCPMS)
    gr <- CPMS_gr_lst[[1]]
    mcols(gr) <- NULL
    gr$cpms <- cbind(CPMS_gr_lst[[1]]$cpms,CPMS_gr_lst[[1]]$cpms) %>% rowMeans
    return(gr)
  }

  tx_cpms_gr <- combine2Cpms(txBamFile1,txBamFile2)
  in_cpms_gr <- combine2Cpms(inBamFile1,inBamFile2)

  yLimits <- c(
    quantile(
      c(tx_cpms_gr$cpms, in_cpms_gr$cpms),
      0.0001
    ),
    quantile(
      c(tx_cpms_gr$cpms, in_cpms_gr$cpms),
      0.9999
    )
  ) %>% {.*1.2}

  outputList <- list(
    "tx_cpms_gr"=tx_cpms_gr,
    "in_cpms_gr"=in_cpms_gr,
    "yLimits"=yLimits
  )

  return(outputList)
}

#' Generate a Genome-Wide Log2 Fold Change Plot
#'
#' This function computes the log2 fold change (log2FC) in read coverage between treatment and input BAM files
#' over a specified genomic region, and visualizes it as an area plot using a sliding window approach.
#'
#' @param chrom Character. Chromosome to analyze (e.g., `"chr3"`).
#' @param trackStart Integer. Genomic start coordinate of the region to analyze.
#' @param trackEnd Integer. Genomic end coordinate of the region to analyze.
#' @param windowSize Integer. Size of the rolling window for smoothing signal (default: 25000).
#' @param stepSize Integer. Step size between sliding windows (default: 5000).
#' @param txBamFile Character. Path to the treatment BAM file.
#' @param inBamFile Character. Path to the control/input BAM file.
#' @param yLimits Numeric vector of length 2, defining the y-axis limits (default: c(-0.8, 0.8)).
#' @param autoCalculateYLimits Logical. If TRUE, compute y-axis limits automatically from data (default: TRUE).
#' @param fillColor Character. Fill color for the area plot (default: "#006494").
#' @param tx_scaleFactor Numeric or NA. Optional library size scaling factor for treatment BAM (default: NA).
#' @param in_scaleFactor Numeric or NA. Optional library size scaling factor for input BAM (default: NA).
#'
#' @return A `ggplot2` object showing the log2 fold change across the specified region.
#'
#' @importFrom Rsamtools idxstatsBam
#' @importFrom bamsignals bamProfile
#' @importFrom zoo rollapply
#' @importFrom GenomicRanges GRanges disjoin
#' @importFrom IRanges IRanges
#' @importFrom ggplot2 ggplot geom_area theme_void scale_y_continuous scale_x_continuous theme margin element_text aes
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' makeLog2FcGenomePlot(
#'   chrom = "chr3",
#'   trackStart = 1e6,
#'   trackEnd = 2e7,
#'   yLimits = c(-1, 1),
#'   fillColor = "steelblue"
#' )
#' }
makeLog2FcGenomePlot <- function(chrom="chr3",
                                  trackStart=5000000,
                                  trackEnd=55000000,
                                  windowSize=25000,
                                  stepSize=5000,
                                  txBamFile="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                  inBamFile="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                  yLimits=c(-0.8,0.8),
                                  autoCalculateYLimits=TRUE,
                                  fillColor="#006494",
                                  tx_scaleFactor=NA,
                                  in_scaleFactor=NA){
  gr <- GRanges(seqnames=chrom,ranges = IRanges(start=trackStart,end=trackEnd))
  windws <- slidingWindows(gr,windowSize,step=stepSize) %>% .[[1]]
  makeCPMS <- function(BamFile=txBamFile1,scaleFactor) {
    counts <- bamProfile(BamFile, gr, paired.end = "midpoint", mapqual = 15, binsize = stepSize)
    counts_rolled <- zoo::rollapply(counts@signals[[1]], width = windowSize / stepSize, sum, align = "left")
    if (is.na(scaleFactor)) {
      totalReadCount <- idxstatsBam(BamFile) %>%
        {.[.$seqnames %in% paste0("chr", c(1:22)), ]} %>%
        {sum(as.numeric(.[, 3]))}
      cpms_rolled <- counts_rolled / totalReadCount * 1000000
    } else {
      cpms_rolled <- counts_rolled / scaleFactor
    }
    cpms_rolled_gr <- windws[1:length(cpms_rolled)] %>%
      {.$cpms <- cpms_rolled;.}
    return(cpms_rolled_gr)
  }

  tx_cpms_gr <- makeCPMS(txBamFile,tx_scaleFactor)
  in_cpms_gr <- makeCPMS(inBamFile,in_scaleFactor)

  makeDataFrameAndYLimits <- function(txBamFile,inBamFile, txScaleFactor,inScaleFactor){
      tx_cpms_gr <- makeCPMS(txBamFile,txScaleFactor)
      in_cpms_gr <- makeCPMS(inBamFile,inScaleFactor)

      makeLg2FcGr <- function(gr_tx,gr_in){
      gr <- disjoin(gr_tx)
      mcols(gr) <- data.frame("score" = log2(gr_tx$cpms/(gr_in$cpms+0.0000000001)))
      return(gr)}

      LgFC_gr <- makeLg2FcGr(tx_cpms_gr,in_cpms_gr)


      yLimits_layer <- c(quantile(LgFC_gr$score,0.0001),quantile(LgFC_gr$score,0.9999)) %>%
      {.*1.2}

      df <- as.data.frame(LgFC_gr) %>%
      {data.frame(x=rowMeans(.[,c(2,3)]),y=.$score)} %>%
      .[is.finite(.$y),]
      return(list(df,yLimits_layer))
    }

  lst <- makeDataFrameAndYLimits(txBamFile,inBamFile,tx_scaleFactor,in_scaleFactor)

  if (autoCalculateYLimits) {
      yLimits <- lst[[2]]
    }

  df <- lst[[1]]

  ggplt <- ggplot(data=df,aes(x=x,y=y)) +
  geom_area(fill=fillColor) +
  scale_y_continuous(limits = yLimits) +
  scale_x_continuous(limits = c(trackStart,trackEnd),expand=c(0,0)) +
  theme_void() +  # Remove all axes, labels, etc.
  theme(legend.position = "none",  # Remove the legend
        plot.margin = margin(0, 0, 0, 0),
        axis.text.y  = element_text(color = "black", size = 8))

  return(ggplt)
}

#' Overlay Log2 Fold Change Genome Plot for Two Sample Pairs
#'
#' This function computes the log2 fold change (log2FC) between treatment and control BAM files
#' for two distinct sample sets across a specified genomic region, and overlays the results in a single
#' genome plot. Each set is visualized as a separate filled area using different colors.
#'
#' @param chrom Character. Chromosome to analyze (default: `"chr3"`).
#' @param trackStart Integer. Genomic start coordinate of the region to analyze (default: 5000000).
#' @param trackEnd Integer. Genomic end coordinate of the region to analyze (default: 55000000).
#' @param windowSize Integer. Size of the rolling window for smoothing signal (default: 25000).
#' @param stepSize Integer. Step size between sliding windows (default: 5000).
#' @param txBamFile_bottom Character. Path to treatment BAM file for the bottom overlay.
#' @param inBamFile_bottom Character. Path to input/control BAM file for the bottom overlay.
#' @param txBamFile_top Character. Path to treatment BAM file for the top overlay.
#' @param inBamFile_top Character. Path to input/control BAM file for the top overlay.
#' @param yLimits Numeric vector. Y-axis limits for the plot (default: c(-0.8, 0.8)).
#' @param autoCalculateYLimits Logical. If TRUE, compute y-axis limits automatically (default: TRUE).
#' @param bottom_fill Character. Fill color for the bottom overlay (default: "#006494").
#' @param top_fill Character. Fill color for the top overlay (default: "#006494").
#' @param tx_bottom_scaleFactor Numeric or NA. Optional scaling factor for bottom treatment BAM file (default: NA).
#' @param in_bottom_scaleFactor Numeric or NA. Optional scaling factor for bottom control BAM file (default: NA).
#' @param tx_top_scaleFactor Numeric or NA. Optional scaling factor for top treatment BAM file (default: NA).
#' @param in_top_scaleFactor Numeric or NA. Optional scaling factor for top control BAM file (default: NA).
#'
#' @return A `ggplot2` object showing overlaid log2 fold change plots for two sample sets.
#'
#' @importFrom Rsamtools idxstatsBam
#' @importFrom bamsignals bamProfile
#' @importFrom zoo rollapply
#' @importFrom GenomicRanges GRanges disjoin
#' @importFrom IRanges IRanges
#' @importFrom ggplot2 ggplot geom_area scale_fill_manual scale_y_continuous scale_x_continuous theme_void theme margin element_text aes
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' makeLog2FcGenomePlotOverlay(
#'   chrom = "chr3",
#'   trackStart = 1e6,
#'   trackEnd = 2e7,
#'   bottom_fill = "steelblue",
#'   top_fill = "darkorange"
#' )
#' }
makeLog2FcGenomePlotOverlay <- function(chrom="chr3",
                                  trackStart=5000000,
                                  trackEnd=55000000,
                                  windowSize=25000,
                                  stepSize=5000,
                                  txBamFile_bottom="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                  inBamFile_bottom="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                  txBamFile_top="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                  inBamFile_top="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                  yLimits=c(-0.8,0.8),
                                  autoCalculateYLimits=TRUE,
                                  bottom_fill="#006494",
                                  top_fill="#006494",
                                  tx_bottom_scaleFactor=NA,
                                  in_bottom_scaleFactor=NA,
                                  tx_top_scaleFactor=NA,
                                  in_top_scaleFactor=NA){
  gr <- GRanges(seqnames=chrom,ranges = IRanges(start=trackStart,end=trackEnd))
  windws <- slidingWindows(gr,windowSize,step=stepSize) %>% .[[1]]
  makeCPMS <- function(BamFile=txBamFile1,scaleFactor) {
    counts <- bamProfile(BamFile, gr, paired.end = "midpoint", mapqual = 15, binsize = stepSize)
    counts_rolled <- zoo::rollapply(counts@signals[[1]], width = windowSize / stepSize, sum, align = "left")
    if (is.na(scaleFactor)) {
      totalReadCount <- idxstatsBam(BamFile) %>%
        {.[.$seqnames %in% paste0("chr", c(1:22)), ]} %>%
        {sum(as.numeric(.[, 3]))}
      cpms_rolled <- counts_rolled / totalReadCount * 1000000
    } else {
      cpms_rolled <- counts_rolled / scaleFactor
    }
    cpms_rolled_gr <- windws[1:length(cpms_rolled)] %>%
      {.$cpms <- cpms_rolled;.}
    return(cpms_rolled_gr)
  }

  makeDataFrameAndYLimits <- function(txBamFile,inBamFile, txScaleFactor,inScaleFactor){
    tx_cpms_gr <- makeCPMS(txBamFile,txScaleFactor)
    in_cpms_gr <- makeCPMS(inBamFile,inScaleFactor)

    makeLg2FcGr <- function(gr_tx,gr_in){
    gr <- disjoin(gr_tx)
    mcols(gr) <- data.frame("score" = log2(gr_tx$cpms/(gr_in$cpms+0.0000000001)))
    return(gr)}

    LgFC_gr <- makeLg2FcGr(tx_cpms_gr,in_cpms_gr)


    yLimits_layer <- c(quantile(LgFC_gr$score,0.0001),quantile(LgFC_gr$score,0.9999)) %>%
    {.*1.2}

    df <- as.data.frame(LgFC_gr) %>%
    {data.frame(x=rowMeans(.[,c(2,3)]),y=.$score)} %>%
    .[is.finite(.$y),]
    return(list(df,yLimits_layer))
  }

  lst_bottom <- makeDataFrameAndYLimits(txBamFile_bottom,inBamFile_bottom,tx_bottom_scaleFactor,in_bottom_scaleFactor)
  lst_top <- makeDataFrameAndYLimits(txBamFile_top,inBamFile_top,tx_top_scaleFactor,in_top_scaleFactor)

  if (autoCalculateYLimits) {
      yLimits <- c(min(c(lst_bottom[[2]][1],lst_top[[2]][1])),
        max(c(lst_bottom[[2]][2],lst_top[[2]][2])))
    }

  df_bottom <- lst_bottom[[1]] %>%
    {.$layer <- "bottom";.}
  df_top <- lst_top[[1]] %>%
    {.$layer <- "top";.}
  df <- rbind(df_bottom,df_top)

  ggplt <- ggplot(data=df,aes(x=x,y=y,fill=layer)) +
  geom_area() +
  scale_fill_manual(values = c(bottom_fill,top_fill)) +
  scale_y_continuous(limits = yLimits,breaks=c(yLimits[1],0,yLimits[2])) +
  scale_x_continuous(limits = c(trackStart,trackEnd),expand=c(0,0)) +
  theme_void() +  # Remove all axes, labels, etc.
  theme(legend.position = "none",  # Remove the legend
        plot.margin = margin(0, 0, 0, 0),
        axis.text.y  = element_text(color = "black", size = 8))

  return(ggplt)
}

#' Plot Log2 Fold Change Between Two Sample Pairs Using Mean CPM
#'
#' This function computes the average CPM (counts per million) across two treatment and two control BAM files
#' within a specified genomic region, calculates the log2 fold change (log2FC), and plots the result.
#' The signal is smoothed using a rolling window and visualized as a filled area plot.
#'
#' @param chrom Character. Chromosome to analyze (default: `"chr3"`).
#' @param trackStart Integer. Start coordinate of the genomic region (default: 5000000).
#' @param trackEnd Integer. End coordinate of the genomic region (default: 55000000).
#' @param windowSize Integer. Size of the window used for rolling sum (default: 25000).
#' @param stepSize Integer. Step size between sliding windows (default: 5000).
#' @param txBamFile1 Character. Path to the first treatment BAM file.
#' @param inBamFile1 Character. Path to the first input/control BAM file.
#' @param txBamFile2 Character. Path to the second treatment BAM file.
#' @param inBamFile2 Character. Path to the second input/control BAM file.
#' @param yLimits Numeric vector. Limits for the y-axis (default: c(-0.8, 0.8)).
#' @param fillColor Character. Fill color for the area plot (default: "#006494").
#'
#' @return A `ggplot2` object visualizing the log2 fold change between average treatment and control CPMs
#'         across the genomic region.
#'
#' @importFrom Rsamtools idxstatsBam
#' @importFrom bamsignals bamProfile
#' @importFrom zoo rollapply
#' @importFrom GenomicRanges GRanges disjoin
#' @importFrom IRanges IRanges
#' @importFrom ggplot2 ggplot geom_area scale_y_continuous scale_x_continuous theme_void theme margin aes
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' makeLog2FcGenomePlot_2SampleMean(
#'   chrom = "chr3",
#'   trackStart = 1e6,
#'   trackEnd = 2e7,
#'   txBamFile1 = "treatment_250.bam",
#'   inBamFile1 = "input_250.bam",
#'   txBamFile2 = "treatment_1000.bam",
#'   inBamFile2 = "input_1000.bam"
#' )
#' }
makeLog2FcGenomePlot_2SampleMean <- function(chrom="chr3",
                                  trackStart=5000000,
                                  trackEnd=55000000,
                                  windowSize=25000,
                                  stepSize=5000,
                                  txBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                  inBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                  txBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_1000.bam",
                                  inBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_1000.bam",
                                  yLimits=c(-0.8,0.8),
                                  fillColor="#006494"
){
  gr <- GRanges(seqnames=chrom,ranges = IRanges(start=trackStart,end=trackEnd))
  windws <- slidingWindows(gr,windowSize,step=stepSize) %>% .[[1]]
  makeCPMS <- function(BamFile=txBamFile1) {
    counts <- bamProfile(BamFile, gr, paired.end = "midpoint", mapqual = 15, binsize = stepSize)
    counts_rolled <- zoo::rollapply(counts@signals[[1]], width = windowSize / stepSize, sum, align = "left")
    totalReadCount <- idxstatsBam(BamFile) %>%
      {
        .[.$seqnames %in% paste0("chr", c(1:22)), ]
      } %>%
      {
        sum(as.numeric(.[, 3]))
      }
    cpms_rolled <- counts_rolled / totalReadCount * 1000000
    cpms_rolled_gr <- windws[1:length(cpms_rolled)] %>%
      {
        .$cpms <- cpms_rolled
        .
      }
    return(cpms_rolled_gr)
  }
  combine2Cpms <- function(BamFile1,BamFile2){
    CPMS_gr_lst <- lapply(c(BamFile1,BamFile2),makeCPMS)
    gr <- CPMS_gr_lst[[1]]
    mcols(gr) <- NULL
    gr$cpms <- cbind(CPMS_gr_lst[[1]]$cpms,CPMS_gr_lst[[1]]$cpms) %>% rowMeans
    return(gr)
  }

  tx_cpms_gr <- combine2Cpms(txBamFile1,txBamFile2)
  in_cpms_gr <- combine2Cpms(inBamFile1,inBamFile2)

  yLimits <- c(
    quantile(
      c(tx_cpms_gr$cpms, in_cpms_gr$cpms),
      0.0001
    ),
    quantile(
      c(tx_cpms_gr$cpms, in_cpms_gr$cpms),
      0.9999
    )
  ) %>% {.*1.2}

  outputList <- list(
    "tx_cpms_gr"=tx_cpms_gr,
    "in_cpms_gr"=in_cpms_gr,
    "yLimits"=yLimits
  )

  makeLg2FcGr <- function(gr_tx,gr_in){
  gr <- disjoin(gr_tx)
  mcols(gr) <- data.frame(
    "score" = log2(gr_tx$cpms/(gr_in$cpms+0.0000000001))
    )
  return(gr)}

  LgFC_gr <- makeLg2FcGr(outputList[[1]],outputList[[2]])

  df <- as.data.frame(LgFC_gr) %>%
  {data.frame(x=rowMeans(.[,c(2,3)]),y=.$score)} %>%
  .[is.finite(.$y),]

  ggplt <- ggplot(data=df,aes(x=x,y=y)) +
  geom_area(fill=fillColor) +
  scale_y_continuous(limits = yLimits) +
  scale_x_continuous(limits = c(min(RepSeqData_subset$X), max(RepSeqData_subset$X)),expand=c(0,0)) +
  theme_void() +  # Remove all axes, labels, etc.
  theme(legend.position = "none",  # Remove the legend
        plot.margin = margin(0, 0, 0, 0))

  return(ggplt)
}


# the function below does not work well
makeAggregatePlotWithRibbons <- function(
  regions,
  useFactor=FALSE,
  mergedBamsTx=list(
    "siControl"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/ReplicatePeakAnalyzer/results/mergedDownSampledBams/controlMTBP_tx.bam",0.4140772)
      ),
    "siTreslin"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/ReplicatePeakAnalyzer/results/mergedDownSampledBams/siTicrrMTBP_tx.bam",0.3118385)
      )
    ),
  mergedBamsIn=list(
    "siControl"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/ReplicatePeakAnalyzer/results/mergedDownSampledBams/controlMTBP_in.bam",0.2851148)
      ),
    "siTreslin"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/ReplicatePeakAnalyzer/results/mergedDownSampledBams/siTicrrMTBP_in.bam",0.2851148)
      )
    ),
  replicateBamsTx=list(
    "siControl"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/MTBPctrl1.bam",0.36982481),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/MTBPctrl2.bam",0.59287014),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/MTBPctrl3.bam",0.27953668)
      ),
    "siTreslin"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/MTBPsiT1.bam",0.37347100),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/MTBPsiT2.bam",0.30103036),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/MTBPsiT3.bam",0.26101428)
      )
    ),
  replicateBamsIn=list(
    "siControl"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/HCT1.bam",0.39012439),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/HCT2.bam",0.25055516),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/HCT3.bam",0.21466478)
      ),
    "siTreslin"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/HCT1.bam",0.39012439),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/HCT2.bam",0.25055516),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/HCT3.bam",0.21466478)
      )
    )
  )
  {
  makeNormalizedMatrix <- function(
    bamPath = mergedBamsTx[[1]],
    Ranges = regions,
    plottingRange = 4001,
    windowSize = 150,
    stepSize = 5,
    coreCount = 8,
    chromosomes = paste0("chr",1:22)){
      countMatrix <- mclapply(paste0("chr", 1:22), function(chrm){
        chrom_ranges <- Ranges[seqnames(Ranges) == chrm]
        result <- bamProfile(bamPath[[1]], chrom_ranges, binsize = stepSize, paired.end = "midpoint", verbose = TRUE) %>%
          alignSignals(.) %>%
          as.matrix(.) %>%
          rollsumr(., k = windowSize / stepSize, align = "center") + 0.0001
        return(result)
      },
      mc.cores = coreCount) %>%
      do.call(cbind, .)

      # Helper function to calculate total read count
      getTotalReadCount <- function(bam){
        idxstatsBam(bam, index = bam) %>%
          {sapply(chromosomes, function(chr){
            .[which(.[,1] == chr), 3]
          })} %>%
          sum / 1000000
      }

      totalCount <- getTotalReadCount(bamPath[[1]])

      if (useFactor){
        normalizedCounts <- countMatrix / bamPath[[2]]
        } else {
          normalizedCounts <- countMatrix / totalCount
        }

      meanNormalizedCounts <- apply(normalizedCounts,1,mean,na.rm=TRUE)
      return(meanNormalizedCounts)
    }

  # Generate matrices of normalized counts
  mergedBamsIn_normCounts_lst <- lapply(mergedBamsIn, function(bamList) {
    lapply(bamList, makeNormalizedMatrix)
  })

  mergedBamsTx_normCounts_lst <- lapply(mergedBamsTx, function(bamList) {
    lapply(bamList, makeNormalizedMatrix)
  })

  replicateBamsIn_normCounts_lst <- lapply(replicateBamsIn, function(bamList) {
    lapply(bamList, makeNormalizedMatrix)
  })

  replicateBamsTx_normCounts_lst <- lapply(replicateBamsTx, function(bamList) {
    lapply(bamList, makeNormalizedMatrix)
  })

  # Define groups
  groups <- names(mergedBamsTx)

  makeNormalizedCountDataframe <- function(group){
    df <- data.frame("x"=NA,"y"=NA,"group"=NA,"min"=NA,"max"=NA)
    df$x <- 1:length(mergedBamsTx_normCounts_lst[[group]][[1]])
    df$y <- log2(mergedBamsTx_normCounts_lst[[group]][[1]]/(mergedBamsIn_normCounts_lst[[group]][[1]]+0.0000000001))
    replicateInMx <- replicateBamsIn_normCounts_lst[[group]] %>%
      do.call(rbind,.) %>%
      as.matrix
    replicateTxMx <- replicateBamsTx_normCounts_lst[[group]] %>%
      do.call(rbind,.) %>%
      as.matrix
    replicateMx <- log2(replicateTxMx/(replicateInMx+0.0000000001))
    df$min <- apply(replicateMx,2,min)
    df$max <- apply(replicateMx,2,max)
    df$group <- group
    return(df)
  }

  # Assemble data frames
  df <- lapply(groups,makeNormalizedCountDataframe) %>%
    do.call(rbind,.)

  ggplot(data=df, aes(x=x, y=y, color=group, ymin=min, ymax=max)) +
  geom_line() +
  theme_minimal() +
  #scale_color_manual(values = c("siControl" = siControl_color, "siTreslin" = MTBP_color)) +
  scale_x_continuous(
    breaks=c(1,386,772),
    labels = c("-2kb","Summit","+2kb"), expand=c(0,0)) +
  geom_ribbon(fill = "grey70")
}


#' Generate an UpSet plot of reproducible and replicate peaks
#'
#' This function creates an UpSet plot showing the intersection between a set of reproducible peaks
#' and multiple replicate peak sets. It pads replicate peak lists with dummy values to ensure that
#' the set size bars in the UpSet plot accurately reflect the total number of peaks in each replicate.
#' The final plot includes only reproducible peaks in the intersection matrix but uses the accurate set size bars.
#'
#' @param reproduciblePeaksPath Character string. File path to a BED file of reproducible peaks.
#' @param replicatePeakPaths Character vector. File paths to BED files containing peak calls for each replicate.
#' @param replicateLabels Character vector. Display labels for the replicate peak sets (must match order of paths).
#'
#' @return A list object generated by \code{UpSetR::upset()}, representing the UpSet plot.
#' @examples
#' makeUpSetPlot(
#'   reproduciblePeaksPath = "Reproducible_MTBP.bed",
#'   replicatePeakPaths = c("rep1.bed", "rep2.bed", "rep3.bed"),
#'   replicateLabels = c("Rep 1", "Rep 2", "Rep 3")
#' )
#' @import UpSetR
#' @export
makeUpSetPlot <- function(
    reproduciblePeaksPath = "../MakeReproduciblePeaks/MTBP_Peaks_Reproducible_6.narrowPeak",
    replicatePeakPaths = paste0("../01_Asynchronous_HCT116/results/macs2_normalPeaks/", asynchronousPeakFiles[grep("MTBP", asynchronousPeakFiles)]),
    replicateLabels = c("1:1000 1", "1:1000 2", "1:1000 3", "1:250 1", "1:250 2", "1:250 3")
) {
  # Read reproducible peaks into a GenomicRanges object
  reproduciblePeaks <- readBed(reproduciblePeaksPath)

  # Read replicate peak files into a named list of GenomicRanges
  replicatePeaksList <- lapply(replicatePeakPaths, readBed) %>%
    setNames(., replicateLabels)

  # Function to extract overlapping peaks between reproducible peaks and each replicate
  getOverlappingPeaks <- function(queryPeaks, peakList) {
    lapply(peakList, function(gr) {
      queryPeaks[queryPeaks %over% gr] %>%
        { paste0(seqnames(.), ":", start(.), "-", end(.)) }
    })
  }

  # Generate named list of reproducible peaks overlapping each replicate
  overlappingPeaks <- getOverlappingPeaks(reproduciblePeaks, replicatePeaksList) %>%
    setNames(., replicateLabels)

  # Function to pad replicate peak lists with dummy peaks so UpSetR plots total set sizes correctly
  padWithDummyPeaks <- function(sampleLabel, overlapList, allPeaksList) {
    totalPeaks <- length(allPeaksList[[sampleLabel]])
    overlapping <- length(overlapList[[sampleLabel]])
    nMissing <- totalPeaks - overlapping
    dummyPeaks <- paste0(sampleLabel, seq_len(nMissing))
    c(overlapList[[sampleLabel]], dummyPeaks)
  }

  # Pad overlapping peaks with dummy values to match full peak counts per replicate
  paddedPeaks <- lapply(names(overlappingPeaks), padWithDummyPeaks, overlappingPeaks, replicatePeaksList) %>%
    set_names(., names(overlappingPeaks))

  # Create list with only reproducible peaks for UpSetR input (for overlap-only plot)
  overlapOnlyList <- list(paste0(seqnames(reproduciblePeaks), ":", start(reproduciblePeaks), "-", end(reproduciblePeaks))) %>%
    append(overlappingPeaks, .)
  names(overlapOnlyList)[length(overlapOnlyList)] <- "Reproducible Peaks"

  # Create list with dummy-padded peak overlaps (for full peak count plot)
  paddedOverlapList <- list(paste0(seqnames(reproduciblePeaks), ":", start(reproduciblePeaks), "-", end(reproduciblePeaks))) %>%
    append(paddedPeaks, .)
  names(paddedOverlapList)[length(paddedOverlapList)] <- "Reproducible Peaks"

  # Convert both lists into UpSetR-compatible data frames
  overlapOnly_df <- UpSetR::fromList(overlapOnlyList)
  paddedOverlap_df <- UpSetR::fromList(paddedOverlapList)

  # Create UpSetR plot showing only reproducible peak overlaps
  overlapOnlyPlot <- UpSetR::upset(
    overlapOnly_df,
    sets = c(replicateLabels, "Reproducible Peaks"),
    keep.order = TRUE,
    order.by = "freq",
    number.angles = 90,
    nintersects = 20
  )

  # Create UpSetR plot including dummy peaks (for accurate set sizes)
  fullSetSizePlot <- UpSetR::upset(
    paddedOverlap_df,
    group.by = "degree",
    order.by = "freq",
    sets = c(replicateLabels, "Reproducible Peaks"),
    keep.order = TRUE,
    number.angles = 90,
    nintersects = 20
  )

  # Replace set size bars in overlap-only plot with accurate ones from the full set size plot
  overlapOnlyPlot$Sizes <- fullSetSizePlot$Sizes

  # Return the final plot
  overlapOnlyPlot
}


#' Identify reproducible peaks across multiple replicates
#'
#' This function identifies peaks from a query set that overlap with a minimum number of replicates.
#' It returns only those peaks from the query set that are present in at least \code{overlapCountThreshold}
#' of the GRanges objects in \code{subjectPeaksList}.
#'
#' @param queryPeaks A \code{GRanges} object representing the set of query peaks (e.g., reproducible peaks).
#' @param subjectPeaksList A list of \code{GRanges} objects, typically representing peak calls from multiple replicates.
#' @param overlapCountThreshold Integer. Minimum number of replicates in which a peak must be present to be considered reproducible. Default is 5.
#'
#' @return A \code{GRanges} object containing only those peaks from \code{queryPeaks} that overlap at least \code{overlapCountThreshold} replicates.
#' @examples
#' # Example usage:
#' # reproduciblePeaks <- GetReproduciblePeaks(queryPeaks, listOfReplicatePeaks, overlapCountThreshold = 4)
#'
#' @import GenomicRanges
#' @export
GetReproduciblePeaks <- function(queryPeaks,subjectPeaksList,overlapCountThreshold=5){
  lapply(subjectPeaksList,function(gr){
    queryPeaks %over% gr %>%
      as.numeric
  }) %>%
    do.call(cbind,.) %>%
    {rowSums(.)>=overlapCountThreshold} %>%
    queryPeaks[.]
}


#' Generate a data frame of overlap counts across replicate peak sets
#'
#' This function calculates how many peaks from a reproducible peak set are supported by at least
#' a given number of replicate peak sets. It returns a data frame where each row corresponds to
#' a different overlap threshold.
#'
#' @param reproduciblePeaksPath Character string. File path to a BED file containing reproducible peaks.
#' @param replicatePeakPaths Character vector. File paths to BED files of replicate peak calls.
#' @param targetName Character string. Label to associate with the resulting data (e.g., \"MTBP\", \"TRESLIN\").
#'
#' @return A data frame with columns: \code{peakCount}, \code{peakOverlapThreshold}, \code{Target}, and \code{peakPercentage}.
#'
#' @examples
#' makeOverlapThresholdSummary(
#'   reproduciblePeaksPath = "MTBP_Peaks_Reproducible_6.narrowPeak",
#'   replicatePeakPaths = c("rep1.bed", "rep2.bed", "rep3.bed"),
#'   targetName = "MTBP"
#' )
#'
#' @export
makeOverlapThresholdSummary <- function(
  reproduciblePeaksPath = "../MakeReproduciblePeaks/MTBP_Peaks_Reproducible_6.narrowPeak",
  replicatePeakPaths = paste0("../01_Asynchronous_HCT116/results/macs2_normalPeaks/", asynchronousPeakFiles[grep("MTBP", asynchronousPeakFiles)]),
  targetName = "MTBP"
) {
  # Load reproducible peaks into a GRanges object
  reproduciblePeaks <- readBed(reproduciblePeaksPath)

  # Load replicate peak sets into a list of GRanges objects
  replicatePeaksList <- lapply(replicatePeakPaths, readBed)

  # Count how many peaks from the reproducible set are found in at least i replicates (for i = 1 to N)
  overlapCountDf <- lapply(seq_along(replicatePeaksList), function(i) {
    GetReproduciblePeaks(reproduciblePeaks, replicatePeaksList, overlapCountThreshold = i) %>%
      length() %>%
      data.frame(peakCount = ., peakOverlapThreshold = i)
  }) %>%
    do.call(rbind, .) %>%
    { .$Target <- targetName; . }

  overlapCountDf$peakPercentage <- overlapCountDf$peakCount/length(reproduciblePeaks)

  return(overlapCountDf)
}


#' Generate a Correlation Matrix of log2 Fold Change in Read Counts
#'
#' This function calculates the log2 fold change (log2FC) of read counts
#' across a set of genomic regions and visualizes the sample-to-sample
#' correlation matrix using `corrplot`.
#'
#' @param peakRegionsFile Path to a BED file containing the genomic regions of interest.
#' @param treatmentBamPaths Character vector of BAM file paths for treatment samples.
#' @param inputBamPaths Character vector of BAM file paths for input/control samples.
#' @param sampleLabels Character vector of sample names to assign to columns in the output matrix.
#' @param correlationMethod Character string specifying the method used in `cor()` ("spearman", "pearson", etc.). Default is "spearman".
#'
#' @return A correlation matrix plot is displayed via `corrplot`.
#' @importFrom bamsignals bamCount
#' @importFrom corrplot corrplot
#' @export
#'
#' @examples
#' generateCorrelationMatrixOfReadLg2FcInRanges(
#'   peakRegionsFile = "path/to/peaks.bed",
#'   treatmentBamPaths = c("tx_rep1.bam", "tx_rep2.bam"),
#'   inputBamPaths = c("input_rep1.bam", "input_rep2.bam"),
#'   sampleLabels = c("rep1", "rep2")
#' )
generateCorrelationMatrixOfReadLg2FcInRanges <- function(
    peakRegionsFile = "../MakeReproduciblePeaks/MTBP_Peaks_Reproducible_6.narrowPeak",
    treatmentBamPaths,
    inputBamPaths,
    sampleLabels,
    correlationMethod = "spearman"
) {
  # Load peak regions as GRanges
  peakRegions <- readBed(peakRegionsFile)

  # Count reads in peaks for treatment and input samples
  treatmentCounts <- lapply(treatmentBamPaths, bamsignals::bamCount, gr = peakRegions, paired.end = c("midpoint"))
  inputCounts <- lapply(inputBamPaths, bamsignals::bamCount, gr = peakRegions, paired.end = c("midpoint"))

  # Combine counts into matrices
  treatmentMatrix <- do.call(cbind, treatmentCounts)
  inputMatrix <- do.call(cbind, inputCounts)

  # Compute log2 fold change: log2(Tx / Input)
  log2FCMatrix <- log2(treatmentMatrix / (inputMatrix + 1e-7))

  # Label columns for easier interpretation
  colnames(log2FCMatrix) <- sampleLabels

  # Scale values
  log2FCMatrix <- scale(log2FCMatrix)

  # Compute correlation matrix
  correlationMatrix <- cor(log2FCMatrix, method = correlationMethod)

  # Visualize correlation matrix
  list(corrplot::corrplot(correlationMatrix, method = "number"),
    correlationMatrix,
    log2FCMatrix)
}

#' Generate a Library-Normalized Read Count Matrix for Genomic Regions
#'
#' This function computes a matrix of read counts across sliding windows
#' centered on specified genomic regions, normalized by library size (CPM).
#' The output is suitable for metaplot or heatmap visualization.
#'
#' @param bamPath Character string. Path to the BAM file.
#' @param Ranges A \code{GRanges} object representing the genomic regions of interest.
#' @param plottingRange Integer. Total width of the region around each peak (default: 4001).
#' @param windowSize Integer. Size of the smoothing window (default: 150).
#' @param stepSize Integer. Bin size for signal extraction (default: 5).
#' @param coreCount Integer. Number of parallel threads to use (default: 8).
#'
#' @return A matrix of CPM-normalized read counts (rows = positions, columns = regions).
#'
#' @importFrom Rsamtools idxstatsBam
#' @importFrom bamsignals bamProfile
#' @importFrom zoo rollsumr
#' @importFrom parallel mclapply
#' @import GenomicRanges
#' @export
#'
#' @examples
#' peaks <- GRanges("chr1", IRanges(start = c(100000, 200000), width = 4001))
#' normMatrix <- makeLibraryNormalizedMatrix("example.bam", peaks)
makeLibraryNormalizedMatrix <- function(
  bamPath = "../01_Asynchronous_HCT116/results/mergedDownSampledBams/MTBP_WT_HCT116_anti-GFP_250.bam",
  Ranges,
  plottingRange = 4001,
  windowSize = 150,
  stepSize = 5,
  coreCount = 8
) {
  # Define chromosomes to analyze
  chromosomes <- paste0("chr", 1:22)

  # Resize regions to the desired plotting range
  Ranges <- resize(Ranges, width = plottingRange, fix = "center")

  # Compute raw count matrix by chromosome
  countMatrix <- mclapply(chromosomes, function(chr) {
    chrRanges <- Ranges[seqnames(Ranges) == chr]
    signalMatrix <- bamProfile(
      bamPath,
      chrRanges,
      binsize = stepSize,
      paired.end = "midpoint",
      verbose = FALSE
    ) %>%
      alignSignals() %>%
      as.matrix() %>%
      zoo::rollsumr(k = windowSize / stepSize, align = "center") + 1e-4
    return(signalMatrix)
  }, mc.cores = coreCount) %>%
    do.call(cbind, .)

  # Get total aligned reads for library size normalization
  getTotalReadCountInMillions <- function(bamFile) {
    idxstatsBam(bamFile) %>%
      subset(seqnames %in% chromosomes) %>%
      with(sum(as.numeric(mapped))) / 1e6
  }

  totalReads <- getTotalReadCountInMillions(bamPath)

  # Normalize by total reads (CPM)
  normalizedMatrix <- countMatrix / totalReads
  return(normalizedMatrix)
}

#' Import a BED file as a GenomicRanges object
#'
#' This function reads a BED file and converts it to a \code{GRanges} object,
#' retaining only the specified chromosomes.
#'
#' @param bedPath Character string. Path to the BED file.
#' @param chromosomesToImport Character vector. A list of chromosomes to keep (default is \code{chromosomes}, assumed to be defined in the environment).
#'
#' @return A \code{GRanges} object containing only the specified chromosomes.
#'
#' @examples
#' importBED("path/to/peaks.bed", chromosomesToImport = paste0("chr", 1:22))
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb `seqlevels<-`
#' @importFrom GenomeInfoDb seqlevels
#' @export
importBED <- function(bedPath, chromosomesToImport = chromosomes) {
  read.table(bedPath) %>%
    set_names(c("seqnames", "start", "end")) %>%
    GenomicRanges::makeGRangesFromDataFrame() %>%
    { GenomeInfoDb::`seqlevels`(., pruning.mode = "coarse") <- chromosomesToImport; . }
}

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
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges `%over%`
#' @importFrom eulerr euler
#' @importFrom eulerr plot
#'
#' @export
plot_two_sample_euler <- function(gr1, gr2, labels = c("Set1", "Set2"), fills = c("grey", "white")) {
  # Ensure GRanges inputs
  stopifnot(inherits(gr1, "GRanges"), inherits(gr2, "GRanges"))

  # Create union of input GRanges and reduce overlaps
  reduced_peaks_gr <- reduce(c(gr1, gr2))

  # Create binary overlap matrix
  attendance_test_df <- data.frame(
    !!labels[1] := reduced_peaks_gr %over% gr1,
    !!labels[2] := reduced_peaks_gr %over% gr2
  )

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

