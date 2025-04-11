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