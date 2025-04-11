
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
