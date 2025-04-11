
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
