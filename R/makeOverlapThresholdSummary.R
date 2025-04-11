
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
