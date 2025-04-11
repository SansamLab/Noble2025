
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
