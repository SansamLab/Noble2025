
#' Calculate Peak Density by Replication Timing Fraction
#'
#' Computes the density of peaks (overlaps) in each replication timing fraction
#' based on overlap with a set of peaks (e.g., ChIP or origin peaks). Peak density
#' is normalized per 50 kb.
#'
#' @param repliseq_gr A \code{GRanges} object containing replication timing windows with 16 S phase fractions.
#' @param peak_gr A \code{GRanges} object containing genomic regions such as peaks.
#'
#' @return A \code{data.frame} with one row per replication timing fraction and the following columns:
#' \itemize{
#'   \item \code{S.fraction}: Fraction label (S1–S16)
#'   \item \code{bp}: Total base pairs in that S fraction
#'   \item \code{peakCount}: Total number of peaks overlapping the S fraction
#'   \item \code{meanPeakCount}: Mean number of peaks per 50kb window in the S fraction
#'   \item \code{peakDensity}: Normalized peak density (peaks per 50 kb)
#'   \item \code{S.fraction.numeric}: Numeric version of the S fraction (1–16)
#' }
#'
#' @import GenomicRanges
#' @import IRanges
#' @export
calculate_peak_density_by_fraction <- function(repliseq_gr, peak_gr) {
  repliseq_clean <- repliseq_gr[complete.cases(mcols(repliseq_gr))] %>%
    sortSeqlevels() %>%
    sort()

  repliseq_clean$OverlapCounts <- countOverlaps(repliseq_clean, peak_gr)
  repliseq_clean$S.fraction <- apply(mcols(repliseq_clean), 1, which.max) %>% paste0("S", .)
  repliseq_clean$bp <- width(repliseq_clean)

  df <- as.data.frame(mcols(repliseq_clean))
  df2 <- aggregate(bp ~ S.fraction, df, sum)
  df2$peakCount <- aggregate(OverlapCounts ~ S.fraction, df, sum)$OverlapCounts
  df2$meanPeakCount <- aggregate(OverlapCounts ~ S.fraction, df, mean)$OverlapCounts
  df2$peakDensity <- df2$peakCount / df2$bp * 50000
  df2$S.fraction.numeric <- gsub("S", "", df2$S.fraction) %>% as.numeric()

  return(df2)
}
