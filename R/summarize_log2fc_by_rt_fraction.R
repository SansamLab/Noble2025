
#' Annotate Log2FC GRanges by Replication Timing Fraction
#'
#' Matches log2 fold change values from a GRanges object to high-resolution replication timing (Repli-seq)
#' fractions, assigning each region to the S-phase fraction where replication signal is maximal.
#'
#' @param log2fc_gr A `GRanges` object with a metadata column named `Log2FC`.
#' @param repliseq_gr A `GRanges` object with Repli-seq data across S-phase fractions. Each metadata column
#'   corresponds to an S-phase fraction (e.g., S1 to S16).
#'
#' @return A `data.frame` with columns:
#'   \describe{
#'     \item{S.fraction}{S-phase fraction with max Repli-seq signal (e.g., "S1", "S2", ...).}
#'     \item{Log2FC}{Log2 fold change value from the input `log2fc_gr`.}
#'   }
#'
#' @import GenomicRanges
#' @import S4Vectors
#' @import GenomeInfoDb
#' @export
summarize_log2fc_by_rt_fraction <- function(log2fc_gr, repliseq_gr) {
  repliseq_clean <- repliseq_gr[complete.cases(mcols(repliseq_gr))] %>%
    sortSeqlevels() %>% sort()

  repliseq_clean$Max <- apply(mcols(repliseq_clean), 1, which.max) %>% paste0("S", .)

  log2fc_gr_filtered <- log2fc_gr[!is.na(match(log2fc_gr, repliseq_clean))] %>%
    sortSeqlevels() %>% sort()

  data.frame(
    "S.fraction" = factor(repliseq_clean$Max, levels = paste0("S", 1:16)),
    "Log2FC" = log2fc_gr_filtered$Log2FC
  )
}
