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
#' @export
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
