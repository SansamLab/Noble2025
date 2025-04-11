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
#' @export
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
