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
