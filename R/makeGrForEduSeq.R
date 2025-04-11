#' Convert EduSeq Data Frame to GRanges Object
#'
#' This function transforms an EduSeq-format data frame into a \code{GRanges} object.
#' It calculates genomic coordinates based on the `bin` column and filters to include
#' only standard human chromosomes (chr1–chr22).
#'
#' @param Eduseqdata_df A data frame containing at least a column named \code{bin}.
#'   Additional columns are preserved in the resulting GRanges object.
#'
#' @return A \code{GRanges} object with chromosome, start, end, and any extra metadata columns.
#'
#' @details The start coordinate is calculated as \code{(bin + 1) * 10000}, and the end
#' coordinate is set to \code{start + 1}. Only chromosomes chr1 through chr22 are retained.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame seqnames seqlevels
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(bin = 0:4, seqnames = "chr1")
#' gr <- makeGrForEduSeq(df)
#' }
makeGrForEduSeq <- function(Eduseqdata_df) {
  # Compute the start positions for the GRanges object.
  # The start position is calculated as (bin + 1) * 10000.
  Eduseqdata_df$start <- (Eduseqdata_df$bin + 1) * 10000

  # Compute the end positions for the GRanges object.
  # The end position is defined as start + 1.
  Eduseqdata_df$end <- Eduseqdata_df$start + 1

  # Create a GRanges object from the dataframe.
  # - `keep.extra.columns = T` ensures additional columns in the dataframe are preserved.
  # - `ignore.strand = T` specifies that strand information is not required.
  standard_chromosomes <- paste0("chr", c(1:22))
  gr <- makeGRangesFromDataFrame(Eduseqdata_df,
                                 keep.extra.columns = TRUE,
                                 ignore.strand = TRUE) %>%
    .[seqnames(.) %in% standard_chromosomes]
  seqlevels(gr) <- standard_chromosomes
  gr
}
