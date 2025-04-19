#' Convert UCSC-style coordinate strings to a GRanges object
#'
#' Parses UCSC-style genome coordinate strings (e.g., "chr1:10000-20000")
#' and returns a corresponding \code{GRanges} object.
#'
#' @param coord_strings A character vector of UCSC-style coordinates, each in the format "chr:start-end".
#' @param genome A character string indicating the genome build (e.g., "hg38", "mm10"). Default is "unknown".
#'
#' @return A \code{GRanges} object with one range per input string.
#'
#' @examples
#' ucsc_coords <- c("chr1:10000-20000", "chr2:15000-15500")
#' gr <- parseUCSCtoGRanges(ucsc_coords, genome = "hg38")
#' gr
#'
#' @import GenomicRanges
#' @import IRanges
#' @export
#'

parseUCSCtoGRanges <- function(coord_strings, genome = "unknown") {
  split_coords <- lapply(coord_strings, function(x) {
    parts <- strsplit(x, ":|-")[[1]]
    data.frame(
      seqnames = parts[1],
      start = as.numeric(parts[2]),
      end = as.numeric(parts[3]),
      stringsAsFactors = FALSE
    )
  })

  coords_df <- do.call(rbind, split_coords)

  GenomicRanges::GRanges(
    seqnames = coords_df$seqnames,
    ranges = IRanges::IRanges(start = coords_df$start, end = coords_df$end),
    genome = genome
  )
}
