#' Import a BED file as a GenomicRanges object
#'
#' This function reads a BED file and converts it to a \code{GRanges} object,
#' retaining only the specified chromosomes.
#'
#' @param bedPath Character string. Path to the BED file.
#' @param chromosomesToImport Character vector. A list of chromosomes to keep (default is \code{chromosomes}, assumed to be defined in the environment).
#'
#' @return A \code{GRanges} object containing only the specified chromosomes.
#'
#' @examples
#' importBED("path/to/peaks.bed", chromosomesToImport = paste0("chr", 1:22))
#'
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @export
importBED <- function(bedPath, chromosomesToImport = paste0("chr",1:22)) {
  chromToImport <- as.character(chromosomesToImport)

  # Read first line to detect if it contains a header
  first_row <- read.table(bedPath, nrows = 1, stringsAsFactors = FALSE)

  # Detect if columns 2 and 3 (start, end) are numeric
  second_col_is_numeric <- suppressWarnings(!is.na(as.numeric(first_row[[2]])))
  third_col_is_numeric <- suppressWarnings(!is.na(as.numeric(first_row[[3]])))

  has_header <- !(second_col_is_numeric && third_col_is_numeric)

  # Read full table with or without header
  bed_df <- read.table(bedPath, header = has_header, stringsAsFactors = FALSE)

  # Standardize required column names
  colnames(bed_df)[1:3] <- c("seqnames", "start", "end")

  # Convert to GRanges and restrict to specified chromosomes
  gr <- GenomicRanges::makeGRangesFromDataFrame(bed_df, keep.extra.columns = TRUE)
  seqlevels(gr, pruning.mode = "coarse") <- chromToImport

  return(gr)
}
