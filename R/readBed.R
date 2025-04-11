#' Read and convert a BED file to a GenomicRanges object
#'
#' This function reads a BED file and converts it into a GenomicRanges object
#' using the GenomicRanges package. The BED file should have columns named
#' "seqnames", "start", and "end" to represent genomic ranges.
#' @import GenomicRanges
#' @export
#' @param bedFile Path to the BED file to be read and converted.
#'
#' @return A GenomicRanges object representing the data in the BED file.
#'
#' @examples
#' bedFile <- "path/to/your/bedfile.bed"
#' gr <- readBed(bedFile)
readBed <- function(bedFile){
  bedData <- read.table(bedFile)
  if(!is.numeric(bedData[1,2])){
    bedData <- bedData[-1,]
  }
  bedData <- set_names(bedData, c("seqnames", "start", "end"))
  return(makeGRangesFromDataFrame(bedData))
}
