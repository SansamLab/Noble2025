#' Create a plot from BED files
#'
#' This function reads a list of BED files, calculates overlaps, and creates a plot
#' based on the specified plot choice (UpSet or Euler). It can be used to visualize
#' genomic overlap patterns between the BED files.
#'
#' @param BedFilenames A character vector containing paths to the BED files.
#' @param plotChoice A character string specifying the type of plot to create ("UpSet" or "Euler").
#' @param sampleLabelsDF A data frame with file labels and corresponding file paths.
#' @param colrs A list of colors for the plot elements.
#' @param alfa Alpha transparency value for plot elements.
#' @param fontScale Scaling factor for text size in the plot.
#' @param minFracOverlap Minimum fraction of overlap required for inclusion in the plot.
#'
#' @return NULL (The function generates and displays the plot without returning a value.)
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce findOverlaps width pintersect seqnames start end queryHits subjectHits
#' @import UpSetR
#' @export
#'
#' @examples
#' BedFilenames <- c("path/to/bedfile1.bed", "path/to/bedfile2.bed")
#' plotChoice <- "UpSet"
#' sampleLabelsDF <- data.frame(files = BedFilenames, labels = c("Label1", "Label2"))
#' colrs <- list(colors = c("#1b9e77", "#d95f02"))
#' alfa <- 0.7
#' fontScale <- 8L
#' minFracOverlap <- 0L
#' makePlotFromBeds(BedFilenames, plotChoice, sampleLabelsDF, colrs, alfa, fontScale, minFracOverlap)
#'
makePlotFromBeds <- function(BedFilenames, plotChoice, sampleLabelsDF, colrs, alfa, fontScale, minFracOverlap) {
  # Read and convert BedFilenames to GenomicRanges
  Beds.gr <- lapply(BedFilenames, function(file) {
    bed.df <- read.table(file)
    GenomicRanges::makeGRangesFromDataFrame(bed.df, seqnames.field = "V1", start.field = "V2", end.field = "V3")
  })

  # Reduce to get AllBeds.gr
  AllBeds.gr <- GenomicRanges::reduce(do.call("c", Beds.gr))

  # Function to get overlaps
  GetOverlapsWithAll <- function(bed.gr) {
    hits <- GenomicRanges::findOverlaps(AllBeds.gr, bed.gr)
    sigHits <- hits[GenomicRanges::width(GenomicRanges::pintersect(AllBeds.gr[queryHits(hits)], bed.gr[subjectHits(hits)])) / GenomicRanges::width(AllBeds.gr[queryHits(hits)]) > minFracOverlap]
    GRovrlps <- AllBeds.gr[queryHits(sigHits)]
    paste(GenomicRanges::seqnames(GRovrlps), GenomicRanges::start(GRovrlps), GenomicRanges::end(GRovrlps), sep = "_")
  }

  # Get overlaps for all BedFilenames
  lst <- lapply(Beds.gr, GetOverlapsWithAll)

  # Match sampleLabelsDF and remove file extensions
  nmes <- gsub("\\.[^.]*$", "", sampleLabelsDF$labels[match(sampleLabelsDF$files, basename(BedFilenames))])

  # Set names for lst
  names(lst) <- nmes

  # Create upsetList and plot
  upsetList <- UpSetR::fromList(lst)
  if (plotChoice == "UpSet") {
    UpSetR::upset(upsetList, nsets = length(lst), order.by = "freq", text.scale = fontScale)
  } else {
    plot(eulerr::euler(upsetList), quantities = TRUE, fills = colrs$colors, alpha = alfa)
  }
}
