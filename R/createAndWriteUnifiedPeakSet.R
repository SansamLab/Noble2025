#' Create a unified peak set from multiple peak files and write it to a file
#'
#' This function takes a list of peak files and creates a unified peak set
#' by combining the second set of peaks with non-overlapping peaks from the first set.
#' The resulting peak set is then written to a specified output file in BED format.
#'
#' @param pkFiles A list of paths to the peak files to be combined.
#' @param outputPath Path to the output file where the unified peak set will be saved.
#'
#' @return NULL (The function writes the result to a file and doesn't return a value.)
#' @export
#' @import GenomicRanges
#' @examples
#' pkFiles <- c("path/to/peakfile1.bed", "path/to/peakfile2.bed")
#' outputPath <- "path/to/output/unified_peaks"
#' createAndWriteUnifiedPeakSet(pkFiles, outputPath)
#'
createAndWriteUnifiedPeakSet <- function(pkFiles, smtFiles,outputPath){
  # Read and convert peak files to GenomicRanges objects
  peaks <- lapply(pkFiles, readBed)
  summits <- lapply(smtFiles,readBed)

  # Combine the second set of peaks with non-overlapping peaks from the first set
  unifiedPeaks <- peaks[[2]][!peaks[[2]] %over% peaks[[1]]] %>% c(., peaks[[1]])
  unifiedSummits <- summits[[2]][!peaks[[2]] %over% peaks[[1]]] %>% c(., summits[[1]])

  # Create the output directory if it doesn't exist
  if (!file.exists(dirname(outputPath))) {
    dir.create(dirname(outputPath), recursive = TRUE)
  }

  # Write the unified peak set to the specified output file
  write.table(unifiedPeaks %>%
                as.data.frame %>%
                .[,1:3],
              paste0(outputPath,".bed"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  write.table(unifiedSummits %>%
                as.data.frame %>%
                .[,1:3],
              paste0(outputPath,"_summits.bed"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  mcols(unifiedPeaks)$summits <- start(unifiedSummits)
  saveRDS(unifiedPeaks,
          paste0(outputPath,".rds"))
}
