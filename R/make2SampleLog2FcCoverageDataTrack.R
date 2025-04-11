
#' Create a DataTrack for 2-Sample Log2 Fold Change (Log2FC) Coverage Data
#'
#' This function generates a DataTrack for 2-sample Log2FC coverage data based on
#' input BAM files. It calculates Log2FC values for specified genomic regions
#' and creates a DataTrack for visualization.
#'
#' @param chromosome The chromosome for the genomic region of interest.
#' @param start The start position of the genomic region.
#' @param end The end position of the genomic region.
#' @param windowSize The size of the sliding windows for Log2FC calculation.
#' @param stepSize The step size for sliding windows.
#' @param txBamFile1 The path to the first treatment BAM file.
#' @param inBamFile1 The path to the first control BAM file.
#' @param txBamFile2 The path to the second treatment BAM file.
#' @param inBamFile2 The path to the second control BAM file.
#' @param trackName The name for the DataTrack.
#' @param HistogramColor The color for the histogram in the DataTrack.
#'
#' @return A DataTrack object for 2-sample Log2FC coverage data.
#'
#' @import GenomicRanges
#' @import zoo
#' @export
#' @examples
#' \dontrun{
#' # Example usage:
#' dt <- make2SampleLog2FcCoverageDataTrack(chromosome="chr3", start=5000000, end=55000000,
#'                                          windowSize=25000, stepSize=5000,
#'                                          txBamFile1="path/to/txBamFile1.bam",
#'                                          inBamFile1="path/to/inBamFile1.bam",
#'                                          txBamFile2="path/to/txBamFile2.bam",
#'                                          inBamFile2="path/to/inBamFile2.bam",
#'                                          trackName="Log2FC", HistogramColor="#1b9e77")
#' }
#'
make2SampleLog2FcCoverageDataTrack <- function(chromosome="chr3",
                                               start=5000000,
                                               end=55000000,
                                               windowSize=25000,
                                               stepSize=5000,
                                               txBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                               inBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                               txBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_1000.bam",
                                               inBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_1000.bam",
                                               trackName="Log2FC",
                                               HistogramColor="#1b9e77",
                                               yLimits=c(-0.5,0.5)
){
  gr <- GRanges(seqnames=chromosome,ranges = IRanges(start=start,end=end))
  windws <- slidingWindows(gr,windowSize,step=stepSize) %>% .[[1]]
  makeCPMS <- function(BamFile=txBamFile1) {
    counts <- bamProfile(BamFile, gr, paired.end = "midpoint", mapqual = 15, binsize = stepSize)
    counts_rolled <- zoo::rollapply(counts@signals[[1]], width = windowSize / stepSize, sum, align = "left")
    totalReadCount <- idxstatsBam(BamFile) %>%
      {
        .[.$seqnames %in% paste0("chr", c(1:22)), ]
      } %>%
      {
        sum(as.numeric(.[, 3]))
      }
    cpms_rolled <- counts_rolled / totalReadCount * 1000000
    cpms_rolled_gr <- windws[1:length(cpms_rolled)] %>%
      {
        .$cpms <- cpms_rolled
        .
      }
    return(cpms_rolled_gr)
  }
  TX_CPMS1 <- makeCPMS()
  IN_CPMS1 <- makeCPMS(inBamFile1)
  TX_CPMS2 <- makeCPMS(txBamFile2)
  IN_CPMS2 <- makeCPMS(inBamFile2)

  TX_CPMS <- cbind(TX_CPMS1$cpms,TX_CPMS2$cpms) %>% rowMeans
  IN_CPMS <- cbind(IN_CPMS1$cpms,IN_CPMS2$cpms) %>% rowMeans

  output_gr <- TX_CPMS1
  mcols(output_gr) <- NULL
  mcols(output_gr) <- log2(TX_CPMS/(IN_CPMS+0.0000000001))
  if(!is.na(yLimits)){
    dt <- DataTrack(
      range = output_gr, genome = "hg38", type = "hist",name = trackName,
      #col = HistogramColor, col.histogram = HistogramColor, fill.histogram = HistogramColor,
      fill.histogram = HistogramColor,col.histogram = "transparent",col="transparent",
      lwd = 0.5, ylim = yLimits
    )
  } else {
    dt <- DataTrack(
      range = output_gr, genome = "hg38", type = "hist",name = trackName,
      #col = HistogramColor, col.histogram = HistogramColor, fill.histogram = HistogramColor,
      fill.histogram = HistogramColor,col.histogram = "transparent",col="transparent",
      lwd = 0.5, ylim=c(quantile(output_gr$X,0.0001), quantile(output_gr$X,0.9999)))
  }

  return(dt)
}
