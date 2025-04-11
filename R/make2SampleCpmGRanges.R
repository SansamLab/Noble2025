#' Compute CPM Tracks and y-axis Limits for Two Samples
#'
#' This function computes CPM (counts per million) values for two treatment and two control BAM files
#' across a specified genomic region, using sliding windows. It returns GRanges objects for each condition
#' along with quantile-based y-axis limits useful for visualization.
#'
#' @param chrom Character. Chromosome to analyze (e.g., "chr3").
#' @param trackStart Integer. Start coordinate of the genomic region.
#' @param trackEnd Integer. End coordinate of the genomic region.
#' @param windowSize Integer. Size of the rolling window for smoothing (default: 25000).
#' @param stepSize Integer. Step size for sliding the window (default: 5000).
#' @param txBamFile1 Character. Path to the first treatment BAM file.
#' @param inBamFile1 Character. Path to the first control/input BAM file.
#' @param txBamFile2 Character. Path to the second treatment BAM file.
#' @param inBamFile2 Character. Path to the second control/input BAM file.
#'
#' @return A list with three components:
#' \itemize{
#'   \item \code{tx_cpms_gr}: A \code{GRanges} object with CPM values for treatment samples.
#'   \item \code{in_cpms_gr}: A \code{GRanges} object with CPM values for control samples.
#'   \item \code{yLimits}: A numeric vector of y-axis limits based on CPM distribution quantiles.
#' }
#'
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom Rsamtools idxstatsBam
#' @importFrom bamsignals bamProfile
#' @importFrom zoo rollapply
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- make2SampleCpmGRanges(chrom = "chr3", trackStart = 1e6, trackEnd = 20e6)
#' tx <- result$tx_cpms_gr
#' in <- result$in_cpms_gr
#' yLims <- result$yLimits
#' }
make2SampleCpmGRanges <- function(chrom="chr3",
                                  trackStart=5000000,
                                  trackEnd=55000000,
                                  windowSize=25000,
                                  stepSize=5000,
                                  txBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                  inBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                  txBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_1000.bam",
                                  inBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_1000.bam"
){
  gr <- GRanges(seqnames=chrom,ranges = IRanges(start=trackStart,end=trackEnd))
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
  combine2Cpms <- function(BamFile1,BamFile2){
    CPMS_gr_lst <- lapply(c(BamFile1,BamFile2),makeCPMS)
    gr <- CPMS_gr_lst[[1]]
    mcols(gr) <- NULL
    gr$cpms <- cbind(CPMS_gr_lst[[1]]$cpms,CPMS_gr_lst[[1]]$cpms) %>% rowMeans
    return(gr)
  }

  tx_cpms_gr <- combine2Cpms(txBamFile1,txBamFile2)
  in_cpms_gr <- combine2Cpms(inBamFile1,inBamFile2)

  yLimits <- c(
    quantile(
      c(tx_cpms_gr$cpms, in_cpms_gr$cpms),
      0.0001
    ),
    quantile(
      c(tx_cpms_gr$cpms, in_cpms_gr$cpms),
      0.9999
    )
  ) %>% {.*1.2}

  outputList <- list(
    "tx_cpms_gr"=tx_cpms_gr,
    "in_cpms_gr"=in_cpms_gr,
    "yLimits"=yLimits
  )

  return(outputList)
}