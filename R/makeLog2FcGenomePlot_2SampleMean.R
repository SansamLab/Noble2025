

#' Plot Log2 Fold Change Between Two Sample Pairs Using Mean CPM
#'
#' This function computes the average CPM (counts per million) across two treatment and two control BAM files
#' within a specified genomic region, calculates the log2 fold change (log2FC), and plots the result.
#' The signal is smoothed using a rolling window and visualized as a filled area plot.
#'
#' @param chrom Character. Chromosome to analyze (default: `"chr3"`).
#' @param trackStart Integer. Start coordinate of the genomic region (default: 5000000).
#' @param trackEnd Integer. End coordinate of the genomic region (default: 55000000).
#' @param windowSize Integer. Size of the window used for rolling sum (default: 25000).
#' @param stepSize Integer. Step size between sliding windows (default: 5000).
#' @param txBamFile1 Character. Path to the first treatment BAM file.
#' @param inBamFile1 Character. Path to the first input/control BAM file.
#' @param txBamFile2 Character. Path to the second treatment BAM file.
#' @param inBamFile2 Character. Path to the second input/control BAM file.
#' @param yLimits Numeric vector. Limits for the y-axis (default: c(-0.8, 0.8)).
#' @param fillColor Character. Fill color for the area plot (default: "#006494").
#'
#' @return A `ggplot2` object visualizing the log2 fold change between average treatment and control CPMs
#'         across the genomic region.
#'
#' @importFrom Rsamtools idxstatsBam
#' @importFrom bamsignals bamProfile
#' @importFrom zoo rollapply
#' @importFrom GenomicRanges GRanges disjoin
#' @importFrom IRanges IRanges
#' @importFrom ggplot2 ggplot geom_area scale_y_continuous scale_x_continuous theme_void theme margin aes
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' makeLog2FcGenomePlot_2SampleMean(
#'   chrom = "chr3",
#'   trackStart = 1e6,
#'   trackEnd = 2e7,
#'   txBamFile1 = "treatment_250.bam",
#'   inBamFile1 = "input_250.bam",
#'   txBamFile2 = "treatment_1000.bam",
#'   inBamFile2 = "input_1000.bam"
#' )
#' }
makeLog2FcGenomePlot_2SampleMean <- function(chrom="chr3",
                                  trackStart=5000000,
                                  trackEnd=55000000,
                                  windowSize=25000,
                                  stepSize=5000,
                                  txBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                  inBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                  txBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_1000.bam",
                                  inBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_1000.bam",
                                  yLimits=c(-0.8,0.8),
                                  fillColor="#006494"
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

  makeLg2FcGr <- function(gr_tx,gr_in){
  gr <- disjoin(gr_tx)
  mcols(gr) <- data.frame(
    "score" = log2(gr_tx$cpms/(gr_in$cpms+0.0000000001))
    )
  return(gr)}

  LgFC_gr <- makeLg2FcGr(outputList[[1]],outputList[[2]])

  df <- as.data.frame(LgFC_gr) %>%
  {data.frame(x=rowMeans(.[,c(2,3)]),y=.$score)} %>%
  .[is.finite(.$y),]

  ggplt <- ggplot(data=df,aes(x=x,y=y)) +
  geom_area(fill=fillColor) +
  scale_y_continuous(limits = yLimits) +
  scale_x_continuous(limits = c(min(RepSeqData_subset$X), max(RepSeqData_subset$X)),expand=c(0,0)) +
  theme_void() +  # Remove all axes, labels, etc.
  theme(legend.position = "none",  # Remove the legend
        plot.margin = margin(0, 0, 0, 0))

  return(ggplt)
}
