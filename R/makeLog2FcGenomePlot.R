#' Generate a Genome-Wide Log2 Fold Change Plot
#'
#' This function computes the log2 fold change (log2FC) in read coverage between treatment and input BAM files
#' over a specified genomic region, and visualizes it as an area plot using a sliding window approach.
#'
#' @param chrom Character. Chromosome to analyze (e.g., `"chr3"`).
#' @param trackStart Integer. Genomic start coordinate of the region to analyze.
#' @param trackEnd Integer. Genomic end coordinate of the region to analyze.
#' @param windowSize Integer. Size of the rolling window for smoothing signal (default: 25000).
#' @param stepSize Integer. Step size between sliding windows (default: 5000).
#' @param txBamFile Character. Path to the treatment BAM file.
#' @param inBamFile Character. Path to the control/input BAM file.
#' @param yLimits Numeric vector of length 2, defining the y-axis limits (default: c(-0.8, 0.8)).
#' @param autoCalculateYLimits Logical. If TRUE, compute y-axis limits automatically from data (default: TRUE).
#' @param fillColor Character. Fill color for the area plot (default: "#006494").
#' @param tx_scaleFactor Numeric or NA. Optional library size scaling factor for treatment BAM (default: NULL).
#' @param in_scaleFactor Numeric or NA. Optional library size scaling factor for input BAM (default: NULL).
#' @param blacklist_granges GenomicRanges object with blacklisted regions (default:  NULL)
#'
#' @return A `ggplot2` object showing the log2 fold change across the specified region.
#'
#' @importFrom Rsamtools idxstatsBam
#' @importFrom bamsignals bamProfile
#' @importFrom zoo rollapply
#' @importFrom GenomicRanges GRanges disjoin
#' @importFrom IRanges IRanges
#' @importFrom ggplot2 ggplot geom_area theme_void scale_y_continuous scale_x_continuous theme margin element_text aes
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' makeLog2FcGenomePlot(
#'   chrom = "chr3",
#'   trackStart = 1e6,
#'   trackEnd = 2e7,
#'   yLimits = c(-1, 1),
#'   fillColor = "steelblue"
#' )
#' }
makeLog2FcGenomePlot <- function(chrom="chr3",
                                  trackStart=5000000,
                                  trackEnd=55000000,
                                  windowSize=25000,
                                  stepSize=5000,
                                  txBamFile="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                  inBamFile="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                  yLimits=c(-0.8,0.8),
                                  autoCalculateYLimits=TRUE,
                                  fillColor="#006494",
                                  tx_scaleFactor=NULL,
                                  in_scaleFactor=NULL,
                                 blacklist_granges = NULL){
  gr <- GRanges(seqnames=chrom,ranges = IRanges(start=trackStart,end=trackEnd))
  windws <- slidingWindows(gr,windowSize,step=stepSize) %>% .[[1]]
  makeCPMS <- function(BamFile=txBamFile1,scaleFactor) {
    counts <- bamProfile(BamFile, gr, paired.end = "midpoint", mapqual = 15, binsize = stepSize)
    counts_rolled <- zoo::rollapply(counts@signals[[1]], width = windowSize / stepSize, sum, align = "left")
    if (is.null(scaleFactor)) {
      totalReadCount <- idxstatsBam(BamFile) %>%
        {.[.$seqnames %in% paste0("chr", c(1:22)), ]} %>%
        {sum(as.numeric(.[, 3]))}
      cpms_rolled <- counts_rolled / totalReadCount * 1000000
    } else {
      cpms_rolled <- counts_rolled / scaleFactor
    }
    cpms_rolled_gr <- windws[1:length(cpms_rolled)] %>%
      {.$cpms <- cpms_rolled;.}
    if (is.null(blacklist_granges)) {
      return(cpms_rolled_gr)
    } else {
      return(cpms_rolled_gr[!cpms_rolled_gr %over% blacklist_granges])
    }

  }

  tx_cpms_gr <- makeCPMS(txBamFile,tx_scaleFactor)
  in_cpms_gr <- makeCPMS(inBamFile,in_scaleFactor)

  makeDataFrameAndYLimits <- function(txBamFile,inBamFile, txScaleFactor,inScaleFactor){
      tx_cpms_gr <- makeCPMS(txBamFile,txScaleFactor)
      in_cpms_gr <- makeCPMS(inBamFile,inScaleFactor)

      makeLg2FcGr <- function(gr_tx,gr_in){
      gr <- disjoin(gr_tx)
      mcols(gr) <- data.frame("score" = log2(gr_tx$cpms/(gr_in$cpms+0.0000000001)))
      return(gr)}

      LgFC_gr <- makeLg2FcGr(tx_cpms_gr,in_cpms_gr)


      yLimits_layer <- c(quantile(LgFC_gr$score,0.0001),quantile(LgFC_gr$score,0.9999)) %>%
      {.*1.2}

      df <- as.data.frame(LgFC_gr) %>%
      {data.frame(x=rowMeans(.[,c(2,3)]),y=.$score)} %>%
      .[is.finite(.$y),]
      return(list(df,yLimits_layer))
    }

  lst <- makeDataFrameAndYLimits(txBamFile,inBamFile,tx_scaleFactor,in_scaleFactor)

  if (autoCalculateYLimits) {
      yLimits <- lst[[2]]
    }

  df <- lst[[1]]

  ggplt <- ggplot(data=df,aes(x=x,y=y)) +
  geom_area(fill=fillColor) +
  scale_y_continuous(limits = yLimits) +
  scale_x_continuous(limits = c(trackStart,trackEnd),expand=c(0,0)) +
  theme_void() +  # Remove all axes, labels, etc.
  theme(legend.position = "none",  # Remove the legend
        plot.margin = margin(0, 0, 0, 0),
        axis.text.y  = element_text(color = "black", size = 8))

  return(ggplt)
}
