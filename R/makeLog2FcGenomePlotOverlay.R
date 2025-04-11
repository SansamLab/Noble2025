#' Overlay Log2 Fold Change Genome Plot for Two Sample Pairs
#'
#' This function computes the log2 fold change (log2FC) between treatment and control BAM files
#' for two distinct sample sets across a specified genomic region, and overlays the results in a single
#' genome plot. Each set is visualized as a separate filled area using different colors.
#'
#' @param chrom Character. Chromosome to analyze (default: `"chr3"`).
#' @param trackStart Integer. Genomic start coordinate of the region to analyze (default: 5000000).
#' @param trackEnd Integer. Genomic end coordinate of the region to analyze (default: 55000000).
#' @param windowSize Integer. Size of the rolling window for smoothing signal (default: 25000).
#' @param stepSize Integer. Step size between sliding windows (default: 5000).
#' @param txBamFile_bottom Character. Path to treatment BAM file for the bottom overlay.
#' @param inBamFile_bottom Character. Path to input/control BAM file for the bottom overlay.
#' @param txBamFile_top Character. Path to treatment BAM file for the top overlay.
#' @param inBamFile_top Character. Path to input/control BAM file for the top overlay.
#' @param yLimits Numeric vector. Y-axis limits for the plot (default: c(-0.8, 0.8)).
#' @param autoCalculateYLimits Logical. If TRUE, compute y-axis limits automatically (default: TRUE).
#' @param bottom_fill Character. Fill color for the bottom overlay (default: "#006494").
#' @param top_fill Character. Fill color for the top overlay (default: "#006494").
#' @param tx_bottom_scaleFactor Numeric or NA. Optional scaling factor for bottom treatment BAM file (default: NA).
#' @param in_bottom_scaleFactor Numeric or NA. Optional scaling factor for bottom control BAM file (default: NA).
#' @param tx_top_scaleFactor Numeric or NA. Optional scaling factor for top treatment BAM file (default: NA).
#' @param in_top_scaleFactor Numeric or NA. Optional scaling factor for top control BAM file (default: NA).
#'
#' @return A `ggplot2` object showing overlaid log2 fold change plots for two sample sets.
#'
#' @importFrom Rsamtools idxstatsBam
#' @importFrom bamsignals bamProfile
#' @importFrom zoo rollapply
#' @importFrom GenomicRanges GRanges disjoin
#' @importFrom IRanges IRanges
#' @importFrom ggplot2 ggplot geom_area scale_fill_manual scale_y_continuous scale_x_continuous theme_void theme margin element_text aes
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' makeLog2FcGenomePlotOverlay(
#'   chrom = "chr3",
#'   trackStart = 1e6,
#'   trackEnd = 2e7,
#'   bottom_fill = "steelblue",
#'   top_fill = "darkorange"
#' )
#' }
makeLog2FcGenomePlotOverlay <- function(chrom="chr3",
                                  trackStart=5000000,
                                  trackEnd=55000000,
                                  windowSize=25000,
                                  stepSize=5000,
                                  txBamFile_bottom="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                  inBamFile_bottom="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                  txBamFile_top="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                  inBamFile_top="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                  yLimits=c(-0.8,0.8),
                                  autoCalculateYLimits=TRUE,
                                  bottom_fill="#006494",
                                  top_fill="#006494",
                                  tx_bottom_scaleFactor=NA,
                                  in_bottom_scaleFactor=NA,
                                  tx_top_scaleFactor=NA,
                                  in_top_scaleFactor=NA){
  gr <- GRanges(seqnames=chrom,ranges = IRanges(start=trackStart,end=trackEnd))
  windws <- slidingWindows(gr,windowSize,step=stepSize) %>% .[[1]]
  makeCPMS <- function(BamFile=txBamFile1,scaleFactor) {
    counts <- bamProfile(BamFile, gr, paired.end = "midpoint", mapqual = 15, binsize = stepSize)
    counts_rolled <- zoo::rollapply(counts@signals[[1]], width = windowSize / stepSize, sum, align = "left")
    if (is.na(scaleFactor)) {
      totalReadCount <- idxstatsBam(BamFile) %>%
        {.[.$seqnames %in% paste0("chr", c(1:22)), ]} %>%
        {sum(as.numeric(.[, 3]))}
      cpms_rolled <- counts_rolled / totalReadCount * 1000000
    } else {
      cpms_rolled <- counts_rolled / scaleFactor
    }
    cpms_rolled_gr <- windws[1:length(cpms_rolled)] %>%
      {.$cpms <- cpms_rolled;.}
    return(cpms_rolled_gr)
  }

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

  lst_bottom <- makeDataFrameAndYLimits(txBamFile_bottom,inBamFile_bottom,tx_bottom_scaleFactor,in_bottom_scaleFactor)
  lst_top <- makeDataFrameAndYLimits(txBamFile_top,inBamFile_top,tx_top_scaleFactor,in_top_scaleFactor)

  if (autoCalculateYLimits) {
      yLimits <- c(min(c(lst_bottom[[2]][1],lst_top[[2]][1])),
        max(c(lst_bottom[[2]][2],lst_top[[2]][2])))
    }

  df_bottom <- lst_bottom[[1]] %>%
    {.$layer <- "bottom";.}
  df_top <- lst_top[[1]] %>%
    {.$layer <- "top";.}
  df <- rbind(df_bottom,df_top)

  ggplt <- ggplot(data=df,aes(x=x,y=y,fill=layer)) +
  geom_area() +
  scale_fill_manual(values = c(bottom_fill,top_fill)) +
  scale_y_continuous(limits = yLimits,breaks=c(yLimits[1],0,yLimits[2])) +
  scale_x_continuous(limits = c(trackStart,trackEnd),expand=c(0,0)) +
  theme_void() +  # Remove all axes, labels, etc.
  theme(legend.position = "none",  # Remove the legend
        plot.margin = margin(0, 0, 0, 0),
        axis.text.y  = element_text(color = "black", size = 8))

  return(ggplt)
}
