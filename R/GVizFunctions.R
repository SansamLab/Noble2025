#' Create a DataTrack for Log2 Fold Change (Log2FC) Coverage Data
#'
#' This function generates a DataTrack for Log2FC coverage data based on
#' input BAM files. It calculates Log2FC values for specified genomic regions
#' and creates a DataTrack for visualization.
#'
#' @param chromosome The chromosome for the genomic region of interest.
#' @param start The start position of the genomic region.
#' @param end The end position of the genomic region.
#' @param windowSize The size of the sliding windows for Log2FC calculation.
#' @param stepSize The step size for sliding windows.
#' @param txBamFile The path to the treatment BAM file.
#' @param inBamFile The path to the control BAM file.
#' @param trackName The name for the DataTrack.
#' @param HistogramColor The color for the histogram in the DataTrack.
#'
#' @return A DataTrack object for Log2FC coverage data.
#'
#' @import GenomicRanges
#' @import zoo
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' dt <- makeLog2FcCoverageDataTrack(chromosome="chr3", start=5000000, end=55000000,
#'                                   windowSize=25000, stepSize=5000,
#'                                   txBamFile="path/to/txBamFile.bam",
#'                                   inBamFile="path/to/inBamFile.bam",
#'                                   trackName="Log2FC", HistogramColor="#1b9e77")
#' }
#'

makeLog2FcCoverageDataTrack <- function(chromosome="chr3",
                                        start=5000000,
                                        end=55000000,
                                        windowSize=25000,
                                        stepSize=5000,
                                        txBamFile="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                        inBamFile="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                        trackName="Log2FC",
                                        HistogramColor="#1b9e77"
){
  gr <- GRanges(seqnames=chromosome,ranges = IRanges(start=start,end=end))
  windws <- slidingWindows(gr,windowSize,step=stepSize) %>% .[[1]]
  makeCPMS <- function(BamFile=txBamFile) {
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
  TX_CPMS <- makeCPMS()
  IN_CPMS <- makeCPMS(inBamFile)
  output_gr <- TX_CPMS
  mcols(output_gr) <- NULL
  mcols(output_gr) <- log2(TX_CPMS$cpms/(IN_CPMS$cpms+0.0000000001))
  dt <- DataTrack(
    range = output_gr, genome = "hg38", type = "hist",name = trackName,
    #col = HistogramColor, col.histogram = HistogramColor, fill.histogram = HistogramColor,
    fill.histogram = HistogramColor,col.histogram = "transparent",col="transparent",
    lwd = 0.5, ylim = c(quantile(output_gr$X,0.001), quantile(output_gr$X,0.999))
  )
  return(dt)
}




make2SampleCpmCoverageDataTrack <- function(chrom="chr3",
                                               trackStart=5000000,
                                               trackEnd=55000000,
                                               windowSize=25000,
                                               stepSize=5000,
                                               txBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                                               inBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                                               txBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_1000.bam",
                                               inBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_1000.bam",
                                               trackName="cpm",
                                               HistogramColor="#1b9e77"
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

  makeOverlayTrack <- function(tx_gr=tx_cpms_gr,in_gr=in_cpms_gr){
    yLimits <- c(
      quantile(
        c(tx_gr$cpms, in_gr$cpms),
        0.0001
      ),
      quantile(
        c(tx_gr$cpms, in_gr$cpms),
        0.9999
      )
    ) %>% {.*1.2}
    dt_tx <- DataTrack(
      range=tx_gr,
      genome="hg38",
      type="hist",
      name=trackName,
      fill.histogram=HistogramColor,
      col.histogram = "transparent",
      col="transparent",
      lwd=0.5, ylim=yLimits)
    dt_in <- DataTrack(
      range=in_gr,
      genome="hg38",
      type="hist",
      name=trackName,
      fill.histogram="grey",
      col.histogram = "transparent",
      col="transparent",
      lwd=0.5, ylim=yLimits)
    ot <- OverlayTrack(trackList=list(dt_tx, dt_in))
    return(ot)
  }


  ot <- makeOverlayTrack()

  return(ot)
}

# make2SampleCpmGRanges

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
