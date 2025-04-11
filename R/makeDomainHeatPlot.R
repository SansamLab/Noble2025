#' Create a Complex Heatmap of Domain Log2 Fold Change
#'
#' Generates a Complex Heatmap of log2 fold change (Lg2FC) values for genomic domains.
#'
#' @param domainsBed Path to the BED file containing genomic domains.
#' @param windowSize The size of the sliding window used for calculating Lg2FC.
#' @param stepSize The step size for sliding the window.
#' @param txBamFile1 Path to the BAM file for treatment sample 1.
#' @param inBamFile1 Path to the BAM file for control sample 1.
#' @param txBamFile2 Path to the BAM file for treatment sample 2.
#' @param inBamFile2 Path to the BAM file for control sample 2.
#' @param span The span used for generating Lg2FC values.
#'
#' @return A ComplexHeatmap object representing the log2 fold change heatmap.
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom rtracklayer readBed idxstatsBam
#' @importFrom GenomicRanges split resize width order rev GRanges
#' @importFrom Rsamtools bamProfile
#' @importFrom BiocParallel mclapply
#' @importFrom zoo rollapply
#'
#' @examples
#' heatmap <- makeDomainHeatPlot(domainsBed = "path/to/domains.bed",
#'                               windowSize = 100000,
#'                               stepSize = 5000,
#'                               txBamFile1 = "path/to/txBam1.bam",
#'                               inBamFile1 = "path/to/inBam1.bam",
#'                               txBamFile2 = "path/to/txBam2.bam",
#'                               inBamFile2 = "path/to/inBam2.bam",
#'                               span = 3000000)
#'
#' @export
makeDomainHeatPlot <- function(domainsBed="../GenerateTimingWindowsFromHCTRepliseq/TreslinDepletedLateRepDomains.bed",
                               windowSize=100000,
                               stepSize=5000,
                               txBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_250.bam",
                               inBamFile1="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_250.bam",
                               txBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/TICRR_WT_HCT116_anti-GFP_1000.bam",
                               inBamFile2="../01_Asynchronous_HCT116/results/mergedDownSampledBams/WT_HCT116_anti-GFP_1000.bam",
                               span=3000000){
  ranges_gr <- readBed(domainsBed) %>% .[order(width(.))] %>%
    .[width(.)>200000] %>%
    .[width(.)<2000000]
  #windws_lst <- slidingWindows(gr,windowSize,step=stepSize)
  makeCPMS <- function(BamFile=txBamFile1) {
    gr_lst <- split(ranges_gr,f=seqnames(ranges_gr)) %>%
      .[lapply(.,length)>0]
    counts_lst <- mclapply(gr_lst,function(gr){
      gr2 <- resize(gr,width=span,fix="center")
      df <- bamProfile(BamFile, gr2, paired.end = "midpoint", mapqual = 15, binsize = stepSize) %>%
        {.@signals} %>%
        do.call(rbind,.) %>%
        as.data.frame
      row.names(df) <- as.character(gr)
      return(df)
    })
    names(counts_lst) <- NULL
    counts_df <- do.call(rbind,counts_lst)
    counts_rolled <- apply(counts_df,1,function(cnts){
      zoo::rollapply(cnts, width = windowSize / stepSize, sum, align = "left")
    })
    sortingVector <- colnames(counts_rolled) %>% gsub(".*\\.","",.) %>%
      GRanges %>% width %>% order %>% rev
    counts_rolled_sorted <- t(counts_rolled)[sortingVector,]
    totalReadCount <- idxstatsBam(BamFile) %>%
      {
        .[.$seqnames %in% paste0("chr", c(1:22)), ]
      } %>%
      {
        sum(as.numeric(.[, 3]))
      }
    cpms_rolled <- counts_rolled_sorted / totalReadCount * 1000000

    return(cpms_rolled)
  }
  TX_CPMS1 <- makeCPMS(txBamFile1)
  IN_CPMS1 <- makeCPMS(inBamFile1)
  TX_CPMS2 <- makeCPMS(txBamFile2)
  IN_CPMS2 <- makeCPMS(inBamFile2)

  TX_CPMS <- Reduce("+",list(TX_CPMS1,TX_CPMS2))/2
  IN_CPMS <- Reduce("+",list(IN_CPMS1,IN_CPMS2))/2
  Lg2FC <- log2(TX_CPMS/(IN_CPMS))

  col_fun = colorRamp2(c(-0.5, 0, 0.5), c("green", "white", "red"))

  hmp <- ComplexHeatmap::Heatmap(Lg2FC,col=col_fun,cluster_columns = F,cluster_rows = F)

  return(hmp)
}
