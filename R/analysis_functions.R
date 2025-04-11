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
#' @param tx_scaleFactor Numeric or NA. Optional library size scaling factor for treatment BAM (default: NA).
#' @param in_scaleFactor Numeric or NA. Optional library size scaling factor for input BAM (default: NA).
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
                                  tx_scaleFactor=NA,
                                  in_scaleFactor=NA){
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


# the function below does not work well
makeAggregatePlotWithRibbons <- function(
  regions,
  useFactor=FALSE,
  mergedBamsTx=list(
    "siControl"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/ReplicatePeakAnalyzer/results/mergedDownSampledBams/controlMTBP_tx.bam",0.4140772)
      ),
    "siTreslin"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/ReplicatePeakAnalyzer/results/mergedDownSampledBams/siTicrrMTBP_tx.bam",0.3118385)
      )
    ),
  mergedBamsIn=list(
    "siControl"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/ReplicatePeakAnalyzer/results/mergedDownSampledBams/controlMTBP_in.bam",0.2851148)
      ),
    "siTreslin"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/ReplicatePeakAnalyzer/results/mergedDownSampledBams/siTicrrMTBP_in.bam",0.2851148)
      )
    ),
  replicateBamsTx=list(
    "siControl"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/MTBPctrl1.bam",0.36982481),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/MTBPctrl2.bam",0.59287014),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/MTBPctrl3.bam",0.27953668)
      ),
    "siTreslin"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/MTBPsiT1.bam",0.37347100),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/MTBPsiT2.bam",0.30103036),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/MTBPsiT3.bam",0.26101428)
      )
    ),
  replicateBamsIn=list(
    "siControl"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/HCT1.bam",0.39012439),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/HCT2.bam",0.25055516),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/HCT3.bam",0.21466478)
      ),
    "siTreslin"=list(
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/HCT1.bam",0.39012439),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/HCT2.bam",0.25055516),
      c("../2024Dec30_siTICRRonMTBP_CutRun/siTICRRonMTBP_CutRun_Analysis/results/aligned_speciesOfInterest/HCT3.bam",0.21466478)
      )
    )
  )
  {
  makeNormalizedMatrix <- function(
    bamPath = mergedBamsTx[[1]],
    Ranges = regions,
    plottingRange = 4001,
    windowSize = 150,
    stepSize = 5,
    coreCount = 8,
    chromosomes = paste0("chr",1:22)){
      countMatrix <- mclapply(paste0("chr", 1:22), function(chrm){
        chrom_ranges <- Ranges[seqnames(Ranges) == chrm]
        result <- bamProfile(bamPath[[1]], chrom_ranges, binsize = stepSize, paired.end = "midpoint", verbose = TRUE) %>%
          alignSignals(.) %>%
          as.matrix(.) %>%
          rollsumr(., k = windowSize / stepSize, align = "center") + 0.0001
        return(result)
      },
      mc.cores = coreCount) %>%
      do.call(cbind, .)

      # Helper function to calculate total read count
      getTotalReadCount <- function(bam){
        idxstatsBam(bam, index = bam) %>%
          {sapply(chromosomes, function(chr){
            .[which(.[,1] == chr), 3]
          })} %>%
          sum / 1000000
      }

      totalCount <- getTotalReadCount(bamPath[[1]])

      if (useFactor){
        normalizedCounts <- countMatrix / bamPath[[2]]
        } else {
          normalizedCounts <- countMatrix / totalCount
        }

      meanNormalizedCounts <- apply(normalizedCounts,1,mean,na.rm=TRUE)
      return(meanNormalizedCounts)
    }

  # Generate matrices of normalized counts
  mergedBamsIn_normCounts_lst <- lapply(mergedBamsIn, function(bamList) {
    lapply(bamList, makeNormalizedMatrix)
  })

  mergedBamsTx_normCounts_lst <- lapply(mergedBamsTx, function(bamList) {
    lapply(bamList, makeNormalizedMatrix)
  })

  replicateBamsIn_normCounts_lst <- lapply(replicateBamsIn, function(bamList) {
    lapply(bamList, makeNormalizedMatrix)
  })

  replicateBamsTx_normCounts_lst <- lapply(replicateBamsTx, function(bamList) {
    lapply(bamList, makeNormalizedMatrix)
  })

  # Define groups
  groups <- names(mergedBamsTx)

  makeNormalizedCountDataframe <- function(group){
    df <- data.frame("x"=NA,"y"=NA,"group"=NA,"min"=NA,"max"=NA)
    df$x <- 1:length(mergedBamsTx_normCounts_lst[[group]][[1]])
    df$y <- log2(mergedBamsTx_normCounts_lst[[group]][[1]]/(mergedBamsIn_normCounts_lst[[group]][[1]]+0.0000000001))
    replicateInMx <- replicateBamsIn_normCounts_lst[[group]] %>%
      do.call(rbind,.) %>%
      as.matrix
    replicateTxMx <- replicateBamsTx_normCounts_lst[[group]] %>%
      do.call(rbind,.) %>%
      as.matrix
    replicateMx <- log2(replicateTxMx/(replicateInMx+0.0000000001))
    df$min <- apply(replicateMx,2,min)
    df$max <- apply(replicateMx,2,max)
    df$group <- group
    return(df)
  }

  # Assemble data frames
  df <- lapply(groups,makeNormalizedCountDataframe) %>%
    do.call(rbind,.)

  ggplot(data=df, aes(x=x, y=y, color=group, ymin=min, ymax=max)) +
  geom_line() +
  theme_minimal() +
  #scale_color_manual(values = c("siControl" = siControl_color, "siTreslin" = MTBP_color)) +
  scale_x_continuous(
    breaks=c(1,386,772),
    labels = c("-2kb","Summit","+2kb"), expand=c(0,0)) +
  geom_ribbon(fill = "grey70")
}


#' Generate an UpSet plot of reproducible and replicate peaks
#'
#' This function creates an UpSet plot showing the intersection between a set of reproducible peaks
#' and multiple replicate peak sets. It pads replicate peak lists with dummy values to ensure that
#' the set size bars in the UpSet plot accurately reflect the total number of peaks in each replicate.
#' The final plot includes only reproducible peaks in the intersection matrix but uses the accurate set size bars.
#'
#' @param reproduciblePeaksPath Character string. File path to a BED file of reproducible peaks.
#' @param replicatePeakPaths Character vector. File paths to BED files containing peak calls for each replicate.
#' @param replicateLabels Character vector. Display labels for the replicate peak sets (must match order of paths).
#'
#' @return A list object generated by \code{UpSetR::upset()}, representing the UpSet plot.
#' @examples
#' makeUpSetPlot(
#'   reproduciblePeaksPath = "Reproducible_MTBP.bed",
#'   replicatePeakPaths = c("rep1.bed", "rep2.bed", "rep3.bed"),
#'   replicateLabels = c("Rep 1", "Rep 2", "Rep 3")
#' )
#' @import UpSetR
#' @export
makeUpSetPlot <- function(
    reproduciblePeaksPath = "../MakeReproduciblePeaks/MTBP_Peaks_Reproducible_6.narrowPeak",
    replicatePeakPaths = paste0("../01_Asynchronous_HCT116/results/macs2_normalPeaks/", asynchronousPeakFiles[grep("MTBP", asynchronousPeakFiles)]),
    replicateLabels = c("1:1000 1", "1:1000 2", "1:1000 3", "1:250 1", "1:250 2", "1:250 3")
) {
  # Read reproducible peaks into a GenomicRanges object
  reproduciblePeaks <- readBed(reproduciblePeaksPath)

  # Read replicate peak files into a named list of GenomicRanges
  replicatePeaksList <- lapply(replicatePeakPaths, readBed) %>%
    setNames(., replicateLabels)

  # Function to extract overlapping peaks between reproducible peaks and each replicate
  getOverlappingPeaks <- function(queryPeaks, peakList) {
    lapply(peakList, function(gr) {
      queryPeaks[queryPeaks %over% gr] %>%
        { paste0(seqnames(.), ":", start(.), "-", end(.)) }
    })
  }

  # Generate named list of reproducible peaks overlapping each replicate
  overlappingPeaks <- getOverlappingPeaks(reproduciblePeaks, replicatePeaksList) %>%
    setNames(., replicateLabels)

  # Function to pad replicate peak lists with dummy peaks so UpSetR plots total set sizes correctly
  padWithDummyPeaks <- function(sampleLabel, overlapList, allPeaksList) {
    totalPeaks <- length(allPeaksList[[sampleLabel]])
    overlapping <- length(overlapList[[sampleLabel]])
    nMissing <- totalPeaks - overlapping
    dummyPeaks <- paste0(sampleLabel, seq_len(nMissing))
    c(overlapList[[sampleLabel]], dummyPeaks)
  }

  # Pad overlapping peaks with dummy values to match full peak counts per replicate
  paddedPeaks <- lapply(names(overlappingPeaks), padWithDummyPeaks, overlappingPeaks, replicatePeaksList) %>%
    set_names(., names(overlappingPeaks))

  # Create list with only reproducible peaks for UpSetR input (for overlap-only plot)
  overlapOnlyList <- list(paste0(seqnames(reproduciblePeaks), ":", start(reproduciblePeaks), "-", end(reproduciblePeaks))) %>%
    append(overlappingPeaks, .)
  names(overlapOnlyList)[length(overlapOnlyList)] <- "Reproducible Peaks"

  # Create list with dummy-padded peak overlaps (for full peak count plot)
  paddedOverlapList <- list(paste0(seqnames(reproduciblePeaks), ":", start(reproduciblePeaks), "-", end(reproduciblePeaks))) %>%
    append(paddedPeaks, .)
  names(paddedOverlapList)[length(paddedOverlapList)] <- "Reproducible Peaks"

  # Convert both lists into UpSetR-compatible data frames
  overlapOnly_df <- UpSetR::fromList(overlapOnlyList)
  paddedOverlap_df <- UpSetR::fromList(paddedOverlapList)

  # Create UpSetR plot showing only reproducible peak overlaps
  overlapOnlyPlot <- UpSetR::upset(
    overlapOnly_df,
    sets = c(replicateLabels, "Reproducible Peaks"),
    keep.order = TRUE,
    order.by = "freq",
    number.angles = 90,
    nintersects = 20
  )

  # Create UpSetR plot including dummy peaks (for accurate set sizes)
  fullSetSizePlot <- UpSetR::upset(
    paddedOverlap_df,
    group.by = "degree",
    order.by = "freq",
    sets = c(replicateLabels, "Reproducible Peaks"),
    keep.order = TRUE,
    number.angles = 90,
    nintersects = 20
  )

  # Replace set size bars in overlap-only plot with accurate ones from the full set size plot
  overlapOnlyPlot$Sizes <- fullSetSizePlot$Sizes

  # Return the final plot
  overlapOnlyPlot
}


#' Identify reproducible peaks across multiple replicates
#'
#' This function identifies peaks from a query set that overlap with a minimum number of replicates.
#' It returns only those peaks from the query set that are present in at least \code{overlapCountThreshold}
#' of the GRanges objects in \code{subjectPeaksList}.
#'
#' @param queryPeaks A \code{GRanges} object representing the set of query peaks (e.g., reproducible peaks).
#' @param subjectPeaksList A list of \code{GRanges} objects, typically representing peak calls from multiple replicates.
#' @param overlapCountThreshold Integer. Minimum number of replicates in which a peak must be present to be considered reproducible. Default is 5.
#'
#' @return A \code{GRanges} object containing only those peaks from \code{queryPeaks} that overlap at least \code{overlapCountThreshold} replicates.
#' @examples
#' # Example usage:
#' # reproduciblePeaks <- GetReproduciblePeaks(queryPeaks, listOfReplicatePeaks, overlapCountThreshold = 4)
#'
#' @import GenomicRanges
#' @export
GetReproduciblePeaks <- function(queryPeaks,subjectPeaksList,overlapCountThreshold=5){
  lapply(subjectPeaksList,function(gr){
    queryPeaks %over% gr %>%
      as.numeric
  }) %>%
    do.call(cbind,.) %>%
    {rowSums(.)>=overlapCountThreshold} %>%
    queryPeaks[.]
}


#' Generate a data frame of overlap counts across replicate peak sets
#'
#' This function calculates how many peaks from a reproducible peak set are supported by at least
#' a given number of replicate peak sets. It returns a data frame where each row corresponds to
#' a different overlap threshold.
#'
#' @param reproduciblePeaksPath Character string. File path to a BED file containing reproducible peaks.
#' @param replicatePeakPaths Character vector. File paths to BED files of replicate peak calls.
#' @param targetName Character string. Label to associate with the resulting data (e.g., \"MTBP\", \"TRESLIN\").
#'
#' @return A data frame with columns: \code{peakCount}, \code{peakOverlapThreshold}, \code{Target}, and \code{peakPercentage}.
#'
#' @examples
#' makeOverlapThresholdSummary(
#'   reproduciblePeaksPath = "MTBP_Peaks_Reproducible_6.narrowPeak",
#'   replicatePeakPaths = c("rep1.bed", "rep2.bed", "rep3.bed"),
#'   targetName = "MTBP"
#' )
#'
#' @export
makeOverlapThresholdSummary <- function(
  reproduciblePeaksPath = "../MakeReproduciblePeaks/MTBP_Peaks_Reproducible_6.narrowPeak",
  replicatePeakPaths = paste0("../01_Asynchronous_HCT116/results/macs2_normalPeaks/", asynchronousPeakFiles[grep("MTBP", asynchronousPeakFiles)]),
  targetName = "MTBP"
) {
  # Load reproducible peaks into a GRanges object
  reproduciblePeaks <- readBed(reproduciblePeaksPath)

  # Load replicate peak sets into a list of GRanges objects
  replicatePeaksList <- lapply(replicatePeakPaths, readBed)

  # Count how many peaks from the reproducible set are found in at least i replicates (for i = 1 to N)
  overlapCountDf <- lapply(seq_along(replicatePeaksList), function(i) {
    GetReproduciblePeaks(reproduciblePeaks, replicatePeaksList, overlapCountThreshold = i) %>%
      length() %>%
      data.frame(peakCount = ., peakOverlapThreshold = i)
  }) %>%
    do.call(rbind, .) %>%
    { .$Target <- targetName; . }

  overlapCountDf$peakPercentage <- overlapCountDf$peakCount/length(reproduciblePeaks)

  return(overlapCountDf)
}


#' Generate a Correlation Matrix of log2 Fold Change in Read Counts
#'
#' This function calculates the log2 fold change (log2FC) of read counts
#' across a set of genomic regions and visualizes the sample-to-sample
#' correlation matrix using `corrplot`.
#'
#' @param peakRegionsFile Path to a BED file containing the genomic regions of interest.
#' @param treatmentBamPaths Character vector of BAM file paths for treatment samples.
#' @param inputBamPaths Character vector of BAM file paths for input/control samples.
#' @param sampleLabels Character vector of sample names to assign to columns in the output matrix.
#' @param correlationMethod Character string specifying the method used in `cor()` ("spearman", "pearson", etc.). Default is "spearman".
#'
#' @return A correlation matrix plot is displayed via `corrplot`.
#' @importFrom bamsignals bamCount
#' @importFrom corrplot corrplot
#' @export
#'
#' @examples
#' generateCorrelationMatrixOfReadLg2FcInRanges(
#'   peakRegionsFile = "path/to/peaks.bed",
#'   treatmentBamPaths = c("tx_rep1.bam", "tx_rep2.bam"),
#'   inputBamPaths = c("input_rep1.bam", "input_rep2.bam"),
#'   sampleLabels = c("rep1", "rep2")
#' )
generateCorrelationMatrixOfReadLg2FcInRanges <- function(
    peakRegionsFile = "../MakeReproduciblePeaks/MTBP_Peaks_Reproducible_6.narrowPeak",
    treatmentBamPaths,
    inputBamPaths,
    sampleLabels,
    correlationMethod = "spearman"
) {
  # Load peak regions as GRanges
  peakRegions <- readBed(peakRegionsFile)

  # Count reads in peaks for treatment and input samples
  treatmentCounts <- lapply(treatmentBamPaths, bamsignals::bamCount, gr = peakRegions, paired.end = c("midpoint"))
  inputCounts <- lapply(inputBamPaths, bamsignals::bamCount, gr = peakRegions, paired.end = c("midpoint"))

  # Combine counts into matrices
  treatmentMatrix <- do.call(cbind, treatmentCounts)
  inputMatrix <- do.call(cbind, inputCounts)

  # Compute log2 fold change: log2(Tx / Input)
  log2FCMatrix <- log2(treatmentMatrix / (inputMatrix + 1e-7))

  # Label columns for easier interpretation
  colnames(log2FCMatrix) <- sampleLabels

  # Scale values
  log2FCMatrix <- scale(log2FCMatrix)

  # Compute correlation matrix
  correlationMatrix <- cor(log2FCMatrix, method = correlationMethod)

  # Visualize correlation matrix
  list(corrplot::corrplot(correlationMatrix, method = "number"),
    correlationMatrix,
    log2FCMatrix)
}

#' Generate a Library-Normalized Read Count Matrix for Genomic Regions
#'
#' This function computes a matrix of read counts across sliding windows
#' centered on specified genomic regions, normalized by library size (CPM).
#' The output is suitable for metaplot or heatmap visualization.
#'
#' @param bamPath Character string. Path to the BAM file.
#' @param Ranges A \code{GRanges} object representing the genomic regions of interest.
#' @param plottingRange Integer. Total width of the region around each peak (default: 4001).
#' @param windowSize Integer. Size of the smoothing window (default: 150).
#' @param stepSize Integer. Bin size for signal extraction (default: 5).
#' @param coreCount Integer. Number of parallel threads to use (default: 8).
#'
#' @return A matrix of CPM-normalized read counts (rows = positions, columns = regions).
#'
#' @importFrom Rsamtools idxstatsBam
#' @importFrom bamsignals bamProfile
#' @importFrom zoo rollsumr
#' @importFrom parallel mclapply
#' @import GenomicRanges
#' @export
#'
#' @examples
#' peaks <- GRanges("chr1", IRanges(start = c(100000, 200000), width = 4001))
#' normMatrix <- makeLibraryNormalizedMatrix("example.bam", peaks)
makeLibraryNormalizedMatrix <- function(
  bamPath = "../01_Asynchronous_HCT116/results/mergedDownSampledBams/MTBP_WT_HCT116_anti-GFP_250.bam",
  Ranges,
  plottingRange = 4001,
  windowSize = 150,
  stepSize = 5,
  coreCount = 8
) {
  # Define chromosomes to analyze
  chromosomes <- paste0("chr", 1:22)

  # Resize regions to the desired plotting range
  Ranges <- resize(Ranges, width = plottingRange, fix = "center")

  # Compute raw count matrix by chromosome
  countMatrix <- mclapply(chromosomes, function(chr) {
    chrRanges <- Ranges[seqnames(Ranges) == chr]
    signalMatrix <- bamProfile(
      bamPath,
      chrRanges,
      binsize = stepSize,
      paired.end = "midpoint",
      verbose = FALSE
    ) %>%
      alignSignals() %>%
      as.matrix() %>%
      zoo::rollsumr(k = windowSize / stepSize, align = "center") + 1e-4
    return(signalMatrix)
  }, mc.cores = coreCount) %>%
    do.call(cbind, .)

  # Get total aligned reads for library size normalization
  getTotalReadCountInMillions <- function(bamFile) {
    idxstatsBam(bamFile) %>%
      subset(seqnames %in% chromosomes) %>%
      with(sum(as.numeric(mapped))) / 1e6
  }

  totalReads <- getTotalReadCountInMillions(bamPath)

  # Normalize by total reads (CPM)
  normalizedMatrix <- countMatrix / totalReads
  return(normalizedMatrix)
}

#' Import a BED file as a GenomicRanges object
#'
#' This function reads a BED file and converts it to a \code{GRanges} object,
#' retaining only the specified chromosomes.
#'
#' @param bedPath Character string. Path to the BED file.
#' @param chromosomesToImport Character vector. A list of chromosomes to keep (default is \code{chromosomes}, assumed to be defined in the environment).
#'
#' @return A \code{GRanges} object containing only the specified chromosomes.
#'
#' @examples
#' importBED("path/to/peaks.bed", chromosomesToImport = paste0("chr", 1:22))
#'
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @export
importBED <- function(bedPath, chromosomesToImport = chromosomes) {
  chromToImport <- as.character(chromosomesToImport)

  # Read first line to detect if it contains a header
  first_row <- read.table(bedPath, nrows = 1, stringsAsFactors = FALSE)

  # Detect if columns 2 and 3 (start, end) are numeric
  second_col_is_numeric <- suppressWarnings(!is.na(as.numeric(first_row[[2]])))
  third_col_is_numeric <- suppressWarnings(!is.na(as.numeric(first_row[[3]])))

  has_header <- !(second_col_is_numeric && third_col_is_numeric)

  # Read full table with or without header
  bed_df <- read.table(bedPath, header = has_header, stringsAsFactors = FALSE)

  # Standardize required column names
  colnames(bed_df)[1:3] <- c("seqnames", "start", "end")

  # Convert to GRanges and restrict to specified chromosomes
  gr <- GenomicRanges::makeGRangesFromDataFrame(bed_df, keep.extra.columns = TRUE)
  seqlevels(gr, pruning.mode = "coarse") <- chromToImport

  return(gr)
}

#' Generate a Two-Sample Euler Plot from GRanges Objects
#'
#' This function creates a reduced union of two GRanges objects, computes overlaps,
#' and plots a two-set Euler diagram using specified fill colors.
#'
#' @param gr1 GRanges. First genomic range set.
#' @param gr2 GRanges. Second genomic range set.
#' @param labels Character vector of length 2. Names for the sets (default: c("Set1", "Set2")).
#' @param fills Character vector of length 2. Fill colors for the Euler plot (default: c("grey", "white")).
#'
#' @return A ggplot2 object showing the Euler diagram.
#' @examples
#' plot <- plot_two_sample_euler(kumagai_peaks_gr, noble_peaks_gr, labels = c("Kumagai", "Noble"))
#' print(plot)
#'
#' @import GenomicRanges
#' @import IRanges
#' @import eulerr
#' @import magrittr
#'
#' @export
plot_two_sample_euler <- function(gr1, gr2, labels = c("Set1", "Set2"), fills = c("grey", "white")) {
  # Ensure GRanges inputs
  stopifnot(inherits(gr1, "GRanges"), inherits(gr2, "GRanges"))

  # Create union of input GRanges and reduce overlaps
  reduced_peaks_gr <- reduce(c(gr1, gr2))

  # Create binary overlap matrix
  attendance_test_df <- data.frame(
    reduced_peaks_gr %over% gr1,
    reduced_peaks_gr %over% gr2
  ) %>% set_names(.,labels)

  # Fit Euler diagram
  euler_fit <- eulerr::euler(attendance_test_df)

  # Plot Euler diagram
  euler_plot <- plot(
    euler_fit,
    fills = fills,
    labels = list(fontsize = 8),
    quantities = list(fontsize = 6)
  )

  return(euler_plot)
}

#' Generate a heatmap and average profile plot of log2 fold changes around summit regions
#'
#' This function computes and visualizes signal profiles (log2 fold change) around summit regions.
#' It outputs a ComplexHeatmap object and a ggplot2 line plot showing average signal per group.
#'
#' @param tx_bam Character. Path to treatment BAM file.
#' @param in_bam Character. Path to input/control BAM file.
#' @param summit_file_paths Named character vector. Paths to summit BED files for each group (e.g., TRESLIN, MTBP, Both).
#' @param peak_file_paths Named character vector. Paths to reproducible peak BED files for each group.
#' @param plotting_range Integer. Width (in bp) around each summit for signal profiling. Default: 4001.
#' @param window_size Integer. Size of the rolling window in bp. Default: 150.
#' @param step_size Integer. Step size for sliding window. Default: 5.
#' @param core_count Integer. Number of CPU cores for parallel processing. Default: 8.
#' @param chromosomes Character vector. Chromosomes to include. Default: paste0("chr", 1:22).
#' @param color_breaks Numeric vector. Values to define color scale (length should match number of colors). Default: c(-0.5, 0, 0.5, 1, 2).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{average_profile_plot}{A ggplot2 line plot showing average log2 fold change per position per group.}
#'   \item{heatmap_plot}{A ComplexHeatmap object showing per-region log2 fold change.}
#' }
#'
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom GenomicRanges resize
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar grid.grabExpr
#' @importFrom viridis viridis
#' @importFrom matrixStats rowMedians
#' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous theme_bw element_line element_blank element_text theme
#' @export
plotLog2fcHeatmap <- function(
    tx_bam,
    in_bam,
    summit_file_paths,
    peak_file_paths,
    plotting_range = 4001,
    window_size = 150,
    step_size = 5,
    core_count = 8,
    chromosomes = paste0("chr", 1:22),
    color_breaks = c(-0.5, 0, 0.5, 1, 2)
) {
  # Import peaks and assign labels
  peaks <- lapply(peak_file_paths, importBED,chromosomes)
  all_peaks <- as(peaks, "GRangesList") %>% unlist()
  peak_sets <- lapply(names(peaks), function(nme) {
    gr <- as(peaks, "GRangesList")[[nme]]
    all_peaks %over% gr %>%
      as.character() %>%
      gsub("TRUE", nme, .) %>%
      gsub("FALSE", "", .)
  }) %>%
    do.call(paste, .) %>%
    gsub("^ *", "", .) %>%
    gsub(" *$", "", .)

  # Import summit GRanges and resize
  summit_granges_list <- lapply(summit_file_paths, importBED,chromosomes)
  all_summits_gr <- as(summit_granges_list, "GRangesList") %>% unlist()
  plotting_regions_gr <- GenomicRanges::resize(all_summits_gr, width = plotting_range, fix = "center")

  # Compute matrices
  tx_matrix <- makeLibraryNormalizedMatrix(
    bamPath = tx_bam,
    Ranges = plotting_regions_gr,
    plottingRange = plotting_range,
    windowSize = window_size,
    stepSize = step_size,
    coreCount = core_count
  )
  in_matrix <- makeLibraryNormalizedMatrix(
    bamPath = in_bam,
    Ranges = plotting_regions_gr,
    plottingRange = plotting_range,
    windowSize = window_size,
    stepSize = step_size,
    coreCount = core_count
  )

  log2fc_matrix <- log2(tx_matrix / in_matrix) %>% t()

  trim_length <- round((ncol(log2fc_matrix) - ncol(log2fc_matrix) * 0.1) / 2)
  log2fc_matrix_sorted <- log2fc_matrix[rev(order(rowMedians(log2fc_matrix[, trim_length:(ncol(log2fc_matrix) - trim_length)]))), ]
  colnames(log2fc_matrix_sorted) <- rep("", ncol(log2fc_matrix_sorted)) %>%
    { .[length(.) / 2] <- "Summit"; . } %>%
    { .[1] <- paste0("-", round(plotting_range / 2000, 1), " kb"); . } %>%
    { .[length(.)] <- paste0("+", round(plotting_range / 2000, 1), " kb"); . }
  peak_set_labels_sorted <- peak_sets[rev(order(rowMedians(log2fc_matrix[, trim_length:(ncol(log2fc_matrix) - trim_length)])))]

  col_fun <- circlize::colorRamp2(
    color_breaks,
    viridis::viridis(length(color_breaks))
  )

  heatmap_plot <- ComplexHeatmap::Heatmap(
    log2fc_matrix_sorted,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    col = col_fun,
    row_split = peak_set_labels_sorted,
    row_title_gp = grid::gpar(fontsize = 8),
    column_title_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 8),
      labels_gp = grid::gpar(fontsize = 7)
    ),
    column_names_rot = 45,
    column_names_centered = FALSE,
    show_heatmap_legend = FALSE
  ) %>% ComplexHeatmap::draw() %>% grid::grid.grabExpr()

  matrix_list_by_set <- split(as.data.frame(log2fc_matrix_sorted), f = peak_set_labels_sorted)
  average_log2fc_df <- lapply(names(matrix_list_by_set), function(set_name) {
    df <- matrix_list_by_set[[set_name]]
    data.frame(
      log2FC = apply(t(df), 1, mean),
      index = seq_len(ncol(df)),
      set = set_name
    )
  }) %>%
    do.call(rbind, .)

  average_profile_plot <- ggplot2::ggplot(average_log2fc_df, aes(x = index, y = log2FC, group = set, color = set)) +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      axis.line = ggplot2::element_line(color = "black"),
      plot.background = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )

  return(list(
    average_profile_plot = average_profile_plot,
    heatmap_plot = heatmap_plot
  ))
}

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


#' Create Genome Track Plot for TRESLIN, MTBP, and Replication Features
#'
#' Generates a genome browser-style figure for a specified genomic region,
#' visualizing log2 fold change signals from TRESLIN and MTBP ChIP-seq,
#' peak annotations, and origin mapping (Repli-seq, SNS-seq, Ini-seq).
#' Allows for customizable track inclusion and sizing via a named vector input.
#'
#' @param range_to_plot UCSC-style genomic region to visualize (e.g. "chr1:40000000-120000000").
#' @param treslin_tx_bam_1_path Path to TRESLIN ChIP replicate 1 BAM file.
#' @param treslin_in_bam_1_path Path to input control for TRESLIN replicate 1.
#' @param treslin_tx_bam_2_path Path to TRESLIN ChIP replicate 2 BAM file.
#' @param treslin_in_bam_2_path Path to input control for TRESLIN replicate 2.
#' @param treslin_color Color for TRESLIN signal and peaks. Default: "#006494"
#' @param treslin_log2fc_ylimits Y-axis limits for TRESLIN log2FC signal. Default: c(-0.8, 0.8)
#' @param treslin_log2fc_step_size Step size for sliding window of TRESLIN signal. Default: 5000
#' @param treslin_log2fc_window_size Window size for summing reads (TRESLIN). Default: 25000
#' @param mtbp_tx_bam_1_path Path to MTBP ChIP replicate 1 BAM file.
#' @param mtbp_in_bam_1_path Path to input control for MTBP replicate 1.
#' @param mtbp_tx_bam_2_path Path to MTBP ChIP replicate 2 BAM file.
#' @param mtbp_in_bam_2_path Path to input control for MTBP replicate 2.
#' @param mtbp_color Color for MTBP signal and peaks. Default: "#55C1FF"
#' @param mtbp_log2fc_ylimits Y-axis limits for MTBP log2FC signal. Default: c(-0.8, 0.8)
#' @param mtbp_log2fc_step_size Step size for sliding window of MTBP signal. Default: 5000
#' @param mtbp_log2fc_window_size Window size for summing reads (MTBP). Default: 25000
#' @param treslin_peaks_bed_path Path to BED file with TRESLIN peaks.
#' @param mtbp_peaks_bed_path Path to BED file with MTBP peaks.
#' @param snsseq_bed_path Path to BED file for SNS-seq origins.
#' @param snsseq_track_color Color for SNS-seq annotation track. Default: "#C1666B"
#' @param iniseq_bed_path Path to BED file for Ini-seq origins.
#' @param iniseq_track_color Color for Ini-seq annotation track. Default: "#D4B483"
#' @param repliseq_gr_rds_path Path to RDS file containing Repli-seq GRanges.
#' @param ideogram_csv_path Path to CSV file with ideogram data (columns: chrom, chromStart, chromEnd, name, gieStain).
#' @param highlight_regions_starts Numeric vector of start coordinates for highlight regions.
#' @param highlight_regions_widths Numeric vector of widths for highlight regions.
#' @param figure_width Width of output figure in inches. Default: 16
#' @param figure_height Height of output figure in inches. Default: 16
#' @param genome_build Genome build identifier (e.g. "hg38"). Default: "hg38"
#' @param tracks_to_plot_with_sizes Named numeric vector. Names must match track objects created inside the function. Values control their heights. Tracks appear in the specified order.
#'
#' @return A grid graphics object representing the genome track plot.
#'
#' @import Gviz
#' @import grid
#' @import magrittr
#' @import utils
#' @export

createAsynchronousTracksWithRepliseq <- function(
    range_to_plot = "chr1:40000000-120000000",
    treslin_tx_bam_1_path,
    treslin_in_bam_1_path,
    treslin_tx_bam_2_path,
    treslin_in_bam_2_path,
    treslin_color = "#006494",
    treslin_log2fc_ylimits = c(-0.8, 0.8),
    treslin_log2fc_step_size = 5000,
    treslin_log2fc_window_size = 25000,
    mtbp_tx_bam_1_path,
    mtbp_in_bam_1_path,
    mtbp_tx_bam_2_path,
    mtbp_in_bam_2_path,
    mtbp_color = "#55C1FF",
    mtbp_log2fc_ylimits = c(-0.8, 0.8),
    mtbp_log2fc_step_size = 5000,
    mtbp_log2fc_window_size = 25000,
    treslin_peaks_bed_path,
    mtbp_peaks_bed_path,
    snsseq_bed_path,
    snsseq_track_color = "#C1666B",
    iniseq_bed_path,
    iniseq_track_color = "#D4B483",
    repliseq_gr_rds_path,
    ideogram_csv_path,
    highlight_regions_starts,
    highlight_regions_widths,
    figure_width = 16,
    figure_height = 16,
    genome_build = "hg38",
    tracks_to_plot_with_sizes = c(
      "itrack" = 0.4,
      "axisTrack" = 0.8,
      "MTBP_Lg2FC_dTrack" = 1,
      "MTBP_Unified_Peaks_Track" = 0.5,
      "TRESLIN_Lg2FC_dTrack" = 1,
      "TRESLIN_Unified_Peaks_Track" = 0.5,
      "IniSeqTrack" = 0.5,
      "SnsSeq_Track" = 0.5,
      "dTrack_Repliseq" = 1
    )
){
  setCustomGvizScheme()
  track_sizes <- as.numeric(tracks_to_plot_with_sizes)

  message("Parsing range: ", range_to_plot)
  range_to_plot_split <- strsplit(range_to_plot, ":|-") %>%
    unlist() %>%
    magrittr::set_names(c("chromosome", "start", "end"))

  chromosome <- range_to_plot_split[["chromosome"]]
  start_pos <- as.numeric(range_to_plot_split[["start"]])
  end_pos <- as.numeric(range_to_plot_split[["end"]])

  message("Chromosome: ", chromosome)
  message("Start: ", start_pos)
  message("End: ", end_pos)

  message("Loading ideogram from: ", ideogram_csv_path)
  ideogram_df <- read.csv(ideogram_csv_path)
  names(ideogram_df) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
  itrack <- IdeogramTrack(genome = genome_build, bands = ideogram_df, ucscChromosomeNames = FALSE, showId = FALSE)

  axisTrack <- GenomeAxisTrack(fontsize = 8, distFromAxis = 9)

  message("Loading Repli-seq from: ", repliseq_gr_rds_path)
  dTrack_Repliseq <- readRDS(repliseq_gr_rds_path) %>%
    DataTrack(
      name = "RepliSeq", type = "heatmap",
      gradient = c("white", "#A1A0A1", "#4B4A4B", "black"),
      ylim = c(5, 20), showTitle = FALSE, genome = genome_build
    )

  message("Importing SNS-seq BED from: ", snsseq_bed_path)
  SnsSeq_Track <- AnnotationTrack(
    range = importBED(snsseq_bed_path, chromosomesToImport = chromosome),
    name = "SNSSeq", col = snsseq_track_color,
    fill = snsseq_track_color, col.line = snsseq_track_color,
    showTitle = FALSE, genome = genome_build
  )

  message("Importing Ini-seq BED from: ", iniseq_bed_path)
  IniSeqTrack <- AnnotationTrack(
    range = importBED(iniseq_bed_path, chromosomesToImport = chromosome),
    name = "IniSeq", col = iniseq_track_color,
    fill = iniseq_track_color, col.line = iniseq_track_color,
    showTitle = FALSE, genome = genome_build
  )

  message("Importing MTBP peaks from: ", mtbp_peaks_bed_path)
  MTBP_Unified_Peaks_Track <- AnnotationTrack(
    range = importBED(mtbp_peaks_bed_path, chromosomesToImport = chromosome),
    name = "MTBP peaks", col = mtbp_color,
    fill = mtbp_color, col.line = mtbp_color,
    showTitle = FALSE, genome = genome_build
  )

  message("Importing TRESLIN peaks from: ", treslin_peaks_bed_path)
  TRESLIN_Unified_Peaks_Track <- AnnotationTrack(
    range = importBED(treslin_peaks_bed_path, chromosomesToImport = chromosome),
    name = "TRESLIN peaks", col = treslin_color,
    fill = treslin_color, col.line = treslin_color,
    showTitle = FALSE, genome = genome_build
  )

  message("Creating TRESLIN signal track...")
  TRESLIN_Lg2FC_dTrack <- make2SampleLog2FcCoverageDataTrack(
    chromosome = chromosome,
    start = start_pos,
    end = end_pos,
    windowSize = treslin_log2fc_window_size,
    stepSize = treslin_log2fc_step_size,
    txBamFile1 = treslin_tx_bam_1_path,
    inBamFile1 = treslin_in_bam_1_path,
    txBamFile2 = treslin_tx_bam_2_path,
    inBamFile2 = treslin_in_bam_2_path,
    trackName = "TRESLIN log2FC",
    HistogramColor = treslin_color,
    yLimits = treslin_log2fc_ylimits
  )

  message("Creating MTBP signal track...")
  MTBP_Lg2FC_dTrack <- make2SampleLog2FcCoverageDataTrack(
    chromosome = chromosome,
    start = start_pos,
    end = end_pos,
    windowSize = mtbp_log2fc_window_size,
    stepSize = mtbp_log2fc_step_size,
    txBamFile1 = mtbp_tx_bam_1_path,
    inBamFile1 = mtbp_in_bam_1_path,
    txBamFile2 = mtbp_tx_bam_2_path,
    inBamFile2 = mtbp_in_bam_2_path,
    trackName = "MTBP log2FC",
    HistogramColor = mtbp_color,
    yLimits = mtbp_log2fc_ylimits
  )

  all_tracks_list <- list(
    "itrack" = itrack, "axisTrack" = axisTrack,
    "MTBP_Lg2FC_dTrack" = MTBP_Lg2FC_dTrack,
    "MTBP_Unified_Peaks_Track" = MTBP_Unified_Peaks_Track,
    "TRESLIN_Lg2FC_dTrack" = TRESLIN_Lg2FC_dTrack,
    "TRESLIN_Unified_Peaks_Track" = TRESLIN_Unified_Peaks_Track,
    "IniSeqTrack" = IniSeqTrack, "SnsSeq_Track" = SnsSeq_Track,
    "dTrack_Repliseq" = dTrack_Repliseq
  )

  tracks_to_plot <- all_tracks_list[names(tracks_to_plot_with_sizes)]

  message("Creating highlight track and rendering figure...")
  ht <- HighlightTrack(
    trackList = tracks_to_plot,
    start = highlight_regions_starts,
    width = highlight_regions_widths,
    chromosome = chromosome,
    fill = rep("transparent", length(highlight_regions_starts)),
    col = c("#ED6A5A"),
    inBackground = FALSE, lwd = 2
  )

  track_plot <- plotTracks(
    ht,
    chromosome = chromosome,
    from = start_pos,
    to = end_pos,
    sizes = track_sizes,
    showTitle = FALSE
  ) %>% grid.grabExpr(width = figure_width, height = figure_height)

  message("Track plot complete.")
  return(track_plot)
}


#' Apply Custom Gviz Scheme for Consistent Track Styling
#'
#' This function defines and registers a customized Gviz visual scheme called `"myScheme"`,
#' which sets consistent font sizes, colors, and graphical scaling for genome tracks.
#' It then sets it as the default Gviz scheme via `options(Gviz.scheme = "myScheme")`.
#'
#' @return Invisibly returns the name of the applied scheme.
#' @import Gviz
#' @export
setCustomGvizScheme <- function() {
  scheme <- Gviz::getScheme()
  cex <- 1

  scheme$GdObject$col <- "black"
  scheme$GdObject$fontsize <- 8
  scheme$DataTrack$fontsize.legend <- 8
  scheme$IdeogramTrack$fontsize <- 8
  scheme$AnnotationTrack$fontsize.group <- 8

  scheme$GenomeAxisTrack$cex.id <- cex
  scheme$GenomeAxisTrack$cex <- cex
  scheme$GdObject$cex.axis <- cex
  scheme$GdObject$cex.title <- cex
  scheme$DataTrack$cex.legend <- cex
  scheme$DataTrack$cex.sampleNames <- cex
  scheme$DataTrack$cex <- cex
  scheme$DataTrack$cex.axis <- cex
  scheme$DataTrack$cex.title <- cex
  scheme$IdeogramTrack$cex.bands <- cex
  scheme$IdeogramTrack$cex <- cex
  scheme$AnnotationTrack$cex <- cex
  scheme$AnnotationTrack$cex.group <- cex

  scheme$GenomeAxisTrack$fontcolor <- "black"
  scheme$DataTrack$fontcolor.legend <- "black"
  scheme$GdObject$fontcolor <- "black"
  scheme$IdeogramTrack$fontcolor <- "black"
  scheme$AnnotationTrack$fontcolor.group <- "black"
  scheme$AnnotationTrack$fontcolor.item <- "black"

  scheme$GdObject$col.axis <- "black"
  scheme$GdObject$col.border.title <- "transparent"
  scheme$GdObject$col.frame <- "transparent"
  scheme$GdObject$col.title <- "black"
  scheme$GdObject$col.id <- "black"
  scheme$GenomeAxisTrack$col <- "black"
  scheme$DataTrack$col.sampleNames <- "black"
  scheme$GdObject$background.title <- "white"
  scheme$GdObject$fontcolor.title <- "black"

  scheme$AnnotationTrack$rotation <- 90
  scheme$AnnotationTrack$rotation.item <- 90
  scheme$AnnotationTrack$rotation.group <- 90
  scheme$GdObject$rotation.title <- 0
  scheme$GdObject$rotation <- 90

  Gviz::addScheme(scheme, "myScheme")
  options(Gviz.scheme = "myScheme")

  invisible("myScheme")
}

#' Calculate log2 fold change of BAM signals over GRanges windows
#'
#' This function calculates the log2 fold change (log2FC) of coverage per million (CPM)
#' between treatment and input BAM files over a set of GRanges, typically corresponding
#' to windows such as those used in Repli-seq.
#'
#' @param tx_bam_files Character vector of paths to treatment BAM files.
#' @param in_bam_files Character vector of paths to input BAM files.
#' @param gr A \code{GRanges} object representing regions over which to calculate CPM and log2FC.
#' @param chromosomes_to_include Character vector of chromosome names to include (default: paste0("chr", 1:22)).
#' @param mapq Integer. Minimum mapping quality to include (default: 15).
#' @param n_cores Number of cores for parallel processing (default: 8).
#'
#' @return A GRanges object with a \code{Log2FC} metadata column.
#'
#' @import GenomicRanges
#' @import Rsamtools
#' @import bamsignals
#' @import parallel
#' @import magrittr
#' @export
calculateLog2FcOverGRanges <- function(
    tx_bam_files,
    in_bam_files,
    gr,
    chromosomes_to_include = paste0("chr", 1:22),
    mapq = 15,
    n_cores = 8
) {
  stopifnot(length(tx_bam_files) >= 1, length(in_bam_files) >= 1)

  # Split GRanges by chromosome to reduce memory usage
  gr_lst <- split(gr, seqnames(gr))

  # Count and normalize to CPM
  makeCPMS <- function(bam_file, gr_split = gr_lst) {
    counts_gr <- parallel::mclapply(gr_split, function(g) {
      g$counts <- bamsignals::bamCount(bam_file, g, paired.end = "midpoint", mapqual = mapq)
      return(g)
    }, mc.cores = n_cores) %>% GenomicRanges::GRangesList() %>% unlist()

    total_reads <- Rsamtools::idxstatsBam(bam_file) %>%
      { .[.$seqnames %in% chromosomes_to_include, ] } %>%
      { sum(as.numeric(.[, 3])) }

    counts_gr$cpms <- counts_gr$counts / total_reads * 1e6
    return(counts_gr)
  }

  # Average CPM across replicates
  makeAllCpms <- function(bam_list) {
    if (length(bam_list) == 1) {
      return(makeCPMS(bam_list[[1]]))
    } else {
      cpm_gr_list <- lapply(bam_list, makeCPMS)
      mean_cpms <- lapply(cpm_gr_list, function(gr) matrix(gr$cpms, ncol = 1)) %>%
        do.call(cbind, .) %>%
        rowMeans()

      cpm_gr_list[[1]]$cpms <- mean_cpms
      return(cpm_gr_list[[1]])
    }
  }


  tx_cpms_gr <- makeAllCpms(tx_bam_files)
  in_cpms_gr <- makeAllCpms(in_bam_files)

  tx_cpms_gr$Log2FC <- log2(tx_cpms_gr$cpms / in_cpms_gr$cpms)
  return(tx_cpms_gr)
}

#' Annotate Log2FC GRanges by Replication Timing Fraction
#'
#' Matches log2 fold change values from a GRanges object to high-resolution replication timing (Repli-seq)
#' fractions, assigning each region to the S-phase fraction where replication signal is maximal.
#'
#' @param log2fc_gr A `GRanges` object with a metadata column named `Log2FC`.
#' @param repliseq_gr A `GRanges` object with Repli-seq data across S-phase fractions. Each metadata column
#'   corresponds to an S-phase fraction (e.g., S1 to S16).
#'
#' @return A `data.frame` with columns:
#'   \describe{
#'     \item{S.fraction}{S-phase fraction with max Repli-seq signal (e.g., "S1", "S2", ...).}
#'     \item{Log2FC}{Log2 fold change value from the input `log2fc_gr`.}
#'   }
#'
#' @import GenomicRanges
#' @import S4Vectors
#' @import GenomeInfoDb
#' @export
summarize_log2fc_by_rt_fraction <- function(log2fc_gr, repliseq_gr) {
  repliseq_clean <- repliseq_gr[complete.cases(mcols(repliseq_gr))] %>%
    sortSeqlevels() %>% sort()

  repliseq_clean$Max <- apply(mcols(repliseq_clean), 1, which.max) %>% paste0("S", .)

  log2fc_gr_filtered <- log2fc_gr[!is.na(match(log2fc_gr, repliseq_clean))] %>%
    sortSeqlevels() %>% sort()

  data.frame(
    "S.fraction" = factor(repliseq_clean$Max, levels = paste0("S", 1:16)),
    "Log2FC" = log2fc_gr_filtered$Log2FC
  )
}

#' Plot Median Log2 Fold Change by Replication Timing Fraction
#'
#' Creates a line plot of the median log2 fold change across S-phase fractions (S1 to S16).
#' This is useful for visualizing how signal enrichment or depletion changes across the replication timing program.
#'
#' @param df A `data.frame` with columns `S.fraction` (factor) and `Log2FC` (numeric),
#'   typically generated by `summarize_log2fc_by_rt_fraction()`.
#' @param title Optional plot title to display.
#'
#' @return A `ggplot` object.
#'
#' @import ggplot2
#' @import magrittr
#' @export
plot_log2fc_median_by_fraction <- function(df, title = NULL) {
  df_median <- aggregate(Log2FC ~ S.fraction, df, FUN = function(x) median(x[is.finite(x)]))
  df_median$S.fraction.numeric <- gsub("S", "", df_median$S.fraction) %>% as.numeric()

  ggplot(df_median, aes(x = S.fraction.numeric, y = Log2FC)) +
    geom_line() +
    geom_point(size = 1) +
    theme_bw() +
    sansam_theme() +
    xlab("Rep. timing") +
    ylab("Log2(Sig/Bkg)") +
    ggtitle(title)
}

#' Calculate Peak Density by Replication Timing Fraction
#'
#' Computes the density of peaks (overlaps) in each replication timing fraction
#' based on overlap with a set of peaks (e.g., ChIP or origin peaks). Peak density
#' is normalized per 50 kb.
#'
#' @param repliseq_gr A \code{GRanges} object containing replication timing windows with 16 S phase fractions.
#' @param peak_gr A \code{GRanges} object containing genomic regions such as peaks.
#'
#' @return A \code{data.frame} with one row per replication timing fraction and the following columns:
#' \itemize{
#'   \item \code{S.fraction}: Fraction label (S1–S16)
#'   \item \code{bp}: Total base pairs in that S fraction
#'   \item \code{peakCount}: Total number of peaks overlapping the S fraction
#'   \item \code{meanPeakCount}: Mean number of peaks per 50kb window in the S fraction
#'   \item \code{peakDensity}: Normalized peak density (peaks per 50 kb)
#'   \item \code{S.fraction.numeric}: Numeric version of the S fraction (1–16)
#' }
#'
#' @import GenomicRanges
#' @import IRanges
#' @export
calculate_peak_density_by_fraction <- function(repliseq_gr, peak_gr) {
  repliseq_clean <- repliseq_gr[complete.cases(mcols(repliseq_gr))] %>%
    sortSeqlevels() %>%
    sort()

  repliseq_clean$OverlapCounts <- countOverlaps(repliseq_clean, peak_gr)
  repliseq_clean$S.fraction <- apply(mcols(repliseq_clean), 1, which.max) %>% paste0("S", .)
  repliseq_clean$bp <- width(repliseq_clean)

  df <- as.data.frame(mcols(repliseq_clean))
  df2 <- aggregate(bp ~ S.fraction, df, sum)
  df2$peakCount <- aggregate(OverlapCounts ~ S.fraction, df, sum)$OverlapCounts
  df2$meanPeakCount <- aggregate(OverlapCounts ~ S.fraction, df, mean)$OverlapCounts
  df2$peakDensity <- df2$peakCount / df2$bp * 50000
  df2$S.fraction.numeric <- gsub("S", "", df2$S.fraction) %>% as.numeric()

  return(df2)
}

#' Plot Peak Density by Replication Timing Fraction
#'
#' Generates a line plot showing peak density across S phase replication timing fractions.
#'
#' @param df A data frame returned from \code{calculate_peak_density_by_fraction()}.
#' @param title Optional title for the plot. Default is \code{NULL}.
#'
#' @return A \code{ggplot} object.
#'
#' @import ggplot2
#' @export
#'
plot_peak_density_by_fraction <- function(df, title = NULL) {
  ggplot(df, aes(x = S.fraction.numeric, y = peakDensity)) +
    geom_line() +
    geom_point(size = 1) +
    theme_bw() +
    sansam_theme() +
    xlab("Rep. timing") +
    ylab("Peaks per 50 kb") +
    ggtitle(title)
}

#' Annotate and label a set of genomic peaks
#'
#' Uses ChIPseeker::annotatePeak to annotate peaks and adds a column
#' indicating the peak set.
#'
#' @param peaks A GRanges object representing peaks to annotate.
#' @param set_name Character string labeling the source of the peaks.
#' @param txdb A TxDb object (e.g. from GenomicFeatures::makeTxDbFromGFF).
#' @param tss_region A numeric vector of length 2 specifying the TSS region.
#' Default: c(-3000, 3000)
#'
#' @return A data.frame of annotation statistics with a new column "Set"
#' @import ChIPseeker
#' @import dplyr
annotate_and_label_peaks <- function(peaks, set_name, txdb, tss_region = c(-3000, 3000)) {
  ChIPseeker::annotatePeak(peaks, tssRegion = tss_region, TxDb = txdb, annoDb = "org.Hs.eg.db")@annoStat %>%
    dplyr::mutate(Set = set_name)
}


#' Combine multiple labeled peak annotation sets into a single data.frame
#'
#' Calls `annotate_and_label_peaks()` for each set and binds the results
#' into a single tidy data frame.
#'
#' @param peak_sets Named list of GRanges objects representing different peak sets.
#' @param txdb A TxDb object for gene annotation.
#' @param tss_region A numeric vector of length 2 for TSS region. Default: c(-3000, 3000)
#'
#' @return A data.frame containing annotation summaries for all peak sets.
#' @import dplyr
#' @export
combine_peak_annotation_sets <- function(peak_sets, txdb, tss_region = c(-3000, 3000)) {
  annotation_dfs <- lapply(names(peak_sets), function(set_name) {
    annotate_and_label_peaks(
      peaks = peak_sets[[set_name]],
      set_name = set_name,
      txdb = txdb,
      tss_region = tss_region
    )
  })

  dplyr::bind_rows(annotation_dfs)
}

#' Plot Proportion of Genomic Features by Peak Set
#'
#' @param annotation_df Combined annotation data frame.
#' @param fill_colors Named vector of fill colors for genomic features.
#'
#' @return A ggplot object showing stacked bar chart of annotation percentages.
#' @import ggplot2
#' @export
plot_peak_annotation_proportions <- function(annotation_df, fill_colors) {
  ggplot(annotation_df, aes(fill = Feature, y = Frequency, x = Set)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = fill_colors) +
    theme_bw(base_size = 8) +
    theme(
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 8),
      legend.key.size = unit(2, 'mm'),
      panel.background = element_rect(fill = 'transparent'),
      plot.background = element_rect(fill = 'transparent', color = NA)
    ) +
    coord_flip()
}
