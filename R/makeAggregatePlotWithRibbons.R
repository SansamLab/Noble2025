
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
