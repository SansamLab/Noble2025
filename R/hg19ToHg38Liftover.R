#' Lift over genomic coordinates from hg19 to hg38
#'
#' This function takes a BED file in hg19 coordinates and performs a liftOver
#' to convert the coordinates to hg38 using the provided chain file.
#'
#' @param hg19BedInput Path to the input BED file in hg19 coordinates.
#' @param hg38BedOutput Path to the output BED file in hg38 coordinates.
#'
#' @details The function downloads the hg19 to hg38 liftover chain file from
#' UCSC if it's not already present in the working directory.
#'
#' @return None (output is saved to hg38BedOutput file).
#'
#' @import rtracklayer
#'
#' @examples
#' hg19ToHg38Liftover("input_hg19.bed", "output_hg38.bed")
#'
#' @export
hg19ToHg38Liftover <- function(hg19BedInput,hg38BedOutput){
  # Define the URL of the chain file
  chain_url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

  # Define the destination file path
  dest_file <- "hg19ToHg38.over.chain.gz"

  # Download the chain file
  if (!file.exists("hg19ToHg38.over.chain.gz")) {
    download.file(chain_url, destfile = dest_file, method = "auto")
    gunzip("hg19ToHg38.over.chain.gz")
  }

  bed_file <- hg19BedInput
  bed_df <- read.table(bed_file,header=F,sep="\t",stringsAsFactors = F)
  bed_gr <- GRanges(bed_df$V1,IRanges(start = bed_df$V2, end = bed_df$V3)) %>%
    {mcols(.)<-bed_df[,-c(1,2,3)];.}

  # Load the hg19 to hg38 liftover chain filename
  hg19_to_hg38_chain <- import.chain("hg19ToHg38.over.chain")

  # Convert hg19 genomic ranges object to hg38 using the chain file
  hg38_gr <- liftOver(bed_gr, hg19_to_hg38_chain) %>% unlist

  export(hg38_gr, hg38BedOutput, format = "bed")
}
