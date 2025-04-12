
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

  makeCPMS <- function(bam_file, gr_split = gr_lst) {
    counts_gr <- parallel::mclapply(gr_split, function(g) {
      tryCatch({
        g$counts <- bamsignals::bamCount(bam_file, g, paired.end = "midpoint", mapqual = mapq)
        return(g)
      }, error = function(e) {
        message(sprintf("bamCount failed on chromosome: %s for BAM file: %s", unique(as.character(seqnames(g))), bam_file))
        message("Error message: ", e$message)
        return(NULL)  # Skip this chromosome if bamCount fails
      })
    }, mc.cores = n_cores)

    # Filter out NULL results
    counts_gr <- counts_gr[!vapply(counts_gr, is.null, logical(1))]

    if (length(counts_gr) == 0) {
      stop("bamCount failed for all chromosomes in ", bam_file)
    }

    counts_gr <- GenomicRanges::GRangesList(counts_gr) %>% unlist()

    total_reads <- Rsamtools::idxstatsBam(bam_file) %>%
      subset(seqnames %in% chromosomes_to_include) %>%
      with(sum(as.numeric(mapped)))

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
