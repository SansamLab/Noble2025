
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
