#' Generate a coverage ribbon plot from a deepTools matrix file
#'
#' This function reads a deepTools matrix file (e.g., from computeMatrix), converts it to an
#' enriched heatmap matrix, computes per-group median signal and range (min/max), and generates
#' a ribbon plot showing the average profile with variability.
#'
#' @param matrix_file Character string. Path to the input deepTools matrix file (gzipped).
#' @param group_filter_pattern Character string. Regex pattern to filter groups (e.g., "TICRR").
#' @param group_colors Named character vector of colors for the fill aesthetic, with names matching groups.
#' @param plot_title Character string. Title to display on the plot.
#' @param group_labels Optional named character vector. Names are group IDs and values are display labels.
#' @param group_title_extraction_fun Function applied to each sample name to extract the group label. 
#'   Default: \code{function(x) gsub("_rep.*", "", x)}.
#' @param x_axis_breaks Numeric vector for x-axis tick positions. Default: \code{c(0, 300, 600)}.
#' @param x_axis_labels Character vector for x-axis labels. Default: \code{c("-3kb", "summit", "+3kb")}.
#'
#' @return A \code{ggplot2} object showing a ribbon plot of group-wise median coverage with min/max shading.
#'
#' @import ggplot2
#' @import profileplyr
#' @importFrom magrittr %>%
#' @export
makeCoverageRibbonPlot <- function(matrix_file,
                                   group_filter_pattern,
                                   group_colors,
                                   plot_title = "Coverage Plot",
                                   group_labels = NULL,
                                   group_title_extraction_fun = function(x) gsub("_rep.*", "", x),
                                   x_axis_breaks = c(0, 300, 600),
                                   x_axis_labels = c("-3kb", "summit", "+3kb")) {
  # Import and convert the deepTools matrix
  proplyr_obj <- import_deepToolsMat(matrix_file)
  matrix_list <- profileplyr::convertToEnrichedHeatmapMat(proplyr_obj)

  # Compute per-sample medians
  medians_all <- lapply(names(matrix_list), function(nm) {
    mat <- matrix_list[[nm]]
    data.frame(
      median = apply(mat, 2, median),
      x = seq_len(ncol(mat)),
      sampleName = nm,
      group = group_title_extraction_fun(nm)
    )
  }) %>%
    do.call(rbind, .)

  # Compute summary statistics
  mean_df <- aggregate(median ~ x + group, data = medians_all, mean)
  min_df <- aggregate(median ~ x + group, data = medians_all, min)
  max_df <- aggregate(median ~ x + group, data = medians_all, max)

  # Assemble final plotting dataframe
  plot_df <- data.frame(
    x = mean_df$x,
    group = mean_df$group,
    y = mean_df$median,
    ymins = min_df$median,
    ymaxs = max_df$median
  )

  # Filter groups for plotting
  plot_df <- plot_df[grep(group_filter_pattern, plot_df$group), ]

  # Create plot
  p <- ggplot(plot_df, aes(x = x, y = y, ymin = ymins, ymax = ymaxs, fill = group, linetype = group)) +
    geom_line() +
    geom_ribbon(alpha = 0.5) +
    sansam_theme2() +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0, "pt"),
      legend.key.height = unit(0.1, 'cm'),
      legend.text = element_text(size = 6),
      axis.title.x = element_blank()
    ) +
    scale_x_continuous(
      breaks = x_axis_breaks,
      labels = x_axis_labels) +
    scale_fill_manual(
      name = NULL,
      labels = group_labels,
      values = group_colors
    ) +
    scale_linetype_manual(
      name = NULL,
      labels = group_labels,
      values = setNames(c(1, 2), names(group_colors))
    ) +
    ggtitle(plot_title) +
    ylab("Log2(Signal/Background)")

  return(p)
}
