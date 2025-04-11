
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
