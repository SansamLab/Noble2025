#' Create a Stacked Bar Chart of Genomic Features Overlapping Peaks
#'
#' This function generates a stacked bar chart representing the frequency of genomic
#' features overlapping with peaks.
#'
#' @param df_of_genomic_features_overlapping_peaks A data frame containing the frequency of genomic features overlapping with peaks.
#'   Expected columns:
#'   - `Feature`: Genomic feature (e.g., "Promoter", "Exon", "Intron").
#'   - `Frequency`: The number of peaks overlapping each feature.
#'   - `Set`: The label for the peak set.
#' @param feature_colors A named vector of colors for each genomic feature.
#'   The default color scheme follows a standardized palette for promoters, UTRs, exons, introns, and intergenic regions.
#'
#' @return A `ggplot` object representing a stacked bar chart of genomic feature frequencies.
#'
#' @import ggplot2
#' @importFrom grid unit
#' @export
#' @seealso \code{\link{make_df_of_genomic_features_overlapping_peaks}}
#'
#' @examples
#' df <- data.frame(
#'   Feature = c("Promoter (<=1kb)", "Exon", "Intron", "Distal Intergenic"),
#'   Frequency = c(100, 200, 300, 150),
#'   Set = c("Peaks Set 1", "Peaks Set 1", "Peaks Set 1", "Peaks Set 1")
#' )
#' make_stacked_barchart_of_genomic_features_overlapping_peaks(df)
make_stacked_barchart_of_genomic_features_overlapping_peaks <- function(
  df_of_genomic_features_overlapping_peaks,
  feature_colors = c(
    "Promoter (<=1kb)" = "#9b2226",
    "Promoter (1-2kb)" = "#ae2012",
    "Promoter (2-3kb)" = "#bb3e03",
    "5' UTR" = "#ca6702",
    "3' UTR" = "#ee9b00",
    "1st Exon" = "#e9d8a6",
    "Other Exon" = "#F7F1DE",
    "1st Intron" = "#94d2bd",
    "Other Intron" = "#0a9396",
    "Downstream (<=300)" = "#005f73",
    "Distal Intergenic" = "#001219"
  )
) {
  # Ensure the data frame has the required columns
  if (!all(c("Feature", "Frequency", "Set") %in% colnames(df_of_genomic_features_overlapping_peaks))) {
    stop("Error: Input data frame must contain 'Feature', 'Frequency', and 'Set' columns.")
  }

  # Generate stacked bar chart
  FeaturePlot <- ggplot(df_of_genomic_features_overlapping_peaks, aes(fill = Feature, y = Frequency, x = Set)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = feature_colors) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      text = element_text(size = 8),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 8),
      legend.key.size = unit(2, 'mm'),
      panel.background = element_rect(fill = 'transparent'),
      plot.background = element_rect(fill = 'transparent', color = NA)
    ) +
    coord_flip()

  return(FeaturePlot)
}
