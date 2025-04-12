#' Create an Enriched Heatmap from deepTools Matrix
#'
#' Generates a composite grid of enriched heatmaps from a deepTools matrix file using a specified order
#' and plot titles. Each heatmap is drawn using consistent settings for color scale, sorting, and layout.
#'
#' @param matrix_file_path Path to a gzipped matrix file output by deepTools (e.g., from `computeMatrix`).
#' @param matrix_order A numeric vector indicating the order in which to plot matrices (e.g., `c(1, 3, 2, 4)`).
#' @param column_titles A character vector of column titles to display above each heatmap (e.g., `c("Sample A", "Sample B")`).
#'
#' @return A combined heatmap and legend grid object rendered with `cowplot::plot_grid()`.
#' @export
#'
#' @import profileplyr
#' @import EnrichedHeatmap
#' @import ComplexHeatmap
#' @import circlize
#' @import grid
#' @import gridGraphics
#' @import cowplot
#' @import magrittr
#' @importFrom SummarizedExperiment assays
make_enriched_heatmap_grid <- function(
  matrix_file_path,
  matrix_order,
  column_titles
) {
  # Import the input file
  proplyr_object <- import_deepToolsMat(matrix_file_path)

  # Extract the central 10% of columns for sorting
  sorting_vector <- lapply(assays(proplyr_object), function(mx) {
    left_column <- round(0.45 * ncol(mx))
    right_column <- round(0.55 * ncol(mx))
    central_columns <- mx[, left_column:right_column]
    z_scores <- scale(rowSums(central_columns))
    return(z_scores)
  }) %>%
    do.call(cbind, .) %>%
    rowMeans() %>%
    order() %>%
    rev()

  # Convert data to an enriched heatmap matrix
  matrix_list <- profileplyr::convertToEnrichedHeatmapMat(proplyr_object)

  # Calculate color range for the heatmap
  color_range <- matrix_list %>%
    do.call(rbind, .) %>%
    sample(length(.) * 0.1) %>%
    quantile(c(0.05, 0.95)) %>%
    as.numeric() %>%
    { c(.[1], mean(.), .[2]) }

  # Define heatmap constructor
  make_heatmap <- function(matrix_data, column_title, axis_label_color = "black") {
    EnrichedHeatmap(
      matrix_data,
      pos_line = FALSE,
      border = FALSE,
      column_title = column_title,
      column_title_gp = gpar(fontsize = 8),
      axis_name = c("-3kb", "summit", "+3kb"),
      axis_name_gp = gpar(fontsize = 6, col = axis_label_color),
      axis_name_rot = 65,
      row_order = sorting_vector,
      use_raster = TRUE,
      raster_by_magick = TRUE,
      raster_magick_filter = "Lanczos",
      raster_quality = 5,
      col = colorRamp2(color_range, viridis::viridis(n = 3)),
      row_title = NULL,
      row_title_gp = gpar(col = "transparent"),
      top_annotation = NULL,
      show_heatmap_legend = FALSE
    )
  }

  # Compose the combined heatmap
  combined_map <- NULL
  for (i in seq_along(matrix_order)) {
    heatmap <- make_heatmap(matrix_list[[matrix_order[i]]], column_titles[i])
    combined_map <- if (is.null(combined_map)) heatmap else combined_map + heatmap
  }

  # Draw legend
  legend <- Legend(
    col_fun = colorRamp2(color_range, viridis::viridis(n = 3)),
    legend_height = unit(1, "in"),
    legend_width = unit(0.05, "in"),
    border = "transparent",
    labels_gp = gpar(fontsize = 8)
  ) %>% draw() %>% grid.grabExpr()

  # Draw the heatmap grid
  heatmap_grid <- draw(combined_map) %>% grid.grabExpr()

  # Plot combined heatmap and legend
  plot_grid(heatmap_grid, NULL, legend,
            nrow = 1, rel_widths = c(10, 0.5, 1))
}
