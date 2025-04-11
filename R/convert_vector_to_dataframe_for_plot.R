#' Convert a Numeric Vector to a Data Frame for Plotting
#'
#' This function converts a numeric vector (e.g., peak counts) into a data frame
#' that can be used for visualization.
#'
#' @param count_vector Numeric vector representing counts (e.g., peak counts per bin).
#'
#' @return A data frame with two columns:
#'   - `counts`: The peak count per bin.
#'   - `x`: The bin index (1-based).
#' @family noble_peak_density_heatmap
#' @export
#' @examples
#' convert_vector_to_dataframe_for_plot(c(1, 3, 5, 2, 0))
convert_vector_to_dataframe_for_plot <- function(count_vector){
  count_data_frame <- data.frame(counts = count_vector, x = 1:length(count_vector))
  return(count_data_frame)
}