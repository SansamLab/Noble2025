#' Sansam Lab Default ggplot2 Theme
#'
#' A minimalist `ggplot2` theme with small font sizes and transparent backgrounds,
#' designed for compact scientific figures.
#'
#' @return A `ggplot2::theme` object with customized appearance:
#' - Transparent background
#' - Small font sizes (6–8 pt)
#' - No outer margins
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   sansam_theme()
#'
#' @import ggplot2
#' @export
sansam_theme <- function(){
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    plot.margin=ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    text = ggplot2::element_text(size = 8),
    axis.text = ggplot2::element_text(size=8),
    axis.text.x=ggplot2::element_text(size=6),
    axis.text.y=ggplot2::element_text(size=6)
  )
}

#' Sansam Lab Theme with Axis Border and Title Styling
#'
#' A `ggplot2` theme similar to `sansam_theme()` but includes a black panel border
#' and smaller title text. Ideal for publication-quality plots with boxed panels.
#'
#' @return A `ggplot2::theme` object with:
#' - Transparent backgrounds
#' - Panel border (black)
#' - Small font sizes (6–8 pt)
#' - No outer margins
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   sansam_theme2()
#'
#' @import ggplot2
#' @export
sansam_theme2 <- function(){
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    plot.margin=ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    text = ggplot2::element_text(size = 8),
    axis.text = ggplot2::element_text(size=8),
    plot.title = ggplot2::element_text(size=8),
    axis.text.x=ggplot2::element_text(size=6),
    axis.text.y=ggplot2::element_text(size=6),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
}
