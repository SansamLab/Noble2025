#' custom ggplot theme for Sansam plots

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