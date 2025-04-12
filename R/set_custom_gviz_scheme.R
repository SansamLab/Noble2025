#' Set a Custom Gviz Scheme for Plotting
#'
#' Applies a consistent visual theme to Gviz genome tracks by modifying font sizes, colors, and other style attributes.
#'
#' @param cex A numeric value for scaling font sizes and line widths (default = 1).
#' @param scheme_name A name for the scheme to register (default = "myScheme").
#'
#' @return Invisibly returns the modified scheme list. Registers scheme globally via `options(Gviz.scheme)`.
#' @export
#'
#' @import Gviz
set_custom_gviz_scheme <- function(cex = 1, scheme_name = "myScheme") {
  custom_scheme <- Gviz::getScheme()

  # GdObject
  custom_scheme$GdObject$col <- "black"
  custom_scheme$GdObject$fontsize <- 8
  custom_scheme$GdObject$cex.axis <- cex
  custom_scheme$GdObject$cex.title <- cex
  custom_scheme$GdObject$col.axis <- "black"
  custom_scheme$GdObject$col.border.title <- "transparent"
  custom_scheme$GdObject$col.frame <- "transparent"
  custom_scheme$GdObject$col.title <- "black"
  custom_scheme$GdObject$col.id <- "black"
  custom_scheme$GdObject$background.title <- "transparent"
  custom_scheme$GdObject$fontcolor <- "black"
  custom_scheme$GdObject$fontcolor.title <- "black"
  custom_scheme$GdObject$rotation.title <- 0
  custom_scheme$GdObject$rotation <- 90

  # DataTrack
  custom_scheme$DataTrack$fontsize.legend <- 8
  custom_scheme$DataTrack$cex.legend <- cex
  custom_scheme$DataTrack$cex.sampleNames <- cex
  custom_scheme$DataTrack$cex <- cex
  custom_scheme$DataTrack$cex.axis <- cex
  custom_scheme$DataTrack$cex.title <- cex
  custom_scheme$DataTrack$fontcolor.legend <- "black"
  custom_scheme$DataTrack$col.sampleNames <- "black"

  # AnnotationTrack
  custom_scheme$AnnotationTrack$fontsize.group <- 8
  custom_scheme$AnnotationTrack$cex <- cex
  custom_scheme$AnnotationTrack$cex.group <- cex
  custom_scheme$AnnotationTrack$rotation <- 90
  custom_scheme$AnnotationTrack$rotation.item <- 90
  custom_scheme$AnnotationTrack$rotation.group <- 90
  custom_scheme$AnnotationTrack$fontcolor.group <- "black"
  custom_scheme$AnnotationTrack$fontcolor.item <- "black"
  custom_scheme$AnnotationTrack$background.title <- "transparent"

  # GenomeAxisTrack
  custom_scheme$GenomeAxisTrack$cex.id <- cex
  custom_scheme$GenomeAxisTrack$cex <- cex
  custom_scheme$GenomeAxisTrack$col <- "black"
  custom_scheme$GenomeAxisTrack$fontcolor <- "black"
  custom_scheme$GenomeAxisTrack$size <- 3
  custom_scheme$GenomeAxisTrack$fontsize <- 6
  custom_scheme$GenomeAxisTrack$cex.title <- cex

  # IdeogramTrack
  custom_scheme$IdeogramTrack$fontsize <- 8
  custom_scheme$IdeogramTrack$cex.bands <- cex
  custom_scheme$IdeogramTrack$cex <- cex
  custom_scheme$IdeogramTrack$fontcolor <- "black"
  custom_scheme$IdeogramTrack$size <- 4
  custom_scheme$IdeogramTrack$cex.title <- cex
  custom_scheme$IdeogramTrack$showTitle <- FALSE
  custom_scheme$IdeogramTrack$showId <- FALSE
  custom_scheme$IdeogramTrack$lwd <- 3  # reasonable default, can be overridden

  # GeneRegionTrack
  custom_scheme$GeneRegionTrack$background.title <- "transparent"


  Gviz::addScheme(custom_scheme, scheme_name)
  options(Gviz.scheme = scheme_name)

  invisible(custom_scheme)
}
