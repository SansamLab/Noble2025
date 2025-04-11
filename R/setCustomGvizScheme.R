
#' Apply Custom Gviz Scheme for Consistent Track Styling
#'
#' This function defines and registers a customized Gviz visual scheme called `"myScheme"`,
#' which sets consistent font sizes, colors, and graphical scaling for genome tracks.
#' It then sets it as the default Gviz scheme via `options(Gviz.scheme = "myScheme")`.
#'
#' @return Invisibly returns the name of the applied scheme.
#' @import Gviz
#' @export
setCustomGvizScheme <- function() {
  scheme <- Gviz::getScheme()
  cex <- 1

  scheme$GdObject$col <- "black"
  scheme$GdObject$fontsize <- 8
  scheme$DataTrack$fontsize.legend <- 8
  scheme$IdeogramTrack$fontsize <- 8
  scheme$AnnotationTrack$fontsize.group <- 8

  scheme$GenomeAxisTrack$cex.id <- cex
  scheme$GenomeAxisTrack$cex <- cex
  scheme$GdObject$cex.axis <- cex
  scheme$GdObject$cex.title <- cex
  scheme$DataTrack$cex.legend <- cex
  scheme$DataTrack$cex.sampleNames <- cex
  scheme$DataTrack$cex <- cex
  scheme$DataTrack$cex.axis <- cex
  scheme$DataTrack$cex.title <- cex
  scheme$IdeogramTrack$cex.bands <- cex
  scheme$IdeogramTrack$cex <- cex
  scheme$AnnotationTrack$cex <- cex
  scheme$AnnotationTrack$cex.group <- cex

  scheme$GenomeAxisTrack$fontcolor <- "black"
  scheme$DataTrack$fontcolor.legend <- "black"
  scheme$GdObject$fontcolor <- "black"
  scheme$IdeogramTrack$fontcolor <- "black"
  scheme$AnnotationTrack$fontcolor.group <- "black"
  scheme$AnnotationTrack$fontcolor.item <- "black"

  scheme$GdObject$col.axis <- "black"
  scheme$GdObject$col.border.title <- "transparent"
  scheme$GdObject$col.frame <- "transparent"
  scheme$GdObject$col.title <- "black"
  scheme$GdObject$col.id <- "black"
  scheme$GenomeAxisTrack$col <- "black"
  scheme$DataTrack$col.sampleNames <- "black"
  scheme$GdObject$background.title <- "white"
  scheme$GdObject$fontcolor.title <- "black"

  scheme$AnnotationTrack$rotation <- 90
  scheme$AnnotationTrack$rotation.item <- 90
  scheme$AnnotationTrack$rotation.group <- 90
  scheme$GdObject$rotation.title <- 0
  scheme$GdObject$rotation <- 90

  Gviz::addScheme(scheme, "myScheme")
  options(Gviz.scheme = "myScheme")

  invisible("myScheme")
}
