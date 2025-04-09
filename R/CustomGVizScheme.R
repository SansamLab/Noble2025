library(Gviz)
custom_cex <- 1
customGvizScheme <- Gviz::getScheme()
customGvizScheme$GdObject$col <- "black"
customGvizScheme$GdObject$fontsize <- 8
customGvizScheme$DataTrack$fontsize.legend <- 8
customGvizScheme$IdeogramTrack$fontsize <- 8
customGvizScheme$AnnotationTrack$fontsize.group <- 8
customGvizScheme$GenomeAxisTrack$cex.id=custom_cex
customGvizScheme$GenomeAxisTrack$cex=custom_cex
customGvizScheme$GdObject$cex.axis=custom_cex
customGvizScheme$GdObject$cex.title=custom_cex
customGvizScheme$DataTrack$cex.legend=custom_cex
customGvizScheme$DataTrack$cex.sampleNames=custom_cex
customGvizScheme$DataTrack$cex=custom_cex
customGvizScheme$DataTrack$cex.axis=custom_cex
customGvizScheme$DataTrack$cex.title=custom_cex
customGvizScheme$IdeogramTrack$cex.bands=custom_cex
customGvizScheme$IdeogramTrack$cex=custom_cex
customGvizScheme$AnnotationTrack$cex=custom_cex
customGvizScheme$AnnotationTrack$cex.group=custom_cex
customGvizScheme$GenomeAxisTrack$fontcolor <- "black"
customGvizScheme$DataTrack$fontcolor.legend <- "black"
customGvizScheme$GdObject$fontcolor <- "black"
customGvizScheme$IdeogramTrack$fontcolor <- "black"
customGvizScheme$AnnotationTrack$fontcolor.group <- "black"
customGvizScheme$AnnotationTrack$fontcolor.item <- "black"
customGvizScheme$GdObject$col.axis="black"
customGvizScheme$GdObject$col.border.title="black"
customGvizScheme$GdObject$col.frame="transparent"
customGvizScheme$GdObject$col.title="black"
customGvizScheme$GdObject$col.id="black"
customGvizScheme$GenomeAxisTrack$col="black"
customGvizScheme$DataTrack$col.sampleNames="black"
customGvizScheme$GdObject$col.frame="transparent"
customGvizScheme$DataTrack$col.sampleNames="black"
customGvizScheme$GdObject$background.title="white"
customGvizScheme$GdObject$fontcolor <- "black"
customGvizScheme$AnnotationTrack$fontcolor.item="black"
customGvizScheme$GdObject$col.border.title="transparent"
customGvizScheme$GdObject$fontcolor.title="black"
customGvizScheme$AnnotationTrack$rotation=90
customGvizScheme$AnnotationTrack$rotation.item=90
customGvizScheme$AnnotationTrack$rotation.group=90
customGvizScheme$GdObject$rotation.title=0
customGvizScheme$GdObject$rotation=90
Gviz::addScheme(customGvizScheme, "myScheme")
options(Gviz.scheme = "myScheme")