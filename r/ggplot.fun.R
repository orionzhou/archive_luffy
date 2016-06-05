require(plyr)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(gtable)
require(gridBase)

get_ggplot_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}