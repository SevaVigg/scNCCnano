#This function constructs DotPlot with cluster coloring taken over all cell types. Written by Leonid Uroshlev, 21.01.2021

if(!require( ggplot2)){
  install.packages("ggplot2")
}
library("ggplot2")	

if(!require( dplyr)){
  install.packages("dplyr")
}
library("dplyr")	

if(!require( tidyr)){
  install.packages("tidyr")
}
library("tidyr")	

library("Seurat")

PercentAbove<- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}


dotPlotBalanced <- function(
  object,
  genes.plot,
  cols.use = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by,
  plot.legend = FALSE,
  do.return = FALSE,
  x.lab.rot = FALSE,
  scale.by = "radius"
) {
  if (! missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  data.to.plot %>% gather(
    key = genes.plot,
    value = expression,
    -c(cell, id)
  ) -> data.to.plot
  data.to.plot %>%
    group_by(id, genes.plot) %>%
    summarize(
      avg.exp = mean(expm1(x = expression)),
      pct.exp = PercentAbove(x = expression, threshold = 0)
    ) -> data.to.plot
  data.to.plot %>%
    ungroup() %>%
    mutate(avg.exp.scale = scale(x = avg.exp)) %>%
    mutate(avg.exp.scale = MinMax(
      data = avg.exp.scale,
      max = col.max,
      min = col.min
    )) ->  data.to.plot
  data.to.plot$genes.plot <- factor(
    x = data.to.plot$genes.plot,
    levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
  )
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
    geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
    scale_radius(range = c(0, dot.scale)) +
    scale_color_gradient(low = cols.use[1], high = cols.use[2]) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  if (! plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}


