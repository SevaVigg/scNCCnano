#function clDotPlotLenya

source("R/dotPlotLenya.r")


clDotPlot <- dotPlotLenya( valUmap, genes.plot = rev(rownames(
valUmap@data)), x.lab.rot = TRUE, dot.scale = 10,
plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return =
TRUE, cols.use = c("cyan", "red"))
clDotPlot <- clDotPlot +
theme(
legend.position="right",
legend.title = element_text( size = 24),
legend.text = element_text( size = 14),

axis.text.y = element_text( size = 24),
axis.text.x = element_text( size = 24, angle = 90),
axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
)

