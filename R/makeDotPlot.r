makeDotPlot	<- function( seuratObj, nLines, balanced = TRUE, plotDPI, orientation = "landscape", name){

#this snippet plots seurat DotPlots in high or low quality. Written by Seva Makeev 4.05.2021

source("R/dotPlotBalanced.r")
	
if (orientation == "landscape") { pageWidth = 18; pageHeight = nLines }
if (orientation == "portrait") { pageWidth = 13; pageHeight = 18 }

WTdotPlot 	<- dotPlotBalanced(seuratObj, genes.plot = rev(rownames(seuratObj@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))

	WTdotPlot	<- WTdotPlot +
		theme(
			plot.margin	= margin(40, 40, 40, 40),
			legend.key.size = unit( 0.6, "cm"),			
			legend.position="right",
			legend.box.margin = margin( l = 40, t = 100),
			legend.margin = margin( t = 20, b = 20),  
			legend.title = element_text( size = 40),
			legend.text = element_text( size = 30), 
 
			axis.text.y = element_text( size = 32, margin = margin( r = 20)),
			axis.text.x = element_text( size = 40, angle = 90, hjust = 0.95),
			axis.title  = element_text( size = 40, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)			

if( plotDPI == 600){
	
  dotPlotDir600dpi <- file.path( dotPlotDir, "600dpi")
  dir.create( dotPlotDir600dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = dotPlotDir600dpi, device = "png" , plot = WTdotPlot, width = pageWidth, height = pageHeight, units = "cm", dpi = 600, scale = 4)

}else if( plotDPI == 100){

  dotPlotDir100dpi <- file.path( dotPlotDir, "100dpi")
  dir.create( dotPlotDir100dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = dotPlotDir100dpi, device = "png" , plot = WTdotPlot, width = pageWidth, height = pageHeight, units = "cm", dpi = 100, scale = 4)

}else{ cat("Select plotDPI 600 or plotDPI 100")}

} #makeDotPlot.r
