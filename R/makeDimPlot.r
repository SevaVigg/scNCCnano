makeDimPlot	<- function( seuratObj, dimRed = "umap",  name = "clusterPlot", reorderClusters = FALSE, plotDPI = 100, orientation = "landscape"){

source("R/setClusterColors.r")

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

clusterPlotDir <- file.path( plotDir, "clusterPlots")
dir.create( clusterPlotDir, showWarnings = FALSE)


if (orientation == "landscape") { pageWidth = 18; pageHeight = 13 }
if (orientation == "portrait") { pageWidth = 18; pageHeight = 13 }


if (reorderClusters) {
	orderedCells	<- ordered( seuratObj@ident, levels = c( "eHMP", "ltHMP", "I", "M", "X", "7", "4"))
	orderedCells	<- orderedCells[ order( orderedCells)]
	seuratObj@ident <- orderedCells}

col = setClusterColors( seuratObj)


dimPlot <-  DimPlot( seuratObj, reduction.use = dimRed, cols.use = col, pt.size = 5) +
		labs( color = "Cell type\n") + 
		theme( 	
			plot.margin  = margin( t = 40, r = 40, b = 40, l = 40), 

			axis.title   = element_text( size = 40, face = "bold"),
			axis.text.x  = element_text( size = 40, margin = margin( t = 20)), 
			axis.text.y  = element_text( size = 40, margin = margin( r = 20)),
			axis.title.x = element_text( size = 40, margin = margin( t = 20, r = 0, b = 0, l = 0)), 
			axis.title.y = element_text( size = 40, margin = margin( r = 20)),
			
			legend.key.size = unit( 0.6, "cm"),			
			legend.position="right",
			legend.box.margin = margin( l = 40, t = 80),
			legend.title = element_text( size = 40),
			legend.text  = element_text( size = 40, margin = margin( b = 20))
			
			)

if( plotDPI == 600){
	
  clusterPlotDir600dpi <- file.path( clusterPlotDir, "600dpi")
  dir.create( clusterPlotDir600dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = clusterPlotDir600dpi, device = "png" , plot = dimPlot, width = pageWidth, height = pageHeight, units = "cm", dpi = 600, scale = 4)

}else if( plotDPI == 100){

  clusterPlotDir100dpi <- file.path( clusterPlotDir, "100dpi")
  dir.create( clusterPlotDir100dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = clusterPlotDir100dpi, device = "png" , plot = dimPlot, width = pageWidth, height = pageHeight, units = "cm", dpi = 100, scale = 4)

}else{ cat("Select plotDPI 600 or plotDPI 100")}




}
