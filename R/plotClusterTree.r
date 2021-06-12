plotClusterTree	<- function( seuratObj, plotDPI = "100", treeName){

require( ape)

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

clusterTreeDir		<- file.path( plotDir, "clusterTrees")
dir.create( clusterTreeDir, showWarnings = FALSE)

clusterTreeDir100dpi	<- file.path( clusterTreeDir, "100dpi")
dir.create( clusterTreeDir100dpi, showWarnings = FALSE)

clusterTreeDir600dpi	<- file.path( clusterTreeDir, "600dpi")
dir.create( clusterTreeDir600dpi, showWarnings = FALSE)

if (plotDPI == 600) {
png( file = file.path( clusterTreeDir600dpi, paste0( treeName, ".png")),
	height = 5,
	width =  5,
	units = "in",
	res = 600, 
	pointsize = 4 
)
	PlotClusterTree( seuratObj,
	type	= "phylogram",
	cex 	= 4,
	adj	= 0,
	label.offset = 0.005
	)	
 	nodelabels( text = " ")
dev.off()
}else if (plotDPI == 100) {
png( file = file.path( clusterTreeDir100dpi, paste0( treeName, ".png")),
	height = 5,
	width =  5,
	units = "in",
	res = 100, 
	pointsize = 3
)
	PlotClusterTree( seuratObj, 
	type	= "phylogram",
	cex 	= 4,
	adj	= 0,
	label.offset = 0.005
	)	
 	nodelabels( text = " ")
dev.off()
}else{ cat("Select plotDPI 600 or plotDPI 100")}
}
