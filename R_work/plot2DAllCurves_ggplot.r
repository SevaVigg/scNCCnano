# Plots colored clusters, trajectories and the legend. Written by Leonid Uroshlev 07.02.2021
# Modified by Vsevolod Makeev 09.05.2021

library(Seurat)
plot2DAllCurves <- function( 	seur2D, 
				slingObjs, 
				dimRed2D = "umap", 
				lineageToDrawEnds = c("M", "I"),
				cellsKeepThresh = 0.95,
				plotDPI = 100,
				name = "slingshotClusterPlot", 
				orientation = "landscape"){

require("slingshot")
library(matrixStats)
library(ggplot2)
if(!require(gtools)){install.packages("gtools")}
library(gtools)

source("R/setClusterColors.r")
source("R/getLineageCoords.r")

resDir		<- file.path(getwd(), "Res")
plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

clusterPlotDir 	<- file.path( plotDir, "clusterPlots")
dir.create( clusterPlotDir, showWarnings = FALSE)
 
if (orientation == "landscape") { pageWidth = 11.7; pageHeight = 8.3 }
if (orientation == "portrait") { pageWidth = 8.3; pageHeight = 11.7 }

lineageIds <- which(unlist(lapply( slingObjs$Dim2@lineages, function(x)
tail(x, 1) %in% lineageToDrawEnds)))
#slingshot to here

cells2D 	<- GetCellEmbeddings( seur2D, reduction.type = dimRed2D)
dimNames 	<- attr( cells2D, "dimnames")[[2]]

cellColors <- setClusterColors( seur2D)[ seur2D@ident]
cl <- as.data.frame(cells2D)
cl <- cbind(cl, color=cellColors, labels=seur2D@ident)
cl$labels<-factor(cl$labels, levels=unique(mixedsort(as.character(cl$label))))
fckg_color_vector<-as.character(unlist(unique(cl[,c("labels", "color")])["color"]))
names(fckg_color_vector) <- as.character(unlist(unique(cl[,c("labels","color")])["labels"]))

z<-ggplot()
z<-z+	geom_point(data = cl, aes( get(dimNames[1]), get(dimNames[2]), color=labels), size = 3)+scale_color_manual(values = fckg_color_vector)
z<-z+	labs( color = "Cell type\n", x = colnames(cl)[1], y = colnames(cl)[2]) + 
	theme(
			plot.margin  = margin( t = 40, r = 40, b = 40, l = 40), 

			axis.title   = element_text( size = 32, face = "bold"),
			axis.text.x  = element_text( size = 32, margin = margin( t = 20)), 
			axis.text.y  = element_text( size = 32, margin = margin( r = 20)),
			axis.title.x = element_text( size = 32, margin = margin( t = 20, r = 0, b = 0, l = 0)), 
			axis.title.y = element_text( size = 32, margin = margin( r = 20)),
			
			legend.key.size = unit( 1, "cm"),			
			legend.position="right",
			legend.box.margin = margin( l = 20, t = 80),
			legend.title = element_text( size = 32),
			legend.text  = element_text( size = 32, margin = margin( b = 10))) +
			
			guides(color = guide_legend(override.aes = list(size = 10)))

for (x in lineageIds) {
curve 		<- slingObjs$Dim2@curves[[x]]
curveOrd 	<- curve$s[ curve$ord, ]
curveWt		<- curve$w[ curve$ord ]
curveData	<- curveOrd[ curveWt > cellsKeepThresh, ]
z<-z+geom_path(data=as.data.frame( curveData), aes(x = get(dimNames[1]), y = get(dimNames[2])))
}

if( plotDPI == 600){
	
  clusterPlotDir600dpi <- file.path( clusterPlotDir, "600dpi")
  dir.create( clusterPlotDir600dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = clusterPlotDir600dpi, device = "png" , plot = z, width = pageWidth, height = pageHeight, units = "cm", dpi = 600, scale = 4)

}else if( plotDPI == 100){

  clusterPlotDir100dpi <- file.path( clusterPlotDir, "100dpi")
  dir.create( clusterPlotDir100dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = clusterPlotDir100dpi, device = "png" , plot = z, width = pageWidth, height = pageHeight, units = "cm", dpi = 100, scale = 4)

}else{ cat("Select plotDPI 600 or plotDPI 100")}



z

}
