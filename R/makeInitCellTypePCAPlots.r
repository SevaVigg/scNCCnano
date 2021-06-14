makeInitCellTypePCAPlots <- function( seuratObj, nComps, plotDPI = 100, orientation = "portrait", name = "CellTypePCAPlot"){

# this snippet uses Seurat to make PCA plots
#
# written by Vsevolod J. Makeev 2017 - 2021

library(cowplot)
source("R/setClusterColors.r")

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

PCAPlotDir <- file.path( plotDir, "PCAPlots")
dir.create( PCAPlotDir, showWarnings = FALSE)

if (orientation == "landscape") { pageWidth = 18; pageHeight = 13 }
if (orientation == "portrait") { pageWidth = 13; pageHeight = 18 }

pca = seuratObj@dr$pca
eigValues = (pca@sdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)*100

plotList <- list()

#First draw one plot from which we get the legend
legPlot 	<- PCAPlot( seuratObj, dim.1 = 1, dim.2 = 1, do.return = TRUE, cols.use = setClusterColors( seuratObj))
plotLegend	<-  get_legend(
  # create some space to the left of the legend
  legPlot + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") 
)

maxLim <- max(apply( seuratObj@dr$pca@cell.embeddings, 2, max))
minLim <- min(apply( seuratObj@dr$pca@cell.embeddings, 2, min)) 


#then make a list of all plots without legends
for (c1 in 1:(nComps-1)){
	for (c2 in (c1+1):nComps){
	curPlot	 <- PCAPlot( seuratObj, dim.1 = c1, dim.2 = c2, pt.size = 3, do.return = TRUE, cols.use = setClusterColors( seuratObj))
	plotList <- c(plotList, list( curPlot +
		#xlim( minLim, maxLim) +
		#ylim( minLim, maxLim) + 
		theme(
			legend.position="none", 
			axis.text 	= element_text( size = 24),
			axis.text.y 	= element_text( margin = margin( r = 10)),
			axis.title  	= element_text( size = 30),
			axis.title.y 	= element_text( vjust = -3),
			axis.title.x 	= element_text( vjust = 3),
			axis.text.x	= element_text( margin = margin( t = 10)),
			panel.background = element_rect(fill = "gray60"),
			plot.margin = unit(c(.3, .3, .3, .3), "inches")
		) + 
		scale_x_continuous( 	name = paste0( "\nPC", c1, " : ", format( varExplained[c1], digits = 0), "% of variance"), 
					limits = c(minLim, maxLim), 
					expand = expand_scale(mult = c(0, .2))) +
		scale_y_continuous( 	name = paste0( "PC", c2, " : ", format( varExplained[c2], digits = 0), "% of variance\n"), 
					limits = c(minLim, maxLim),
					expand = expand_scale(mult = c(0, .2)))
	))	
	}
}

#and finally add the legend

nc 	<- if( nComps %% 2) (nComps-1)/2 else nComps/2
nr	<- if( nComps %% 2) nComps else nComps-1
pcaGrid <- plot_grid( plotlist = plotList, nrow = nr, ncol = nc, labels = "AUTO", label_size = 30) + theme( plot.margin = unit( c(1, 1, 1, 1), "inches"))

PCAplots <- plot_grid( pcaGrid) # without any legend
						
if( plotDPI == 600){
	
  PCAPlotDir600dpi <- file.path( PCAPlotDir, "600dpi")
  dir.create( PCAPlotDir600dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = PCAPlotDir600dpi, device = "png" , plot = PCAplots, width = pageWidth, height = pageHeight, units = "cm", dpi = 600, scale = 4)

}else if( plotDPI == 100){

  PCAPlotDir100dpi <- file.path( PCAPlotDir, "100dpi")
  dir.create( PCAPlotDir100dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = PCAPlotDir100dpi, device = "png" , plot = PCAplots, width = pageWidth, height = pageHeight, units = "cm", dpi = 100, scale = 4)

}else{ cat("Select plotDPI 600 or plotDPI 100")}

} #makeInitCellTypePCAPlots





