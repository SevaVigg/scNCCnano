makeInitCellTypePCAPlots <- function( seuratObj, Ncomps){

library(cowplot)
source("R/setClusterColors.r")

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
for (c1 in 1:(Ncomps-1)){
	for (c2 in (c1+1):Ncomps){
	curPlot	 <- PCAPlot( seuratObj, dim.1 = c1, dim.2 = c2, pt.size = 2, do.return = TRUE, cols.use = setClusterColors( seuratObj))
	plotList <- c(plotList, list( curPlot +
		xlim( minLim, maxLim) +
		ylim( minLim, maxLim) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			panel.background = element_rect(fill = "gray60")
		)
	))	
	}
}

#and finally add the legend

nc 	<- if( Ncomps %% 2) (Ncomps-1)/2 else Ncomps/2
nr	<- if( Ncomps %% 2) Ncomps else Ncomps-1
pcaGrid <- plot_grid( plotlist = plotList, nrow = nr, ncol = nc, labels = "AUTO", label_size = 35)

PCAplots <- plot_grid( pcaGrid)							# without any legend
return(PCAplots)

}


