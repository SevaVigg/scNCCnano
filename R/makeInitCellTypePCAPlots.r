makeInitCellTypePCAPlots <- function( seuratObj, Ncomps, doReturn, doPrint = TRUE){

library(cowplot)
source("R/setClusterColors.r")

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
for (c1 in 1:(Ncomps-1)){
	for (c2 in (c1+1):Ncomps){
	curPlot	 <- PCAPlot( seuratObj, dim.1 = c1, dim.2 = c2, pt.size = 3, do.return = TRUE, cols.use = setClusterColors( seuratObj))
	plotList <- c(plotList, list( curPlot +
		#xlim( minLim, maxLim) +
		#ylim( minLim, maxLim) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 24),
			axis.title  = element_text( size = 30),
			panel.background = element_rect(fill = "gray60"),
			plot.margin = unit(c(.4, .4, .4, .4), "inches")
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

nc 	<- if( Ncomps %% 2) (Ncomps-1)/2 else Ncomps/2
nr	<- if( Ncomps %% 2) Ncomps else Ncomps-1
pcaGrid <- plot_grid( plotlist = plotList, nrow = nr, ncol = nc, labels = "AUTO", label_size = 40) + theme( plot.margin = unit( c(1, 1, 1, 1), "inches"))

PCAplots <- plot_grid( pcaGrid) # without any legend

if( doReturn) return(PCAplots)

}


