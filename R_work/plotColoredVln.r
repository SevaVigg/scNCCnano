plotColoredVln <- function( seuratObj, gene){

source("R/setClusterColors.r")
vln <- VlnPlot(seuratObj, gene, do.return = TRUE)

plot(vln 
#	+ scale_x_discrete(limits = slingObj@lineages[[lineageId]])
	+ scale_fill_manual(values = setClusterColors( seuratObj))
)
}
