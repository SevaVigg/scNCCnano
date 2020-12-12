plot2DAllCurves	<- function( seur2D, clIdent, lineageStart = "eNCC", lineageEnds = c("M", "I", "X"), dimRed){

#seuratObject must contain results of 2D dimension reduction for plotting and tSNE and UMAP 2D pre-calculated
#clIdent contains clustering of cells for making rough lineages. seurat@ident formate (factor) would work
#cluster colors are those from seur2D

require("slingshot")
library(matrixStats)

source("R/setClusterColors.r")
source("R/getLineageCoords.r")

cells		<- GetCellEmbeddings( seur2D, reduction.type = dimRed)
#slingLins	<- getLineages( cells, clIdent, start.clus = lineageStart, end.clus = lineageEnds, dist.fun = cosine_dist)
slingLins	<- getLineages( cells, clIdent, start.clus = lineageStart, end.clus = lineageEnds)
slingCurves	<- getCurves( slingLins, extend = "n", reassign = TRUE)

lineageIds 	<- which(unlist(lapply( slingLins@lineages, function(x) tail(x, 1) %in% lineageEnds))) 

cellColors 	<- setClusterColors( seur2D)[ seur2D@ident]
plot( cells, col = cellColors, cex = 1, pch = 16)

lapply( lineageIds, function(x) lines( slingCurves@curves[[x]] ))

return( slingCurves)

}


