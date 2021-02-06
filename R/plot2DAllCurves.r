plot2DAllCurves	<- function( seur2D, seurHD, lineageStart = "eNCC", dims = 1:8, lineageEnds = c("M", "I", "X"), dimRed2D = "umap", dimRedHiD = "umap", distFun = cosineClusterDist){

#seuratObject must contain results of 2D dimension reduction for plotting and tSNE and UMAP 2D pre-calculated
#clIdent contains clustering of cells for making rough lineages. seurat@ident formate (factor) would work
#cluster colors are those from seur2D

require("slingshot")
library(matrixStats)

source("R/setClusterColors.r")
source("R/getLineageCoords.r")

cells2D		<- GetCellEmbeddings( seur2D, reduction.type = dimRed2D)
cellsHD		<- GetCellEmbeddings( seurHD, reduction.type = dimRedHiD)[ , dims ]
#slingLins	<- getLineages( cells, clIdent, start.clus = lineageStart, end.clus = lineageEnds, dist.fun = cosine_dist)
slingLinsHD	<- getLineages( cellsHD, seurHD@ident, start.clus = lineageStart, end.clus = lineageEnds, dist.fun = distFun)
slingLins2D	<- getLineages( cells2D, seurHD@ident, start.clus = lineageStart, end.clus = lineageEnds, dist.fun = distFun)

#now replace the lineage thee. The correct tree must be estimated in HD
slingLins2D@lineages <- slingLinsHD@lineages  
slingCurves	<- getCurves( slingLins2D, extend = "n", reassign = TRUE, stretch = 0, thresh = 0.05)

lineageIds 	<- which(unlist(lapply( slingLins2D@lineages, function(x) tail(x, 1) %in% c("M", "I")))) 

cellColors 	<- setClusterColors( seurHD)[ seurHD@ident]
plot( cells2D, col = cellColors, cex = 1, pch = 16)

lapply( lineageIds, function(x) lines( slingCurves@curves[[x]] ))

return( slingCurves)

}


