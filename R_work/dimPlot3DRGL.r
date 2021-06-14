# This script requires seurat object prepared by Seurat with three axes obtained by dimernsion reduction

dimPlot3DRGL	<- function( seuratObj, dimRed){

require("Seurat")
require("rgl")
source("R/setClusterColors.r")

open3d()
bg3d("grey")

x <- seuratObj@dr[[dimRed]]@cell.embeddings[,1]
y <- seuratObj@dr[[dimRed]]@cell.embeddings[,2]
z <- seuratObj@dr[[dimRed]]@cell.embeddings[,3]


cellColors <- setClusterColors( seuratObj)[ seuratObj@ident]

rgl.points(x,y,z, color = cellColors)
}




