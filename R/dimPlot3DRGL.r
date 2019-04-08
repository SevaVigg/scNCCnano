# This script requires seurat object prepared by Seurat with three axes obtained by dimernsion reduction

dimPlot3DRGL	<- function( seuratObj, dimRed){

require("Seurat")
require("rgl")
source("R/setClusterColors.r")

open3d()
bg3d("darkblue")

x <- seuratObj@dr[[dimRed]]@cell.embeddings[,1]
y <- seuratObj@dr[[dimRed]]@cell.embeddings[,2]
z <- seuratObj@dr[[dimRed]]@cell.embeddings[,3]


cellColors <- setClusterColors( seuratObj)[ seuratObj@ident]
names(cellColors) <- levels( seuratObj@ident)

cellColors[ "I" ] 	<- "cyan"
cellColors[ "M" ]   	<- "black"
cellColors[ "E" ]  	<- "red"

names(cellColors) <- NULL

rgl.points(x,y,z, color = cellColors)
}




