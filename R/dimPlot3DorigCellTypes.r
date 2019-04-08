# This script requires seurat object prepared by Seurat with three axes obtained by dimernsion reduction

dimPlot3DorigCellTypes	<- function( seuratObj, dimRed){

require("Seurat")
require("rgl")
source("R/setClusterColors.r")
source("R/setCellTypeColors.r")

open3d()
bg3d("darkblue")

x <- seuratObj@dr[[dimRed]]@cell.embeddings[,1]
y <- seuratObj@dr[[dimRed]]@cell.embeddings[,2]
z <- seuratObj@dr[[dimRed]]@cell.embeddings[,3]


cellColors <- setCellTypeColors( seuratObj)[ seuratObj@ident]

rgl.points(x,y,z, color = cellColors)
}




