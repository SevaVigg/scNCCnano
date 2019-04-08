# This script requires seurat object prepared by Seurat with three axes obtained by dimernsion reduction

dimPlot3DFeature	<- function( seuratObj, gene, dimRed){

require("Seurat")
require("rgl")
source("R/setClusterColors.r")

open3d()

bg3d("darkblue")
title3d( main = gene, color = "white")
rgl.bringtotop()

x <- seuratObj@dr[[dimRed]]@cell.embeddings[,1]
y <- seuratObj@dr[[dimRed]]@cell.embeddings[,2]
z <- seuratObj@dr[[dimRed]]@cell.embeddings[,3]

featureValues <- seuratObj@data[ gene, ]
cellColors <- cut( featureValues/max(featureValues),
breaks= seq(0, 1, 0.1))

colors <- sequential_hcl(11, palette = "YlGnBu")
levels(cellColors) <- colors
levels(cellColors) <- c(levels(cellColors), "cyan", "black", "red")

rgl.points(x,y,z, color = cellColors)
}



