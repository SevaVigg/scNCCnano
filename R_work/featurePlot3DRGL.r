# This script requires seurat object prepared by Seurat with three axes obtained by dimernsion reduction

featurePlot3DRGL	<- function( seuratObj, gene, dimRed){

if(!require("colorspace")){
install.packages("colorspace")
library(colorspace)}

require("Seurat")
require("rgl")
require("colorspace")
source("R/setClusterColors.r")

#rgl.open()
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




