# Plots colored clusters, trajectories and the legend. Written by Leonid Uroshlev 07.02.2021

library(Seurat)
plot2DAllCurves <- function( seur2D, seurHD, lineageStart = "eHMP",
dims = 1:8, genes.use = NULL, lineageEnds = c("M", "I", "X"), dimRed2D
= "umap", dimRedHiD = "umap", distFun = cosineClusterDist, curveThresh  = 0.95){

require("slingshot")
library(matrixStats)
library(ggplot2)
if(!require(gtools)){install.packages("gtools")}
library(gtools)

source("R/setClusterColors.r")
source("R/getLineageCoords.r")

cells2D 	<- GetCellEmbeddings( seur2D, reduction.type = dimRed2D)
dimNames 	<- attr( cells2D, "dimnames")[[2]]

if ( !is.null( genes.use)) {
cellsHD <- t( as.matrix(seurHD@data[ genes.use, ]))
}else{ cellsHD <- GetCellEmbeddings( seurHD, reduction.type = dimRedHiD)[ , dims ]}

slingLinsHD <- getLineages( cellsHD, seurHD@ident, start.clus = lineageStart, end.clus = lineageEnds, dist.fun = distFun)
slingLins2D <- getLineages( cells2D, seurHD@ident, start.clus = lineageStart, end.clus = lineageEnds, dist.fun = distFun)

#now replace the lineage thee. The correct tree must be estimated in HD
slingLins2D@lineages <- slingLinsHD@lineages
slingCurves <- getCurves( slingLins2D, extend = "n", reassign = FALSE,
stretch = 0, thresh = 0.005, shrink = 0.2)

lineageIds <- which(unlist(lapply( slingLins2D@lineages, function(x)
tail(x, 1) %in% c("M", "I"))))

cellColors <- setClusterColors( seurHD)[ seurHD@ident]
cl <- as.data.frame(cells2D)
cl <- cbind(cl, color=cellColors, labels=seurHD@ident)
cl$labels<-factor(cl$labels, levels=unique(mixedsort(as.character(cl$label))))
fckg_color_vector<-as.character(unlist(unique(cl[,c("labels", "color")])["color"]))
names(fckg_color_vector) <- as.character(unlist(unique(cl[,c("labels","color")])["labels"]))
z<-ggplot()
z<-z+geom_point(data=as.data.frame(cl), aes( get(dimNames[1]), get(dimNames[2]), color=labels))+scale_color_manual(values = fckg_color_vector)
z<-z+theme( axis.text.x = element_text( size = 20), axis.text.y = element_text(size = 20), axis.title.x = element_text( size = 20, margin = margin( t = 5, r = 0, b = 0, l = 0)), axis.title.y = element_text( size = 20), legend.key.size = unit( 1.2, "cm"), legend.position="right",
legend.title =element_blank(),
legend.text = element_text( size = 18)) + xlab(dimNames[1]) + ylab(dimNames[2])

for (x in lineageIds) {
curve 		<- slingCurves@curves[[x]]
curveOrd 	<- curve$s[ curve$ord, ]
curveWt		<- curve$w[ curve$ord ]
curveData	<- curveOrd[ curveWt > curveThresh, ]
z<-z+geom_path(data=as.data.frame( curveData), aes(x = get(dimNames[1]), y = get(dimNames[2])))
}
z

}
