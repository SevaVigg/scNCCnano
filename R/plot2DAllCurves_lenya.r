plot2DAllCurves <- function( seur2D, seurHD, lineageStart = "eNCC", dims = 1:8, lineageEnds = c("M", "I", "X"), dimRed2D = "umap", dimRedHiD = "umap", distFun = cosineClusterDist){

#seuratObject must contain results of 2D dimension reduction for plotting and tSNE and UMAP 2D pre-calculated
#clIdent contains clustering of cells for making rough lineages. seurat@ident formate (factor) would work
#cluster colors are those from seur2D

require("slingshot")
library(matrixStats)
library(ggplot2)

source("R/setClusterColors.r")
source("R/getLineageCoords.r")

cells2D <- GetCellEmbeddings( seur2D, reduction.type = dimRed2D)
cellsHD <- GetCellEmbeddings( seurHD, reduction.type = dimRedHiD)[ , dims ]
#slingLins <- getLineages( cells, clIdent, start.clus = lineageStart, end.clus = lineageEnds, dist.fun = cosine_dist)
slingLinsHD <- getLineages( cellsHD, seurHD@ident, start.clus = lineageStart, end.clus = lineageEnds, dist.fun = distFun)
slingLins2D <- getLineages( cells2D, seurHD@ident, start.clus = lineageStart, end.clus = lineageEnds, dist.fun = distFun)

#now replace the lineage thee. The correct tree must be estimated in HD
slingLins2D@lineages <- slingLinsHD@lineages  
slingCurves <- getCurves( slingLins2D, extend = "n", reassign = TRUE, stretch = 0, thresh = 0.05)

lineageIds <- which(unlist(lapply( slingLins2D@lineages, function(x) tail(x, 1) %in% c("M", "I"))))

cellColors <- setClusterColors( seurHD)[ seurHD@ident]
cl <- as.data.frame(cells2D)
cl <- cbind(cl, color=cellColors, labels=seurHD@ident)
ll<-sort(as.integer(unique(as.character(cl$labels))))
ll2<-sort(unique(as.character(cl$labels)))[length(ll)+1:length(unique(as.character(cl$labels)))]
fact_ord <- c(ll, ll2[!is.na(ll2)])
cl$labels<-factor(cl$labels, levels=fact_ord)
fckg_color_vector<-unlist(unique(cl[,c("labels", "color")])["color"])
names(fckg_color_vector) <- as.character(unlist(unique(cl[,c("labels", "color")])["labels"]))
z<-ggplot()
z<-z+geom_point(data=as.data.frame(cl), aes(UMAP1, UMAP2, color=labels))+scale_color_manual(values = fckg_color_vector) + theme_bw()
z<-z+theme( axis.text.x = element_text( size = 20), axis.text.y = element_text(size = 20),
      axis.title.x = element_text( size = 20, margin = margin( t = 5, r = 0, b = 0, l = 0)), axis.title.y = element_text( size = 20),
   legend.key.size = unit( 1.2, "cm"),
legend.position="right",
legend.title =element_blank(),
legend.text = element_text( size = 18))

for (x in lineageIds) {
a<-slingCurves@curves[[x]]$s[slingCurves@curves[[x]]$ord,]
z<-z+geom_path(data=as.data.frame(a), aes(x=UMAP1, y=UMAP2))
}
z

}
