createSlingShotObject <- function( seuratObj, dimRed, genes.use = NULL, useDims = 1:10, startClust = "eHMP", endClust = c("I", "M"), distFun = cosineClusterDist){

if(!require(slingshot)){
  install.packages("slingshot")
}

library("slingshot")	
library("matrixStats")

source("R/cosineClusterDist.r")
source("R/getFinalClusterTypes.r")

if ( !is.null( genes.use)) {

coordMatrixMD 	<- t( as.matrix(seuratObj@data[ genes.use, ]))

}else{ coordMatrixMD  <-  as.matrix(seuratObj@dr[[dimRed]]@cell.embeddings)[ , useDims]}

clTypes		<- getFinalClusterTypes( seuratObj)
slingShotObj 	<- slingshot(coordMatrixMD, seuratObj@ident, start.clus = startClust ,end.clus = endClust, reassign = TRUE, dist.fun = distFun )

}


