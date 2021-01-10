createSlingShotObject_Lenya <- function( seuratObj, dimRed, startClust = "eNCC", endClust = c("I", "M", "X")){

if(!require(slingshot)){
  install.packages("slingshot")
}

library("slingshot")	
library("matrixStats")

cosine_cluster <- function(X, w1, w2){
    mu1 <- colWeightedMeans(X, w = w1)
    mu2 <- colWeightedMeans(X, w = w2)
    return(sum(mu1*mu2)/sqrt(sum(mu1^2)*sum(mu2^2)))
}

source("R/getFinalClusterTypes.r")

coordMatrixMD  <-  as.matrix(seuratObj@dr[[dimRed]]@cell.embeddings)

clTypes		<- getFinalClusterTypes( seuratObj)
slingShotObj 	<- slingshot(coordMatrixMD, seuratObj@ident, start.clus = startClust ,end.clus = endClust, reassign = TRUE, dist.fun=cosine_cluster )

}


