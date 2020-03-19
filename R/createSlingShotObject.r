createSlingShotObject <- function( seuratObj, dimRed, maxDim, startClust = "E"){

source("R/getClusterTypes.r")

if(!require(slingshot)){
  install.packages("slingshot")
  library("slingshot")	
}

coordMatrixMD  <-  seuratObj@dr[[dimRed]]@cell.embeddings[,1:maxDim ]

clTypes		<- getClusterTypes( seuratObj)
slingShotObj 	<- slingshot(coordMatrixMD, seuratObj@ident, start.clus = startClust ,end.clus=c(clTypes["I"], clTypes["M"]), reassign = TRUE)

}


