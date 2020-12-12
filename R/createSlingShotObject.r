createSlingShotObject <- function( seuratObj, dimRed, startClust = "eNCC", endClust = c("I", "M", "X")){

source("R/getFinalClusterTypes.r")

if(!require(slingshot)){
  install.packages("slingshot")
  library("slingshot")	
}

coordMatrixMD  <-  as.matrix(seuratObj@dr[[dimRed]]@cell.embeddings)

clTypes		<- getFinalClusterTypes( seuratObj)
slingShotObj 	<- slingshot(coordMatrixMD, seuratObj@ident, start.clus = startClust ,end.clus = endClust, reassign = TRUE)

}


