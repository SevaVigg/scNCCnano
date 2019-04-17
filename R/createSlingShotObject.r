createSlingShotObject <- function( seuratObj, dimRed){

source("R/getClusterTypes.r")

if(!require(slingshot)){
  install.packages("slingshot")
  library("slingshot")	
}

coordMatrixMD  <-  seuratObj@dr[[dimRed]]@cell.embeddings

clTypes		<- getClusterTypes( seuratObj)
slingShotObj 	<- slingshot(coordMatrixMD, seuratObj@ident, start.clus = clTypes["E"],end.clus=c(clTypes["I"], clTypes["M"]))

}


