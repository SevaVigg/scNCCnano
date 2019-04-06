createSlingShotObject <- function( coordMatrixMD, seuratObj){

source("R/getClusterTypes.r")

if(!require(slingshot)){
  install.packages("slingshot")
  library("slingshot")	
}

clTypes		<- getClusterTypes( seuratObj)
slingShotObj 	<- slingshot(coordMatrixMD, seuratObj@ident, start.clus = clTypes["E"],end.clus=c(clTypes["I"], clTypes["M"]))

}


