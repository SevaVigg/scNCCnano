plot2DallLineages <- function( slingShotObj, seuratObj, dimRed){

#seuratObject must contain results of clustering and tSNE pre-calculated

source("R/setClusterColors.r")
source("R/getLineageCoords.r")

plotVals	<- seuratObj@dr$[[ dimRed ]]@cell.embeddings

source("R/getLineageCoords.r")
LineageTree	<- getLineageCoords( seuratObj, slingShotObj, dimRed)  
 
cellColors <- setClusterColors( seuratObj)[ seuratObj@ident]

plot( plotVals[,2]~plotVals[,1], cex = 0.7, pch = 16, col = cellColors, xlab = paste0(dimRed, "_1"), ylab = paste0(dimRed, "_2")
sapply(LineageTree, function(x) {lines(t(x)); points(t(x), pch = 16)})
}


