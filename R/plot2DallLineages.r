plot2DallLineages <- function( slingShotObj, seuratObj){

#seuratObject must contain results of clustering and tSNE pre-calculated

source("R/setClusterColors.r")
source("R/getClusterTypes.r")
source("R/getLineageCoords.r")

clTypes		<- getClusterTypes( seuratObj)
plotVals	<- seuratObj@dr$tsne@cell.embeddings

source("R/getLineageCoords.r")
LineageTree	<- getLineageCoords( seuratObj, slingShotObj)  
 
cellColors <- setClusterColors(clTypes)[ seuratObj@ident]

plot( plotVals[,2]~plotVals[,1], cex = 0.7, pch = 16, col = cellColors, xlab = "tSNE1", ylab = "tSNE2")
sapply(LineageTree, function(x) {lines(t(x)); points(t(x), pch = 16)})
}


