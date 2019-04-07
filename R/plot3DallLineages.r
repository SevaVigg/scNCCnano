plot3DallLineages <- function( slingShotObj, seuratObj, dimRed){

#seuratObject must contain results of clustering and tSNE pre-calculated

source("R/setClusterColors.r")
source("R/getLineageCoords.r")
source("R/dimPlot3DRGL.r")
source("R/getLineageCoords.r")

LineageTree	<- getLineageCoords( seuratObj, slingShotObj, dimRed)  
 
dimPlot3DRGL( seuratObj, "umap")
sapply(LineageTree, function(x) {rgl.linestrips(t(x)); rgl.points(t(x), size = 4)})
}


