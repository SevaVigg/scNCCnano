plot3DallLineagesFeature <- function( slingShotObj, seuratObj, gene, dimRed){

#seuratObject must contain results of clustering and tSNE pre-calculated

source("R/setClusterColors.r")
source("R/getLineageCoords.r")
source("R/dimPlot3DFeature.r")
source("R/getLineageCoords.r")

LineageTree	<- getLineageCoords( seuratObj, slingShotObj, dimRed)  
 
dimPlot3DFeature( seuratObj, gene, "umap")
sapply(LineageTree, function(x) {rgl.linestrips(t(x)); rgl.points(t(x), size = 4)})
}


