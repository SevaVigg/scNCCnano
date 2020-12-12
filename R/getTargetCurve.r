getTargetCurve <- function( seuratObj, lineageStart = "eNCC", lineageEnds = c("M", "I", "X"), target = "I", dimRed = "umap"){

require( "slingshot")
source("R/setCellTypeColors.r")

if(!require(ComplexHeatmap)){
	source("https://bioconductor.org/biocLite.R")
	biocLite("ComplexHeatmap")
}

cells 		<- GetCellEmbeddings( seuratObj, reduction.type = dimRed)
slingLins	<- getLineages( cells, seuratObj@ident, start.clus = lineageStart, end.clus = lineageEnds) 
slingLins	<- getCurves( slingLins, extend = "n")

LineageId <- which(unlist(lapply( slingLins@lineages, function(x) tail(x, 1) == target)))  

#slingShotObj	<- getCurves(slingShotObj)
prinCurveDF	<- slingPseudotime( slingLins)
curveWeightDF	<- slingCurveWeights( slingLins)
prinCurve_F	<- prinCurveDF[ which(curveWeightDF[ ,LineageId] > 0.99) , LineageId]
curveCells 	<- sort( prinCurve_F)
targetCurve 	<- list( target = target, cells = curveCells)

return( targetCurve)
}
 
