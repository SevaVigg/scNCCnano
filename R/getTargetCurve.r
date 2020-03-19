getTargetCurve <- function( slingShotObj, LineageType){

source("R/setCellTypeColors.r")

if(!require(ComplexHeatmap)){
	source("https://bioconductor.org/biocLite.R")
	biocLite("ComplexHeatmap")
}
source("R/getLineageCoords.r")

LineageId <- which(unlist(lapply( slingShotObj@lineages, function(x) tail(x, 1) == LineageType)))  

#slingShotObj	<- getCurves(slingShotObj)
prinCurveDF	<- slingPseudotime( slingShotObj)
curveWeightDF	<- slingCurveWeights( slingShotObj)
prinCurve_F	<- prinCurveDF[ which(curveWeightDF[ ,LineageId] > 0.995) , LineageId]
curveCells 	<- sort( prinCurve_F)
targetCurve 	<- list( target = LineageType, cells = curveCells)

return( targetCurve)
}
 
