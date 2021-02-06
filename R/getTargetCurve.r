getTargetCurve <- function( seuratObj, lineageStart = "eNCC", lineageEnds = c("M", "I", "X"), target = "I", dimRed = "umap", distFun = cosineClusterDist){

require( "slingshot")
source("R/setCellTypeColors.r")
source("R/cosineClusterDist.r")

cells 		<- GetCellEmbeddings( seuratObj, reduction.type = dimRed)

slingLins	<- getLineages( cells, seuratObj@ident, start.clus = lineageStart, end.clus = lineageEnds, dist.fun = distFun) 
slingCurvs	<- getCurves( slingLins, extend = "n", stretch = 0, thresh = 0.05)

LineageId <- which(unlist(lapply( slingLins@lineages, function(x) tail(x, 1) == target)))  

#slingShotObj	<- getCurves(slingShotObj)
prinCurveDF	<- slingPseudotime( slingCurvs)
curveWeightDF	<- slingCurveWeights( slingCurvs)
prinCurve_F	<- prinCurveDF[ which(curveWeightDF[ ,LineageId] > 0.) , LineageId]
curveCells 	<- sort( prinCurve_F)
targetCurve 	<- list( target = target, cells = curveCells)

return( targetCurve)
}
 
