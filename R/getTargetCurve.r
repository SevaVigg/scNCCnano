getTargetCurve <- function( seuratObj, lineageStart = "eHMP", genes.use = NULL, lineageEnds = c("M", "I", "X"), target = "I", dimRed = "umap", dimsUseHD = 1:2, distFun = cosineClusterDist){

require( "slingshot")
source("R/setCellTypeColors.r")
source("R/cosineClusterDist.r")

if ( !is.null( genes.use)) {
	cells 	<- t( as.matrix(seuratObj@data[ genes.use, ]))

}else{ cells 	<- GetCellEmbeddings( seuratObj, reduction.type = dimRed)[ , dimsUseHD]}

slingLins	<- getLineages( cells, seuratObj@ident, start.clus = lineageStart, end.clus = lineageEnds, dist.fun = distFun) 
slingCurvs 	<- getCurves( slingLins, extend = "n", reassign = TRUE, stretch = 0, thresh = 0.0001, shrink = 0.2)


LineageId <- which(unlist(lapply( slingLins@lineages, function(x) tail(x, 1) == target)))  

#slingShotObj	<- getCurves(slingShotObj)
prinCurveDF	<- slingPseudotime( slingCurvs)
curveWeightDF	<- slingCurveWeights( slingCurvs)
prinCurve_F	<- prinCurveDF[ which(curveWeightDF[ ,LineageId] > 0.) , LineageId]
curveCells 	<- sort( prinCurve_F)
targetCurve 	<- list( target = target, cells = curveCells)

return( targetCurve)
}
 
