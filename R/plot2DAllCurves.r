#plot2DAllCurves	<- function( seuratObj, slingShotObj, lineageTypes = c("M", "I"), dimRed){

#seuratObject must contain results of clustering and tSNE and UMAP 2D pre-calculated

source("R/setClusterColors.r")
source("R/getLineageCoords.r")

lineageIds 	<- which(unlist(lapply( slingShotObj@lineages, function(x) tail(x, 1) %in% lineageTypes))) 

prinCurveDF	<- slingPseudotime( slingShotObj)
curveWeightDF	<- slingCurveWeights( slingShotObj)
prinCurve_F	<- prinCurveDF[ which(curveWeightDF[ ,LineageId] > 0.995) , LineageId]
lineageLines	<- as.data.frame(seuratObj@dr[[ dimRed ]]@cell.embeddings[ names( prinCurve_F), ])

#make 2D umap for visulaization
seuratObj	<- RunUMAP( seuratObj, genes.use = rownames( seuratObj@data), max.dim = 2, reduction.name = 'umap', n_neighbors = 15L, 
				min_dist = 0.39, metric = "cosine", seed.use = visSeed, spread = 1 )

plotVals	<- seuratObj@dr[[ dimRed ]]@cell.embeddings
cellColors 	<- setClusterColors( seuratObj)[ seuratObj@ident]

lineageLines	<- as.data.frame(seuratObj@dr[[ dimRed ]]@cell.embeddings[ names( sort(prinCurve_F)), ])

plot( plotVals, col = cellColors, cex = 1, pch = 16)

for (LineageId in LineageIds){


sapply( lineageLines, function(x) {lines(t(x)); points(t(x), pch = 16)})
}

#}
