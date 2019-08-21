getCurveHeatMap <- function( seuratObj, slingShotObj, LineageId, dimRed){

source("R/setCellTypeColors.r")

if(!require(ComplexHeatmap)){
	source("https://bioconductor.org/biocLite.R")
	biocLite("ComplexHeatmap")
}


tanyaGenes <- c("foxd3", "impdh1b",  "kita",     "ltk",      "mbpa", "mitfa", "mlphb", "oca2", "phox2b", "pnp4a", "silva", "slc24a5", "snail2",
                 "sox10", "sox9b", "tfec","tyrp1b",  "pax7b", "tfap2e", "tfap2a")

source("R/getLineageCoords.r")

#slingShotObj	<- getCurves(slingShotObj)
prinCurveDF	<- slingPseudotime( slingShotObj)
curveWeightDF	<- slingCurveWeights( slingShotObj)
prinCurve_F	<- prinCurveDF[ which(curveWeightDF[ ,LineageId] > 0.995) , LineageId]
targetCurve	<- logExps[ tanyaGenes, names( sort( prinCurve_F))]
clusterColors	<- setClusterColors( seuratObj)
LineageType	<- tail(slingShotObj@lineages[LineageId][[1]],1)

curveClust 	<- seuratObj@ident[names( targetCurve)]
curveDF 	<- data.frame( clust = curveClust)

annotColors <- as.character(unique(clusterColors[curveDF$clust]))
names(annotColors) <- as.character(unique(curveDF$clust))


topbar		<- columnAnnotation( curveDF,
				    col = list(clust = annotColors), 
				    height = unit(30, "points"),
				    annotation_legend_param = list( nrow = 1)
					)

hplot <- Heatmap( 	targetCurve, 
			name = "log2Exp", 
			cluster_columns = FALSE, 
			show_column_names = FALSE, 
			row_names_gp = gpar(fontsize = 24), 
			column_names_gp = gpar(fontsize = 10), 
			bottom_annotation = topbar,
			column_title = paste0("Expression of ", LineageType, " lineage"), 
			clustering_distance_rows = "euclidean", 
			use_raster = TRUE, raster_device = "png", 
			)
draw(hplot, annotation_legend_side = "bottom")
} 
