getCurveHeatMap <- function( seuratObj, slingShotObj, LineageType, dimRed, do.print = FALSE){

source("R/setCellTypeColors.r")

if(!require(ComplexHeatmap)){
	source("https://bioconductor.org/biocLite.R")
	biocLite("ComplexHeatmap")
}


tanyaGenes <- c("sox10","sox9b","snail2","foxd3","phox2b",
		"tfap2a", "tfap2e",
		"impdh1b","kita","mbpa",
		"mitfa","ltk","tfec",
		"pax7a", "pax7b","pnp4a",
		"tyrp1b","mlphb","oca2","silva","slc24a5")

source("R/getLineageCoords.r")

LineageId <- which(unlist(lapply( slingShotObj@lineages, function(x) tail(x, 1) == LineageType)))  

#slingShotObj	<- getCurves(slingShotObj)
prinCurveDF	<- slingPseudotime( slingShotObj)
curveWeightDF	<- slingCurveWeights( slingShotObj)
prinCurve_F	<- prinCurveDF[ which(curveWeightDF[ ,LineageId] > 0.995) , LineageId]
targetCurve	<- logExps[ tanyaGenes, names( sort( prinCurve_F))]
clusterColors	<- setClusterColors( seuratObj)
LineageType	<- tail(slingShotObj@lineages[LineageId][[1]],1)

curveClust 	<- seuratObj@ident[ colnames( targetCurve)]
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
			cluster_rows = FALSE, 
			show_column_names = FALSE, 
			row_names_gp = gpar(fontsize = 24), 
			column_names_gp = gpar(fontsize = 10), 
			bottom_annotation = topbar,
			column_title = paste0("Expression of ", LineageType, " lineage"), 
			clustering_distance_rows = "euclidean", 
			use_raster = TRUE, raster_device = "png", 
			)

if (do.print){

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypePlotDir <- file.path(plotDir, seuratObj@project.name)
dir.create( experimentTypePlotDir, showWarnings = FALSE)

geneSpacePlotDir <- file.path( experimentTypePlotDir, "geneSpacePlots")
dir.create( geneSpacePlotDir, showWarnings = FALSE)

heatMapDir <- file.path( geneSpacePlotDir, "heatMap")
dir.create( heatMapDir, showWarnings = FALSE)

png( file.path( heatMapDir, paste0( "HeatMap_", LineageType, "_Lineage.png" )))
	draw(hplot, annotation_legend_side = "bottom")
dev.off()

}else{ draw(hplot, annotation_legend_side = "bottom")}
} 
