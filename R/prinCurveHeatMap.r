getCurveHeatMap <- function( seuratObj, slingShotObj, LineageId, dimRed){

source("R/setCellTypeColors.r")

if(!require(ComplexHeatmap)){
	source("https://bioconductor.org/biocLite.R")
	biocLite("ComplexHeatmap")
}


tanyaGenes <- c("foxd3", "impdh1b",  "kita", "ltk", "mbpa", "mitfa", "mlphb", "oca2", "phox2b", "pnp4a", "silva", "slc24a5", "snail2",
                 "sox10", "sox9b", "tfec","tyrp1b",  "pax7b", "tfap2e", "tfap2a")

source("R/getLineageCoords.r")
LineageTree	<- getLineageCoords( seuratObj, slingShotObj, dimRed) 
slingShotObj	<- getCurves(slingShotObj)
prinCurveDF	<- slingPseudotime( slingShotObj)
targetCurve	<- sort(prinCurveDF[ !is.na( prinCurveDF[ , LineageId]), LineageId])
cellColors	<- setClusterColors( seuratObj)[ seuratObj@ident]
names(cellColors) <- names(seuratObj@ident)
LineageType	<- tail(slingShotObj@lineages[LineageId][[1]],1)

heatmapPlotDir 	<- file.path(resolDir, "HeatMaps")
dir.create(heatmapPlotDir, showWarnings = FALSE)

curveClust 	<- seuratObj@ident[names( targetCurve)]
curveDF 	<- data.frame( clust = curveClust)

plot(0,0, main = paste0("Expression of ", LineageType, " lineage"))

topbar		<- columnAnnotation( curveDF, 
				    col = list( clust = cellColors(rownames(curveDF))), 
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
			column_title = "Lineage_expression", 
			clustering_distance_rows = "euclidean", 
			use_raster = TRUE, raster_device = "png", 
			)
draw(hplot, annotation_legend_side = "bottom")
} 
