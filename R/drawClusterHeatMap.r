drawClusterHeatMap	<- function( seuratObj , heatMapHeight, heatMapWidth){

#This snippet makes a heatmap of gene covariation. 

require( ComplexHeatmap)
require( circlize)
require( viridis)
require( proxy)
require( dendsort)

source( "R/setClusterColors.r")


heatMapGenes	<- rownames( seuratObj@data) 
dataMatrix 	<- seuratObj@data

nCol		<- 1024

clusterColors		<- setClusterColors( seuratObj)
names(clusterColors)	<- levels( seuratObj@ident)

orderedCellTypes	<- ordered( seuratObj@ident, levels = c( "eHMP", "ltHMP", "I", "M", "X", "7", "4"))
orderedCellTypes	<- orderedCellTypes[ order( orderedCellTypes)]

#annotColors 		<- as.character(unique( clusterColors[ seuratObj@ident[ orderedCellTypes]]))
#names(annotColors) 	<- as.character(unique( seuratObj@ident[ orderedCellTypes]))


row_dend 	<- dendsort(hclust( dist(  dataMatrix, method = "cosine", pairwise = TRUE)))
colPanelFun	 = colorRamp2( quantile( dataMatrix, seq(0, 1, by = 1/(nCol - 1))), viridis( nCol))

annotBar		<- HeatmapAnnotation( Cluster = orderedCellTypes,
				col = list( Cluster = clusterColors) 
				#height = unit(30, "points"),
				#annotation_legend_param = list( ncol = 10)
				#annotation_legend_side = "bottom"
					)
hMap		<- Heatmap( 	dataMatrix[ , names(orderedCellTypes)], 
			name 	= "gene expression",
			width 	= unit( heatMapWidth,  "inches"), 
			height 	= unit( heatMapHeight, "inches"),
			heatmap_legend_param = list(
				title = "gene exp",
				legend_height = unit( 2, "inches"),
				grid_width = unit( 0.5, "inches"),
				title_gp = gpar( fontsize = 16, fontface = "bold")		
						),
			col	=  colPanelFun,
			bottom_annotation = annotBar,  
			column_order = names(orderedCellTypes),
			column_split = orderedCellTypes,
#			cluster_column_slices = FALSE, 
			cluster_rows = row_dend,
			cluster_columns = FALSE,
			row_dend_reorder = FALSE,
			column_dend_reorder = FALSE,
			show_column_names = FALSE, 
			row_names_gp = gpar(fontsize = 12),
#			column_names_gp = gpar(fontsize = 16),
			clustering_distance_rows = function(x) dist( x, method = "cosine", pairwise = TRUE), 
#			clustering_distance_columns = function(x) dist( x, method = "cosine", pairwise = TRUE), 
			use_raster = TRUE, raster_device = "png"
			)


return( hMap)

}
