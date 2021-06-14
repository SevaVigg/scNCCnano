drawBiclustHeatMap	<- function( seuratObj , heatMapHeight, heatMapWidth, showCellNames = FALSE){

# This snippet makes a heatmap with gene expression biclustering 
# written by Vsevolod J. Makeev 2017 - 2021
# inspired by GitHub code of hrbrmstr in hrbrmstr/viridis_chords.R



require( ComplexHeatmap)
require( circlize)
require( viridis)
require( proxy)
library( dendsort)


heatMapGenes	<- rownames( seuratObj@data) 
dataMatrix <- seuratObj@data

nCol		<- 1024

row_dend 	<- dendsort(hclust( dist(dataMatrix, method = "cosine", pairwise = TRUE)))
column_dend	<- dendsort(hclust( dist( t(dataMatrix), method = "cosine", pairwise = TRUE)))

colPanelFun	 = colorRamp2( quantile( dataMatrix, seq(0, 1, by = 1/(nCol - 1))), viridis( nCol))

hMap		<- Heatmap( 	dataMatrix, 
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
			cluster_columns = column_dend,
			cluster_rows = row_dend,
			row_dend_reorder = TRUE,
			column_dend_reorder = TRUE,
			show_column_names = showCellNames, 
			row_names_gp = gpar(fontsize = 12), 
			column_names_gp = gpar(fontsize = 12),
			clustering_distance_rows = function(x) dist( x, method = "cosine", pairwise = TRUE), 
			clustering_distance_columns = function(x) dist( x, method = "cosine", pairwise = TRUE), 
			use_raster = TRUE, raster_device = "png"
			)


return( hMap)

}
