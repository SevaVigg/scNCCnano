drawGeneCovarHeatMap	<- function( seuratObj , heatMapHeight, heatMapWidth){

# This snippet makes a heatmap of gene covariation matrix.
# written by Vsevolod J. Makeev 2017 - 2021
# inspired by GitHub code of hrbrmstr in hrbrmstr/viridis_chords.R

require( ComplexHeatmap)
require( circlize)
require( viridis)
require( proxy)

heatMapGenes	<- rownames( seuratObj@data) 
covar <- cor(t(seuratObj@data))

nCol		<- 1024

colPanelFun	 = colorRamp2( quantile( covar, seq(0, 1, by = 1/(nCol - 1))), viridis( nCol))

hMap		<- Heatmap( 	covar, 
			name 	= "gene covariance",
			width 	= unit( heatMapWidth,  "inches"), 
			height 	= unit( heatMapHeight, "inches"),
			heatmap_legend_param = list(
				title = "gene covar",
				legend_height = unit( 2, "inches"),
				grid_width = unit( 0.5, "inches"),
				title_gp = gpar( fontsize = 16, fontface = "bold")		
						),
			col	=  colPanelFun,
			cluster_columns = TRUE,
			cluster_rows = TRUE,
			show_column_names = TRUE, 
			row_names_gp = gpar(fontsize = 16), 
			column_names_gp = gpar(fontsize = 16),
			clustering_distance_rows = function(x) dist( x, method = "cosine", pairwise = TRUE), 
			clustering_distance_columns = function(x) dist( x, method = "cosine", pairwise = TRUE), 
			use_raster = TRUE, raster_device = "png"
			)


return( hMap)

}
