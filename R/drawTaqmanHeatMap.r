drawTaqmanHeatMap	<- function( seuratObj, clusterOrderedCells, heatMapHeight, heatMapWidth, showCellNames = FALSE){

#This snippet makes a heatmap. 

require( ComplexHeatmap)
require( circlize)
require( viridis)
require( proxy)
library( dendsort)


source("R/setClusterColors.r")

taqmanGenes <-  c( 	"sox10"		, "mitfa"	, "neurog1"	,"phox2b"	,
			"ltk"		, "pax7b"	, "pnp4a"	, "sox9b"	, 
			"snail2"	, "mbpa"	, "tyrp1b"	, "xdh"		, 
			"elavl3")

dataMatrix 	<- seuratObj@data[taqmanGenes, clusterOrderedCells] 

nCol		<- 1024

row_dend 	<- dendsort(hclust( dist(dataMatrix, method = "cosine", pairwise = TRUE)))
column_dend	<- dendsort(hclust( dist( t(dataMatrix), method = "cosine", pairwise = TRUE)))

#listClDF	<- lapply( levels( seuratObj@ident) , function(x) seuratObj@data[ taqmanGenes, WhichCells( seuratObj, x)])


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
			row_names_gp = gpar(fontsize = 16), 
			column_names_gp = gpar(fontsize = 16),
			clustering_distance_rows = function(x) dist( x, method = "cosine", pairwise = TRUE), 
			clustering_distance_columns = function(x) dist( x, method = "cosine", pairwise = TRUE), 
			use_raster = TRUE, raster_device = "png"
			)



return( hMap)

}
