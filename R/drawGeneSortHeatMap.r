drawGeneSortHeatMap	<- function( seuratObj,  heatMapHeight, heatMapWidth, showCellNames = FALSE){

#This snippet makes a heatmap of gene covariation. 

require( ComplexHeatmap)
require( circlize)
require( viridis)
require( proxy)
library( dendsort)


heatMapGenes	<- rownames( seuratObj@data) 

nCol		<- 1024


geneOrder	<- 	c( 	"ltk"		,"sox9b"	, "foxo1a"	, "tfap2e"	, "tfap2a"	, "her9"	, 
				"foxg1b" 	, "snai1b"	, "alx4b"	, "hmx1"	, 
				"otx2b"		, "sox10"	, "impdh1b"	, "foxo1b"	, "tyr"		, 
				"pax7b"		, "mc1r"	, "id2a"	, "hmx4"	, "foxd3"	, 
				"ednrba"	, "kita"	, "mbpa"	, "phox2bb"	, "tfec"	, 
				"mitfa"		, "foxp4" 	, "hbp1"	,  "mycla"	, "pax7a"	, 
				"tyrp1b"	, "slc24a5"	, "oca2"	, "mlphb"	, "pmela"	, 
				"myo5aa"	, "pnp4a"	, "ets1"	, "fgfr3"	, "pax3a"	, 
				"dpf3"		, "smad9")


dataMatrix 	<- seuratObj@data[ geneOrder , order( seuratObj@data[ "ltk", ],  decreasing = TRUE)]

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
			cluster_columns = FALSE,
			cluster_rows = FALSE,
			row_dend_reorder = FALSE,
			show_column_names = showCellNames, 
			row_names_gp = gpar(fontsize = 12), 
			column_names_gp = gpar(fontsize = 12),
#			clustering_distance_rows = function(x) dist( x, method = "cosine", pairwise = TRUE), 
			use_raster = TRUE, raster_device = "png"
			)


return( hMap)

}
