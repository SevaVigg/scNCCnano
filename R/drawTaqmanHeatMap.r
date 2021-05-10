drawTaqmanHeatMap	<- function( seuratObj,  clusterOrderedCells){

#This snippet makes a heatmap. 

require( ComplexHeatmap)

source("R/setClusterColors.r")

taqmanGenes <-  c( 	"sox10"		, "mitfa"	, "neurog1"	,"phox2b"	,
			"ltk"		, "pax7b"	, "pnp4a"	, "sox9b"	, 
			"snail2"	, "mbpa"	, "tyrp1b"	, "xdh"		, 
			"elavl3")


lineageCells	<- clusterOrderedCells
clusterColors	<- setClusterColors( seuratObj)
curveClust 	<- seuratObj@ident[ lineageCells]
targetExpsDF	<- seuratObj@data[ taqmanGenes, lineageCells ]
#listClDF	<- lapply( levels( seuratObj@ident) , function(x) seuratObj@data[ taqmanGenes, WhichCells( seuratObj, x)])

curveDF 	<- data.frame( clust = curveClust)

annotColors 		<- as.character(unique( clusterColors[ seuratObj@ident[ lineageCells]]))
names(annotColors) 	<- as.character(unique( seuratObj@ident[ lineageCells]))


annotBar		<- columnAnnotation( CellTypes = curveClust,
				col = list( CellTypes = annotColors), 
				height = unit(30, "points"),
				annotation_legend_param = list( ncol = 10)
			#annotation_legend_side = "bottom"
					)

hMap		<- Heatmap( 	targetExpsDF, 
			name = "logExps",
#			column_split = data.frame( ident = seuratObj@ident[ clusterOrderedCells]),
#			column_gap = unit(2, "mm"),
			cluster_columns = FALSE,
			cluster_rows = FALSE,
			show_column_names = FALSE, 
			row_names_gp = gpar(fontsize = 24), 
			column_names_gp = gpar(fontsize = 10),
#			bottom_annotation = annotBar,  
			column_title = paste0("sox10_mitfa_neurog1_phox2b_ltk_pax7b"), 
			clustering_distance_rows = "euclidean", 
			use_raster = TRUE, raster_device = "png", 
			)


return( hMap)

}
