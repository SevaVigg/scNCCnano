drawTargetHeatMap	<- function( seuratObj,  targetCurve){

#This snippet makes a heatmap. 

require( ComplexHeatmap)

source("R/setClusterColors.r")

tanyaGenes 	<- c( 	"tfap2e"	, "tfap2a"	, "her9"	, "sox9b"	, "foxg1b" 	,
			"snai1b"	, "alx4b"	, "hmx1"	, "otx2b"	, "sox10"	, 
			"impdh1b"	, "foxo1b"	, "tyr"		, "pax7b"	, "mc1r"	, 
			"id2a"		, "hmx4"	, "foxd3"	, "ednrba"	, "kita"	, 
			"mbpa"		, "phox2bb"	, "tfec"	, "mitfa"	, "foxp4" 	, 
			"foxo1a"	, "hbp1"	, "ltk"		, "mycla"	, "pax7a"	, 
			"tyrp1b"	, "slc24a5"	, "oca2"	, "mlphb"	, "pmela"	, 
			"myo5aa"	, "pnp4a"	, "ets1"	, "fgfr3"	, "pax3a"	, 
			"dpf3"		, "smad9")

 
lineageCells	<- names(targetCurve$cells)
lineageType	<- targetCurve$target
clusterColors	<- setClusterColors( seuratObj)
names(clusterColors)	<- levels( seuratObj@ident)
curveClust 	<- seuratObj@ident[ lineageCells]
targetExpsDF	<- seuratObj@data[ tanyaGenes, lineageCells ]

curveDF 	<- data.frame( clust = curveClust)

#annotColors 		<- as.character(unique( clusterColors[ seuratObj@ident[ lineageCells]]))
#names(annotColors) 	<- as.character(unique( seuratObj@ident[ lineageCells]))


annotbar		<- HeatmapAnnotation( Cluster = seuratObj@ident[lineageCells],
				    col = list( Cluster = clusterColors), 
				    height = unit(30, "points")
				    #annotation_legend_param = list( nrow = 10)
					)

hplot <- Heatmap( 	targetExpsDF, 
			name = "log2Exp", 
			cluster_columns = FALSE,
			cluster_rows = FALSE, 
			show_column_names = FALSE, 
			row_names_gp = gpar(fontsize = 24), 
			column_names_gp = gpar(fontsize = 10), 
			bottom_annotation = annotbar,
			column_title = paste0("Expression of ", lineageType, " lineage"), 
			clustering_distance_rows = "euclidean", 
			use_raster = TRUE, raster_device = "png", 
			)

if (do.print){return( hplot)}else{ draw(hplot, annotation_legend_side = "bottom")}

}
