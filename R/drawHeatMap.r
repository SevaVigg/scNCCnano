drawHeatMap	<- function( seuratObj,  targetCurve, do.print = FALSE){

#This snippet makes a heatmap. 

require( ComplexHeatmap)

source("R/setClusterColors.r")

tanyaGenes <- c("sox10","sox9b","snai1b","foxd3","phox2bb",
		"tfap2a", "tfap2e",
		"impdh1b","kita","mbpa",
		"mitfa","ltk","tfec", "foxo1a", 
		"pax7a", "pax7b","pnp4a",
		"myo5aa", "tyrp1b","mlphb","oca2","pmela","slc24a5")

lineageCells	<- names(targetCurve$cells)
lineageType	<- targetCurve$target
clusterColors	<- setClusterColors( seuratObj)
curveClust 	<- seuratObj@ident[ lineageCells]
targetExpsDF	<- seuratObj@data[ tanyaGenes, lineageCells ]

curveDF 	<- data.frame( clust = curveClust)

annotColors 		<- as.character(unique( clusterColors[ seuratObj@ident[ lineageCells]]))
names(annotColors) 	<- as.character(unique( seuratObj@ident[ lineageCells]))


topbar		<- columnAnnotation( curveDF,
				    col = list(clust = annotColors), 
				    height = unit(30, "points"),
				    annotation_legend_param = list( nrow = 10)
					)

hplot <- Heatmap( 	targetExpsDF, 
			name = "log2Exp", 
			cluster_columns = FALSE,
			cluster_rows = FALSE, 
			show_column_names = FALSE, 
			row_names_gp = gpar(fontsize = 24), 
			column_names_gp = gpar(fontsize = 10), 
			bottom_annotation = topbar,
			column_title = paste0("Expression of ", lineageType, " lineage"), 
			clustering_distance_rows = "euclidean", 
			use_raster = TRUE, raster_device = "png", 
			)

if (do.print){return( hplot)}else{ draw(hplot, annotation_legend_side = "bottom")}

}
