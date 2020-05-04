drawHeatMap	<- function( seuratObj,  targetCurve, do.print = FALSE){

#This snippet makes a heatmap. 

source("R/setClusterColors.r")

tanyaGenes <- c("sox10","sox9b","snail2","foxd3","phox2b",
		"tfap2a", "tfap2e",
		"impdh1b","kita","mbpa",
		"mitfa","ltk","tfec",
		"pax7a", "pax7b","pnp4a",
		"myo5aa", "tyrp1b","mlphb","oca2","silva","slc24a5")

lineageCells	<- names(targetCurve$cells)
lineageType	<- targetCurve$target
clusterColors	<- setClusterColors( seuratObj)
curveClust 	<- seuratObj@ident[ lineageCells]
targetExpsDF	<- seuratObj@data[ tanyaGenes, lineageCells ]
clusterColors	<- setClusterColors( seuratObj)

curveDF 	<- data.frame( clust = curveClust)

annotColors 		<- as.character(unique( clusterColors[ seuratObj@ident[ lineageCells]]))
names(annotColors) 	<- as.character(unique( seuratObj@ident[ lineageCells]))


topbar		<- columnAnnotation( curveDF,
				    col = list(clust = annotColors), 
				    height = unit(30, "points"),
				    annotation_legend_param = list( nrow = 1)
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

png( file.path( heatMapDir, paste0( "HeatMap_", lineageType, "_Lineage.png" )))
	draw(hplot, annotation_legend_side = "bottom")
dev.off()

}else{ draw(hplot, annotation_legend_side = "bottom")}

}
