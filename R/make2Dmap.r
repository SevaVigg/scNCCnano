make2Dmap <- function( seurObj){

#UMAPSeed <- as.numeric(as.POSIXct(Sys.time()))
UMAPSeed <- 42

seurObj <- calcUmapGeneSpace( seurObj, experimentType = "allCells", Dim = 2, minDist = 4, mySpread = 4, UMAPRandSeed = UMAPSeed)$All
save( seurObj, file = file.path( clusterDataDir, "visualisation2Dumap"))


seurObj <- SetAllIdent( seurObj, id = "genCellTypeIdent")
png( file.path( clusterDataDir, "control2DinitCellPlot.png"))
	DimPlot( seurObj, reduction.use = "umap", cols.use = setClusterColors( seurObj))
dev.off()
png( file.path( clusterDataDir, "controlLtkFeature.png"))
	FeaturePlot( seurObj, reduction.use = "umap", features.plot = "ltk")
dev.off()

return()
}

