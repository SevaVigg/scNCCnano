#makeUMAPClusteringPlots 	<- function( seuratObj, umapDim, clResolution){

source("R/setClusterColors.r")
source("R/setCellTypeColors.r")
source("R/plot2DallLineages.r")
source("R/createSlingShotObject.r")

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypeDir <- file.path(plotDir, seuratObj@project.name)
dir.create( experimentTypeDir, showWarnings = FALSE)

pcaPlotDir 	<- file.path( experimentTypeDir, "PCAdimReduction")
dir.create( pcaPlotDir, showWarnings = FALSE)

PCAcomps	<- seuratObj@calc.params$RunPCA$pcs.compute

compsDir 	<- file.path( pcaPlotDir, paste0("comps", comps))
dir.create(compsDir, showWarnings = FALSE)

umapDir 	<- file.path( compsDir, "UMAPdimRed")
dir.create( umapDir, showWarnings = FALSE)

umapDim		<- 17

umapDimDir	<- file.path( umapDir, paste0("UMAP_", umapDim))
dir.create( umapDimDir, showWarnings = FALSE)

#Clustering with dimenstion reduction

clResolution	<- 0.2

resolDir	<- file.path( umapDimDir, paste0( "resol_", clResolution))
dir.create( resolDir, showWarnings = FALSE)

umapDimSeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( resolDir, "umapDimSeed.txt"), umapDimSeed, "\n")

seuratObj	<- RunUMAP( seuratObj, dims.use = 1 : PCAcomps, reduction.use = "pca", max.dim = umapDim, seed.use = umapDimSeed)


clustSeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( resolDir, "clustSeed.txt"), clustSeed, "\n")

seuratObj	<- BuildSNN( seuratObj, dims.use = 1 : umapDim, prune.SNN = 1/15, reduction.type = "umap")
seuratObj 	<- FindClusters( seuratObj, reuse.SNN = TRUE, resolution = clResolution, random.seed = clustSeed)
seuratObj	<- ValidateClusters( seuratObj, top.genes = 7, pc.use = 1 : PCAcomps, acc.cutoff = 0.9, min.connectivity = 0.05, verbose = TRUE)

clTypes <- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj 	<- BuildClusterTree( seuratObj, pcs.use = 1 : PCAcomps, do.plot = FALSE, do.reorder = TRUE) #This functions renames clusters, so we need to assign cluster types again

png( file.path( resolDir, "ClusterTreePCASpace.png"))
	PlotClusterTree( seuratObj)
dev.off()

clTypes <- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)

#TSNEPlot for clustered data
#png( file.path( resolDir, paste0("TSNEClustersPCASpace_c", comps, "_res", clResolution, ".png")))
#	TSNEPlot( seuratObj, colors.use = setClusterColors( seuratObj))
#dev.off()
seuratSling	<- createSlingShotObject( seuratObj, "umap")
 
seuratObj2D	<- calcUMAP_PCASpace( seuratObj, PCAcomps, umapDimSeed, 2) 

png( file.path( resolDir, "UMAPClustersPCASpace.png"))
	pl <- DimPlot(object = seuratObj2D, reduction.use = 'umap', pt.size = 2, cols.use = setClusterColors( seuratObj), do.return = TRUE, vector.friendly = TRUE)
	plot(pl) + labs( title = "Clusters")
dev.off()

png( file.path( resolDir, paste0("UMAP_", umapDim, "clRes", clResolution, ".png" )))
	plot2DallLineages( seuratSling, seuratObj2D, "umap")
dev.off()

seuratObj2D <- SetAllIdent( seuratObj2D, id = "originalCellTypes")
png( file.path( resolDir, "UmapInitCellTypesPCASpace.png"))
	pl <- DimPlot(object = seuratObj2D, reduction.use = "umap", pt.size = 2, cols.use = setCellTypeColors( seuratObj2D), do.return = TRUE, vector.friendly = TRUE)
	plot(pl) + labs( title = "InitialCellTypes")
dev.off()


#remove values, that are too close to zero

noiseTol	<- log2(19)
denoiseObj	<- seuratObj
denoiseObj@data	<- apply( denoiseObj@data, c(1,2), function(x) if(x>noiseTol) x else 0)

png( file.path( resolDir, paste0("DotPlotPCASpace_c", comps, "_res", clResolution,".png")), width = 800, height = 600)
	DotPlot(denoiseObj, genes.plot = rownames( denoiseObj@data), x.lab.rot = TRUE, dot.scale = 5, plot.legend = TRUE, dot.min = 0, scale.by = "radius")
dev.off()

return( seuratObj)
#}

