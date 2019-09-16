#Clustering in the initial gene space reduction. 
#Needs makeGeneSpacePlots.r to be run first
#returns seuart Object with clustering

makeGeneSpaceClusteringPlots <- function( seuratObj, clResolutin)

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypeDir <- file.path(plotDir, seuratObj@project.name)
dir.create( experimentTypeDir, showWarnings = FALSE)

geneSpacePlotDir <- file.path( experimentTypeDir, "geneSpace")
dir.create( geneSpacePlotDir, showWarnings = FALSE)

resolDir	<- file.path( geneSpacePlotDir, paste0( "resol_", clResolution))
dir.create( resolDir, showWarnings = FALSE)

clustSeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( resolDir, "clustSeed.txt"), clustSeed, "\n")

seuratObj	<- BuildSNN( seuratObj, genes.use = rownames(seuratObj@data), k.param = 20, prune.SNN = 0)
seuratObj 	<- FindClusters( seuratObj, reuse.SNN = TRUE, resolution = clResolution, random.seed = clustSeed)
seuratObj	<- ValidateClusters( seuratObj, top.genes = 7, pc.use = NULL, acc.cutoff = 0.005, min.connectivity = 0.05, verbose = TRUE)

clTypes <- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj 	<- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = TRUE) #This functions renames clusters, so we need to assign cluster types again

png( file.path( resolDir, "ClusterTreeGeneSpace.png"))
	PlotClusterTree( seuratObj)
dev.off()

clTypes <- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)

#End construct TSNEPlot for clustered data
png( file.path( geneSpacePlotDir, "TSNEClustersGeneSpace.png"))
	TSNEPlot( seuratObj, colors.use = setClusterColors( seuratObj))
dev.off()

png( file.path( geneSpacePlotDir, "UMAPClustersGeneSpace.png"))
	DimPlot(object = seuratObj, reduction.use = 'umap', pt.size = 1, cols.use = setClusterColors( seuratObj))
dev.off()


#make DotPlot

#remove values close to zero (to calculate radii)
noiseTol	<- log2(19)

seuratObjDenoised	<- seuratObj
seuratObjDenoised@data <- apply( seuratObjDenoised@data, c(1,2), function(x) if(x>noiseTol) x else 0)

png( file.path( geneSpacePlotDir, "DotPlotClustersGeneSpace.png"), width = 800, height = 600)
	DotPlot(seuratObjDenoised, genes.plot = rownames(seuratObjDenoised@data), x.lab.rot = TRUE, dot.scale = 5, plot.legend = TRUE, dot.min = 0, scale.by = "radius")
dev.off()

return( seuratObj)
}

