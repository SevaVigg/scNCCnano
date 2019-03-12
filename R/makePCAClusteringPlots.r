makePCAClusteringPlots 	<- function( ipmc, clResolution){

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypeDir <- file.path(plotDir, ipmc@project.name)
dir.create( experimentTypeDir, showWarnings = FALSE)

pcaPlotDir 	<- file.path( experimentTypeDir, "PCAdimReduction")
dir.create( pcaPlotDir, showWarnings = FALSE)

comps		<- ipmc@calc.params$RunPCA$pcs.compute

compsDir 	<- file.path( pcaPlotDir, paste0("comps", comps))
dir.create(compsDir, showWarnings = FALSE)


#Clustering with dimenstion reduction
#clResolution	<- 1.2
resolDir	<- file.path( compsDir, paste0( "resol_", clResolution))
dir.create( resolDir, showWarnings = FALSE)

clustSeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( resolDir, "clustSeed.txt"), clustSeed, "\n")

ipmc	<- BuildSNN( ipmc, dims.use = 1:comps, prune.SNN = 0.15)
ipmc 	<- FindClusters( ipmc, reuse.SNN = TRUE, resolution = clResolution, random.seed = clustSeed, algorithm = 2)
ipmc	<- ValidateClusters( ipmc, top.genes = 7, pc.use = 1:comps, acc.cutoff = 0.9, min.connectivity = 0.05, verbose = TRUE)

clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)
ipmc 	<- BuildClusterTree( ipmc, pcs.use = 1:comps, do.plot = FALSE, do.reorder = TRUE) #This functions renames clusters, so we need to assign cluster types again

png( file.path( resolDir, "ClusterTreePCASpace.png"))
	PlotClusterTree( ipmc)
dev.off()

clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)

#End construct TSNEPlot for clustered data
png( file.path( resolDir, paste0("TSNEClustersPCASpace_c", comps, "_res", clResolution, ".png")))
	TSNEPlot( ipmc, colors.use = setClusterColors( clTypes))
dev.off()

source("R/makeLineagesPCASpace.r")
makeLineagesPCASpace( ipmc, comps, clResolution)

#remove values, that are too close to zero
noiseTol	<- log2(19)
ipmc@data	<- apply( ipmc@data, c(1,2), function(x) if(x>noiseTol) x else 0)

png( file.path( resolDir, paste0("DotPlotPCASpace_c", comps, "_res", clResolution,".png")), width = 800, height = 600)
	DotPlot(ipmc, genes.plot = rownames(ipmc@data), x.lab.rot = TRUE, dot.scale = 5, plot.legend = TRUE, dot.min = 0, scale.by = "radius")
dev.off()


}

