
#Clustering with dimenstion reduction
clResolution	<- 1.
resolDir	<- file.path( compsDir, paste0( "resol_", clResolution))
dir.create( resolDir, showWarnings = FALSE)

clustSeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( resolDir, "clustSeed.txt"), TSNESeed, "\n")

ipmc	<- BuildSNN( ipmc, dims.use = 1:6, prune.SNN = 0.15)
ipmc 	<- FindClusters( ipmc, reuse.SNN = TRUE, resolution = clResolution, random.seed = clustSeed)
ipmc	<- ValidateClusters( ipmc, top.genes = 7, pc.use = 1:6, acc.cutoff = 0.005, min.connectivity = 0.05, verbose = TRUE)

clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)
ipmc 	<- BuildClusterTree( ipmc, pcs.use = 1:comps, do.plot = FALSE, do.reorder = TRUE) #This functions renames clusters, so we need to assign cluster types again

png( file.path( resolDir, "ClusterTreePCASpace.png"))
	PlotClusterTree( ipmc)
dev.off()

clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)

#End construct TSNEPlot for clustered data
png( file.path( resolDir, "TSNEClusters_PCASpace.png"))
	TSNEPlot( ipmc, colors.use = setClusterColors( clTypes))
dev.off()


