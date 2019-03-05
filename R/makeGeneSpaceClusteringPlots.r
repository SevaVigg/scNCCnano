#Clustering in the initial gene space reduction

clResolution	<- 1 

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypeDir <- file.path(plotDir, ipmc@project.name)
dir.create( experimentTypeDir, showWarnings = FALSE)

geneSpacePlotDir <- file.path( experimentTypeDir, "geneSpace")
dir.create( geneSpacePlotDir, showWarnings = FALSE)

resolDir	<- file.path( geneSpacePlotDir, paste0( "resol_", clResolution))
dir.create( resolDir, showWarnings = FALSE)

clustSeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( resolDir, "clustSeed.txt"), TSNESeed, "\n")

ipmc	<- BuildSNN( ipmc, genes.use = rownames(ipmc@data), k.param = 10, prune.SNN = 0.15)
ipmc 	<- FindClusters( ipmc, reuse.SNN = TRUE, resolution = clResolution, random.seed = clustSeed)
#ipmc	<- ValidateClusters( ipmc, top.genes = 7, pc.use = NULL, acc.cutoff = 0.005, min.connectivity = 0.05, verbose = TRUE)

clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)
ipmc 	<- BuildClusterTree( ipmc, genes.use = rownames(ipmc@data), do.plot = FALSE, do.reorder = TRUE) #This functions renames clusters, so we need to assign cluster types again

png( file.path( resolDir, "ClusterTreeGeneSpace.png"))
	PlotClusterTree( ipmc)
dev.off()

clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)

#End construct TSNEPlot for clustered data
png( file.path( resolDir, "TSNEClusters_GeneSpace.png"))
	TSNEPlot( ipmc, colors.use = setClusterColors( clTypes))
dev.off()


