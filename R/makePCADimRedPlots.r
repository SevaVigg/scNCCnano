#this snippet useis to make TSNE and PCA plots for initial cell types
#ipmc must contain the precomputed PCAs
#it requires to run first seuratNorm.r, which among others sets the appropriate experimentType ( allCells, WT, WT_sox10)
#in ipmc@project.name

source("R/getClusterTypes.r")
source("R/calcTSNE_PCASpace.r")
source("R/setCellTypeColors.r")
source("R/setClusterColors.r")
source("R/plotInitCellTypePCAs.r")

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypeDir <- file.path(plotDir, ipmc@project.name)
dir.create( experimentTypeDir, showWarnings = FALSE)

pcaPlotDir 	<- file.path( experimentTypeDir, "PCAdimReduction")
dir.create( pcaPlotDir, showWarnings = FALSE)

comps		<- 6

compsDir 	<- file.path( pcaPlotDir, paste0("comps", 6))
dir.create(compsDir, showWarnings = FALSE)

ipmc	 	<- RunPCA(ipmc, pc.genes = rownames( ipmc@data), pcs.compute = comps, do.print = FALSE)

ipmc <- BuildClusterTree( ipmc, pcs.use = 1:comps, do.plot = FALSE, do.reorder = FALSE)
png( file.path( compsDir, "InitCellTypePCATree.png"))
	PlotClusterTree( ipmc)
dev.off()

TSNESeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( compsDir, "TSNESeed.txt"), TSNESeed, "\n")

ipmc <- calcTSNE_PCASpace( ipmc, comps, TSNESeed)
png( file.path( compsDir, "tSNE_PCA_initialCellTypes.png"))
	TSNEPlot( ipmc, colors.use = setCellTypeColors( ipmc))
dev.off()

#clusters in initial gene space
ipmc 	<- FindClusters( ipmc, dims.use = 1:comps, resolution = 0.8, prune.SNN = 0.15)
clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)
ipmc 	<- BuildClusterTree( ipmc, pcs.use = 1:comps, do.plot = FALSE, do.reorder = TRUE) #This functions renames clusters, so we need to assign cluster types again

png( file.path( compsDir, "ClusterTreePCASpace.png"))
	PlotClusterTree( ipmc)
dev.off()

clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)

png( file.path( compsDir, "TSNEClusters_PCASpace.png"))
	TSNEPlot( ipmc, colors.use = setClusterColors( clTypes))
dev.off()




