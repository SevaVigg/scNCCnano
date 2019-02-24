#this snippet is use to make TSNE and PCA plots for initial cell types
#it requires to run first seuratNorm.r, which among others sets the appropriate experimentType ( allCells, WT, WT_sox10)
#in ipmc@project.name

source("R/getClusterTypes.r")
source("R/calcTSNE_GeneSpace.r")
source("R/setCellTypeColors.r")
source("R/setClusterColors.r")
source("R/plotInitCellTypePCAs.r")

resDir		<- file.path(getwd(), "Res")
plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypeDir <- file.path(plotDir, ipmc@project.name)
dir.create( experimentTypeDir, showWarnings = FALSE)

initialCellTypeDir <- file.path( experimentTypeDir, "initialCellTypes")
dir.create( initialCellTypeDir, showWarnings = FALSE)

ipmc <- BuildClusterTree( ipmc, genes.use = rownames(ipmc@data), do.plot = FALSE, do.reorder = FALSE)
png( file.path( initialCellTypeDir, "InitCellTypeTree.png"))
	PlotClusterTree( ipmc)
dev.off()

TSNESeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( initialCellTypeDir, "TSNESeed.txt"), TSNESeed, "\n")
ipmc <- calcTSNE_GeneSpace( ipmc, TSNESeed)
png( file.path( initialCellTypeDir, "TSNEinitialCellTypes.png"))
	TSNEPlot( ipmc, colors.use = setCellTypeColors( ipmc))
dev.off()

ipmc <- RunPCA( ipmc, pc.genes = rownames(ipmc@data), weight.by.var = FALSE)
plotInitCellTypePCAs(ipmc, 6, initialCellTypeDir)	#Plot PCA diagrams with cell colors, uses its own directorial structure

#clusters in initial gene space
ipmc 	<- FindClusters( ipmc, genes.use = rownames(ipmc@data), k.param = 10, resolution = 0.8, prune.SNN = 0.15)
clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)
ipmc 	<- BuildClusterTree( ipmc, genes.use = rownames(ipmc@data), do.plot = FALSE, do.reorder = TRUE) #This functions renames clusters, so we need to assign cluster types again

png( file.path( initialCellTypeDir, "ClusterTreeGeneSpace.png"))
	PlotClusterTree( ipmc)
dev.off()

clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)

png( file.path( initialCellTypeDir, "TSNEClusters_geneSpace.png"))
	TSNEPlot( ipmc, colors.use = setClusterColors( clTypes))
dev.off()




