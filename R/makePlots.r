#this snippet is use to make plots for 

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
ipmc 	<- FindClusters( ipmc, genes.use = rownames(ipmc@data), k.param = 10, resolution = 0.8)
ipmc 	<- BuildClusterTree( ipmc, genes.use = rownames(ipmc@data), do.plot = FALSE, do.reorder = TRUE)
clTypes <- getClusterTypes(ipmc@ident)
ipmc@cluster.tree[[1]]$tip.label <- names(clTypes)

png( file.path( initialCellTypeDir, "ClusterTreeGeneSpace.png"))
	PlotClusterTree( ipmc)
dev.off()

png( file.path( initialCellTypeDir, "TSNEClusters_geneSpace.png"))
	TSNEPlot( ipmc, colors.use = setClusterColors( clTypes))
dev.off()


