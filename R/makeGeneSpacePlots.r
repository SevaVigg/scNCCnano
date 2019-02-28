#this snippet is use to make TSNE and PCA plots for initial cell types and clusters in the gene space
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

geneSpacePlotDir <- file.path( experimentTypeDir, "geneSpace")
dir.create( geneSpacePlotDir, showWarnings = FALSE)

ipmc <- BuildClusterTree( ipmc, genes.use = rownames(ipmc@data), do.plot = FALSE, do.reorder = FALSE)
png( file.path( geneSpacePlotDir, "InitCellTypeTree.png"))
	PlotClusterTree( ipmc)
dev.off()

TSNESeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( geneSpacePlotDir, "TSNESeed.txt"), TSNESeed, "\n")
ipmc <- calcTSNE_GeneSpace( ipmc, TSNESeed)
png( file.path( geneSpacePlotDir, "TSNEinitialCellTypes.png"))
	TSNEPlot( ipmc, colors.use = setCellTypeColors( ipmc))
dev.off()

ipmc <- RunPCA( ipmc, pc.genes = rownames(ipmc@data), weight.by.var = FALSE)
plotInitCellTypePCAs(ipmc, 6, geneSpacePlotDir)	#Plot PCA diagrams with cell colors, uses its own directorial structure

#clusters in initial gene space
ipmc 	<- FindClusters( ipmc, genes.use = rownames(ipmc@data), k.param = 10, resolution = 0.8, prune.SNN = 0.15)
clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)
ipmc 	<- BuildClusterTree( ipmc, genes.use = rownames(ipmc@data), do.plot = FALSE, do.reorder = TRUE) #This functions renames clusters, so we need to assign cluster types again

png( file.path( geneSpacePlotDir, "ClusterTreeGeneSpace.png"))
	PlotClusterTree( ipmc)
dev.off()

clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)

png( file.path( geneSpacePlotDir, "tSNEClustersGeneSpace.png"))
	TSNEPlot( ipmc, colors.use = setClusterColors( clTypes))
dev.off()




