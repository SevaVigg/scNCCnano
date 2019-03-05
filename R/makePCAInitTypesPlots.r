#this snippet useis to make TSNE and PCA plots for initial cell types
#
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
compsCorr	<- if (comps == 45) comps-1 else comps
compsDir 	<- file.path( pcaPlotDir, paste0("comps", comps))
dir.create(compsDir, showWarnings = FALSE)

ipmc	 	<- RunPCA(ipmc, pc.genes = rownames( ipmc@data), pcs.compute = compsCorr, do.print = FALSE, fastpath = FALSE)

#Tree of initial cell types with PCA dimension reduction but without clustering
ipmc 		<- BuildClusterTree( ipmc, pcs.use = 1:compsCorr, do.plot = FALSE, do.reorder = FALSE)
png( file.path( compsDir, "InitCellTypePCATree.png"))
	PlotClusterTree( ipmc)
dev.off()

#Now calculate TSNE for initial cell types with dimension reduction
TSNESeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( compsDir, "TSNESeed.txt"), TSNESeed, "\n")
ipmc <- calcTSNE_PCASpace( ipmc, compsCorr, TSNESeed)
png( file.path( compsDir, "tSNE_PCA_initialCellTypes.png"))
	TSNEPlot( ipmc, colors.use = setCellTypeColors( ipmc))
dev.off()





