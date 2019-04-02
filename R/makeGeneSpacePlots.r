#this snippet is use to make TSNE and PCA plots for initial cell types and clusters in the gene space
#it requires to run first seuratNorm.r, which among others sets the appropriate experimentType ( allCells, WT, WT_sox10)
#in ipmc@project.name

source("R/getClusterTypes.r")
source("R/calcTSNEGeneSpace.r")
source("R/setCellTypeColors.r")
source("R/setClusterColors.r")
source("R/plotInitCellTypePCAs.r")
source("R/calcUMAPGeneSpace.r")

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
ipmc <- calcTSNEGeneSpace( ipmc, TSNESeed)

UMAPSeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( geneSpacePlotDir, "UMAPSeed.txt"), TSNESeed, "\n")
ipmc <- calcUMAPGeneSpace( ipmc, UMAPSeed)

ipmc <- SetAllIdent( ipmc, id = 'originalCellTypes')
png( file.path( geneSpacePlotDir, "TSNEInitCellTypes.png"))
	TSNEPlot( ipmc, colors.use = setCellTypeColors( ipmc))
dev.off()

png( file.path( geneSpacePlotDir, "UMAPInitCellTypes.png"))
	DimPlot(object = ipmc, reduction.use = 'umap', pt.size = 1, cols.use = setCellTypeColors( ipmc))
dev.off()


#ipmc <- RunPCA( ipmc, pc.genes = rownames(ipmc@data), weight.by.var = FALSE)
#plotInitCellTypePCAs(ipmc, 6, geneSpacePlotDir)	#Plot PCA diagrams with cell colors, uses its own directorial structure

#remove values, that are too close to zero
noiseTol		<- log2(19)
ipmcDenoise		<- ipmc
ipmcDenoise@data	<- apply( ipmcDenoise@data, c(1,2), function(x) if(x>noiseTol) x else 0)

png( file.path( geneSpacePlotDir, "DotPlotGeneSpace.png"), width = 800, height = 600)
	DotPlot(ipmcDenoise, genes.plot = rownames(ipmcDenoise@data), x.lab.rot = TRUE, dot.scale = 5, plot.legend = TRUE, dot.min = 0, scale.by = "radius")
dev.off()





