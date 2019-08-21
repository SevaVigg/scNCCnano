#this snippet useis to make TSNE and PCA plots for initial cell types
#
#it requires to run first seuratNorm.r, which among others sets the appropriate experimentType ( allCells, WT, WT_sox10)
#in ipmc@project.name

makePCAInitTypesPlots <- function(ipmc, comps){

source("R/getClusterTypes.r")
source("R/calcTSNE_PCASpace.r")
source("R/setCellTypeColors.r")
source("R/setClusterColors.r")
source("R/plotInitCellTypePCAs.r")
source("R/calcUMAP_PCASpace.r")

resDir		<- file.path(getwd(), "Res")

cat(resDir)

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypeDir <- file.path(plotDir, ipmc@project.name)
dir.create( experimentTypeDir, showWarnings = FALSE)

pcaPlotDir 	<- file.path( experimentTypeDir, "PCAdimReduction")
dir.create( pcaPlotDir, showWarnings = FALSE)

#comps		<- 17
compsCorr	<- if (comps == 45) comps-1 else comps
compsDir 	<- file.path( pcaPlotDir, paste0("comps", comps))
dir.create(compsDir, showWarnings = FALSE)

ipmc		<- SetAllIdent( ipmc, id = "originalCellTypes")
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

UMAPSeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( compsDir, "UMAPSeed.txt"), UMAPSeed, "\n")
ipmc <- calcUMAP_PCASpace( ipmc, compsCorr, TSNESeed, 2)
png( file.path( compsDir, "UMAP_PCA_initialCellTypes.png"))
	DimPlot( ipmc, reduction.use = "umap", cols.use = setCellTypeColors( ipmc))
dev.off()



#remove values, that are too close to zero to make zero-expressed genes for DotPlot; we need to keep initial ipmc for good
noiseTol	<- log2(19)

ipmcDenoised		<- ipmc #for not spoiling initial clustering data
ipmcDenoised@data	<- apply( ipmcDenoised@data, c(1,2), function(x) if(x>noiseTol) x else 0)

png( file.path( compsDir, paste0("DotPlotInitTypesPCASpace_c", comps, ".png")), width = 800, height = 600)
	DotPlot(ipmcDenoised, genes.plot = rownames(ipmcDenoised@data), x.lab.rot = TRUE, dot.scale = 5, plot.legend = TRUE, dot.min = 0, scale.by = "radius")
dev.off()

#source("R/makePCAClusteringPlots.r")
#for (clResolution in seq(0.8)) {makePCAClusteringPlots( ipmc, clResolution)} 

return(ipmc)
}




