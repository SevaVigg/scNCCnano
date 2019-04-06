# Making seurat plots from SC3 clusters

require(SC3)
require(SingleCellExperiment)

source("R/getClusterTypes.r")
source("R/setClusterColors.r")
source("R/calcTSNEGeneSpace.r")
source("R/calcUMAPGeneSpace.r")
source("R/setCellTypeColors.r")
source("R/makeFeaturePlots.r")


SC3Clusters <- colData(ipmc_sce)[, grep("clusters", colnames(colData(ipmc_sce)))]

SC3ClustDf <- as.data.frame(sapply(SC3Clusters, function(x) as.numeric(levels(x)[x])))
rownames(SC3ClustDf) <- rownames(SC3Clusters)

IPCls 		<- sapply(SC3ClustDf[grep("I", rownames(SC3ClustDf)), ], table, simplify = FALSE)
MCCls 		<- sapply(SC3ClustDf[grep("M", rownames(SC3ClustDf)), ], table, simplify = FALSE)
TlCls 		<- sapply(SC3ClustDf[grep("Tl", rownames(SC3ClustDf)), ], table, simplify = FALSE)
SoxCls 		<- sapply(SC3ClustDf[grep("sox", rownames(SC3ClustDf)), ], table, simplify = FALSE)
mitfaCls 	<- sapply(SC3ClustDf[grep("mitfa", rownames(SC3ClustDf)), ], table, simplify = FALSE)
ClCls 		<- sapply(SC3ClustDf, table, simplify = FALSE)

argmax 	<- function( x){ which(x == max(x))}

IPShrs		<- mapply(function(x,y) x[names(argmax(x))]/y[names(argmax(x))], IPCls, ClCls)
MCShrs		<- mapply(function(x,y) x[names(argmax(x))]/y[names(argmax(x))], MCCls, ClCls)
TlShrs		<- mapply(function(x,y) x[names(argmax(x))]/y[names(argmax(x))], TlCls, ClCls)
soxShrs		<- mapply(function(x,y) x[names(argmax(x))]/y[names(argmax(x))], SoxCls, ClCls)
mitfaShrs	<- mapply(function(x,y) x[names(argmax(x))]/y[names(argmax(x))], mitfaCls, ClCls)

source("R/seuratNorm.r")

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypeDir <- file.path(plotDir, ipmc@project.name)
dir.create( experimentTypeDir, showWarnings = FALSE)

SC3ClustNumDir <- file.path( experimentTypeDir, "geneSpace")
dir.create( SC3ClustNumDir, showWarnings = FALSE)

SC3ClustNum	<- 15

SC3ClustDir 	<- file.path( SC3ClustNumDir, "SC3clsuts")
dir.create( SC3ClustDir, showWarnings = FALSE)

SC3ClustNumDir	<- file.path( SC3ClustDir, paste0( "Clust_", SC3ClustNum))
dir.create( SC3ClustNumDir, showWarnings = FALSE)

ipmc@ident 		<- SC3Clusters[[paste0("sc3_", SC3ClustNum, "_clusters")]] 
names(ipmc@ident) 	<- rownames(SC3Clusters)

clTypes 		<- getClusterTypes(ipmc)
levels(ipmc@ident) 	<- names(clTypes)

ipmc <- StashIdent( ipmc, save.name = paste0("SC3_", SC3ClustNum, "_clusters"))

ipmc 	<- BuildClusterTree( ipmc, genes.use = rownames(ipmc@data), do.plot = FALSE, do.reorder = TRUE) #This functions renames clusters, so we need to assign cluster types again

png( file.path( SC3ClustNumDir, "ClusterTreeGeneSpace.png"))
	PlotClusterTree( ipmc)
dev.off()

clTypes <- getClusterTypes(ipmc)
levels(ipmc@ident) <- names(clTypes)

TSNESeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( SC3ClustNumDir, "TSNESeed.txt"), TSNESeed, "\n")
ipmc <- calcTSNEGeneSpace( ipmc, TSNESeed)

UMAPSeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( SC3ClustNumDir, "UMAPSeed.txt"), TSNESeed, "\n")
ipmc <- calcUMAPGeneSpace( ipmc, UMAPSeed)

png( file.path( SC3ClustNumDir, "UMAPSC3ClusterTypes.png"))
	DimPlot(object = ipmc, reduction.use = 'umap', pt.size = 1, cols.use = setClusterColors( ipmc))
dev.off()

FeaturePlotsDir <- file.path( SC3ClustNumDir, "FeaturePlots")
dir.create( FeaturePlotsDir, showWarnings = FALSE)

for (gene in rownames(ipmc@data)){
	png( file.path( FeaturePlotsDir, paste0( gene, ".png")))  
		FeaturePlot( ipmc, gene, cols.use = c("blue", "yellow"), reduction.use = "umap" ) 
	dev.off()
}
ViolinPlotsDir <- file.path( SC3ClustNumDir, "ViolinPlots")
dir.create( ViolinPlotsDir, showWarnings = FALSE)

for ( gene in rownames(ipmc@data)){
	png( file.path( ViolinPlotsDir, paste0( gene, ".png")))  
		VlnP <- VlnPlot( ipmc, gene, do.return = TRUE)
		plot(VlnP) 
	dev.off()
}
	

#remove values, that are too close to zero
noiseTol		<- log2(19)
ipmcDenoise		<- ipmc
ipmcDenoise@data	<- apply( ipmcDenoise@data, c(1,2), function(x) if(x>noiseTol) x else 0)

png( file.path( SC3ClustNumDir, "DotPlotGeneSpace.png"), width = 800, height = 600)
	DotPlot(ipmcDenoise, genes.plot = rownames(ipmcDenoise@data), x.lab.rot = TRUE, dot.scale = 5, plot.legend = TRUE, dot.min = 0, scale.by = "radius")
dev.off()

#make Lineages

coordsMD	<- as.matrix( t(ipmc@data))
coordsMD	<- coordsMD + jitter(coordsMD)
#these are initial coordinates, jitterred to avoid singular values

source("R/createSlingShotObject.r")
ipmcSling 	<- createSlingShotObject( coordsMD, ipmc)

source("R/plot2DallLineages.r")

png( file.path( SC3ClustNumDir, "clustersLineagesTSNE.png"))
plot2DallLineages( ipmcSling, ipmc, "tsne")
dev.off() 

png( file.path( SC3ClustNumDir, "clustersLineagesUMAP.png"))
plot2DallLineages( ipmcSling, ipmc, "umap")
dev.off() 



#ipmc <- RunPCA( ipmc, pc.genes = rownames(ipmc@data), weight.by.var = FALSE)
#plotInitCellTypePCAs(ipmc, 6, SC3ClustNumDir)	#Plot PCA diagrams with cell colors, uses its own directorial structure

#remove values, that are too close to zero
noiseTol		<- log2(19)
ipmcDenoise		<- ipmc
ipmcDenoise@data	<- apply( ipmcDenoise@data, c(1,2), function(x) if(x>noiseTol) x else 0)

png( file.path( SC3ClustNumDir, "DotPlotSC3ClustGeneSpace.png"), width = 800, height = 600)
	DotPlot(ipmcDenoise, genes.plot = rownames(ipmcDenoise@data), x.lab.rot = TRUE, dot.scale = 5, plot.legend = TRUE, dot.min = 0, scale.by = "radius")
dev.off()



ipmc <- SetAllIdent( ipmc, id = 'originalCellTypes')
png( file.path( SC3ClustNumDir, "TSNEInitCellTypes.png"))
	TSNEPlot( ipmc, colors.use = setCellTypeColors( ipmc))
dev.off()

png( file.path( SC3ClustNumDir, "UMAPInitCellTypes.png"))
	DimPlot(object = ipmc, reduction.use = 'umap', pt.size = 1, cols.use = setCellTypeColors( ipmc))
dev.off()


#ipmc <- RunPCA( ipmc, pc.genes = rownames(ipmc@data), weight.by.var = FALSE)
#plotInitCellTypePCAs(ipmc, 6, SC3ClustNumDir)	#Plot PCA diagrams with cell colors, uses its own directorial structure

#remove values, that are too close to zero
noiseTol		<- log2(19)
ipmcDenoise		<- ipmc
ipmcDenoise@data	<- apply( ipmcDenoise@data, c(1,2), function(x) if(x>noiseTol) x else 0)

png( file.path( SC3ClustNumDir, "DotPlotGeneSpace.png"), width = 800, height = 600)
	DotPlot(ipmcDenoise, genes.plot = rownames(ipmcDenoise@data), x.lab.rot = TRUE, dot.scale = 5, plot.legend = TRUE, dot.min = 0, scale.by = "radius")
dev.off()

ipmc <- SetAllIdent( ipmc, id = paste0("SC3_", SC3ClustNum, "_clusters"))



