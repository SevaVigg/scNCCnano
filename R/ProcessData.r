#The main process of treatment NCC nanostring data; execute all other scripts and makes all plots
#
# Finished by Vsevolod J. Makeev 15.05.2021, developed 2017 - 2021 


require( tidyr )
require( ape)
require( ComplexHeatmap)

source("R/ReadSourceFiles.r")
source("R/writeDupFile.r")

plotDPI		<- 600

#in the introduction we make the directory structure
workDir <- getwd()

rawPath 	 <- file.path( workDir, "SourceData")
resDir		 <- file.path( workDir, "Res")

dir.create( resDir, showWarnings = FALSE)

initialTablesPath 	<- file.path( resDir, "InitialTables")
dir.create( initialTablesPath, showWarnings = FALSE)

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

PCAPlotDir	<- file.path( plotDir, "PCAPlots")
dir.create( PCAPlotDir, showWarnings = FALSE)

heatMapDir <- file.path( plotDir, "heatMaps")
dir.create( heatMapDir, showWarnings = FALSE)

dotPlotDir <- file.path( plotDir, "dotPlots")
dir.create( dotPlotDir, showWarnings = FALSE)

clusterDataDir		<- file.path(resDir, "clusterData")
dir.create( clusterDataDir, showWarnings = FALSE)

clusterPlotDir		<- file.path( plotDir, "clusterPlots")
dir.create( clusterPlotDir, showWarnings = FALSE)

clusterTreeDir		<- file.path( plotDir, "clusterTrees")
dir.create( clusterTreeDir, showWarnings = FALSE)

vlnPlotDir <- file.path( plotDir, "vlnPlots")
dir.create( vlnPlotDir, showWarnings = FALSE)

#end of the introduction


#Read and process files
CellTable 	<- ReadSourceFiles(rawPath)

write.table(CellTable$Genes, file = file.path(initialTablesPath, "CellTableUnDup_ET.csv"), sep = "\t")
write.table(CellTable$Cells, file = file.path(initialTablesPath, "CellTableUnDup_CD.csv"), sep = "\t")

#Find duplicates

dupTable  	<- findDuplicated(CellTable)						#Save duplicated entries
writeDupFile(CellTable, dupTable, initialTablesPath)						#my function

			
CellTable$Genes 	<- CellTable$Genes[-dupTable[2,]]
CellTable$Cells 	<- CellTable$Cells[-dupTable[2,]]

genesMissing_I <- which(!complete.cases(CellTable$Genes))
write.table( CellTable$Probes[genesMissing_I,], file = file.path( initialTablesPath, "missingGenes.csv") )

cellsWMissingGenes_I <- unlist( lapply(genesMissing_I, function(gene_I) {cell_I <- which(is.na(CellTable$Genes[gene_I,]))
				geneFileName <- file.path( initialTablesPath, paste0("Missing_", CellTable$Probes[gene_I, "Gene.Name"], "_cells.csv"))
				write.table( CellTable$Cells[,cell_I], file = geneFileName )
					return( cell_I)}))
#remove cells with missing genes

CellTable$Genes	<- CellTable$Genes[-cellsWMissingGenes_I]
CellTable$Cells <- CellTable$Cells[-cellsWMissingGenes_I]

cellNames	<- paste0(CellTable$Cells["CellType", ], "_", "C", 1:ncol(CellTable$Genes), "_",  CellTable$Cells["hpf", ])
rownames( CellTable$Genes) <- make.names( rownames(CellTable$Genes))
colnames( CellTable$Genes) <- cellNames
rownames( CellTable$Cells) <- make.names( rownames(CellTable$Cells))
colnames( CellTable$Cells) <- cellNames
CellTable$Probes$Gene.Name <- make.names( CellTable$Probes$Gene.Name)


#Write deduplicated results

write.csv( CellTable$Genes,  file = file.path( initialTablesPath, "expressionTableDedup.csv"   ) )
write.csv( CellTable$Cells,  file = file.path( initialTablesPath, "cellDescripitonsDedup.csv"  ) )
write.csv( CellTable$Probes, file = file.path( initialTablesPath, "ProbesDescripitonsDedup.csv") )

#qualityControl 
source("R/qualityControl600dpi.r")

#imputation
source("R/makeScTables.r")

#start Seurat analysis

source("R/seuratNorm.r")
seuratWT 	<- seuratNorm("WT")
seuratAll	<- seuratNorm("allCells")

seuratAll <- SetAllIdent( seuratAll, id = "genCellTypeIdent")
seuratWT <- SetAllIdent( seuratWT, id = "genCellTypeIdent")

allGenes <- rownames( seuratWT@data)

#gene heatmaps

heatMapDir100dpi	<- file.path( heatMapDir, "100dpi")
dir.create( heatMapDir100dpi, showWarnings = FALSE)

heatMapDir600dpi	<- file.path( heatMapDir, "600dpi")
dir.create( heatMapDir600dpi, showWarnings = FALSE)

#gene covariance
source( "R/drawGeneCovarHeatMap.r")

heatMapHeight 	<- 10
heatMapWidth 	<- 10
Margin		<- 2

if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "geneCovarHeatMap.png"),
	height = heatMapHeight + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawGeneCovarHeatMap( seuratWT, heatMapHeight, heatMapWidth))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "geneCovarHeatMap.png"),
	height = heatMapHeight + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawGeneCovarHeatMap( seuratWT, heatMapHeight, heatMapWidth))
dev.off()}

#gene biclusters
source( "R/drawBiclustHeatMap.r")

heatMapHeight 	<- 7
heatMapWidth 	<- 10
Margin		<- 2

if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "biclustWTHeatMap.png"),
	height = heatMapHeight + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawBiclustHeatMap( seuratWT, heatMapHeight, heatMapWidth))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "biclustWTHeatMap.png"),
	height = heatMapHeight + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawBiclustHeatMap( seuratWT, heatMapHeight, heatMapWidth))
dev.off()}

source("R/makeDotPlot.r")

nClust = length( levels( seuratWT@ident))
makeDotPlot( seuratWT, balanced = TRUE, nLines = length( levels( seuratWT@ident)), plotDPI = plotDPI, orientation = "landscape", name = "WTdotPlot")

#init cell type distributions, also prepares seuratWT and seuratAll objects
#source("R/makeGeneSpacePlots.r")  

#PCA analysis
source("R/makeInitCellTypePCAPlots.r")

makeInitCellTypePCAPlots( seuratWT,  nComps = 5, plotDPI = plotDPI, name = "WT_PCAcomps")
makeInitCellTypePCAPlots( seuratAll, nComps = 5, plotDPI = plotDPI, name = "All_PCAcomps")

#find best parameters for UMAP clustering
source("R/findBestUmapClusters.r")
Spread		<- 10

#bestUmap 	<- findBestUmapClusters( seuratWT, Spread)
#save( bestUmap, file = file.path( clusterDataDir, "bestUmap_regular_renamed.rObj"))

#OR load the data file contains validated clusters
load(file = file.path( clusterDataDir, "bestUmap_regular_renamed.rObj"))
#


bestDim		<- dim(bestUmap@dr$umap@cell.embeddings)[2]

source("R/getFinalClusterTypes.r")
levels(bestUmap@ident) <- names(getFinalClusterTypes( bestUmap))

#build coarse grain cluster tree
bestUmap <- BuildClusterTree( bestUmap, pcs.use = 1:7, do.reorder = TRUE, reorder.numeric = TRUE, do.plot = FALSE)

levels(bestUmap@ident) <- names(getFinalClusterTypes( bestUmap))

bestUmap <- BuildClusterTree( bestUmap, pcs.use = 1:7, do.reorder = FALSE, reorder.numeric = FALSE, do.plot = FALSE)

source("R/cosineClusterDist.r")

clusterTreeDir100dpi	<- file.path( clusterTreeDir, "100dpi")
dir.create( clusterTreeDir100dpi, showWarnings = FALSE)

clusterTreeDir600dpi	<- file.path( clusterTreeDir, "600dpi")
dir.create( clusterTreeDir600dpi, showWarnings = FALSE)

source("R/plotClusterTree.r")
plotClusterTree( bestUmap, plotDPI = plotDPI, treeName = "coarseGrainClusterTree")

makeDotPlot( bestUmap, balanced = TRUE, nLines = length( levels( bestUmap@ident)), plotDPI = plotDPI, orientation = "landscape", name = "coarseGrainDotPlot")
bestUmap <- StashIdent( bestUmap, save.name = "bestClustersIdent")

#source("R/make2Dmap.r")
#umap2Dinit <- make2Dmap( seuratWT)
#save( file = file.path( clusterDataDir, "visualisation2Dumap_renamed.rObj"), umap2Dinit)

# OR
load( file.path( clusterDataDir, "visualisation2Dumap_renamed.rObj")); 

source("R/setClusterColors.r")
source("R/calcUmapGeneSpace.r")
source("R/makeDimPlot.r")
makeDimPlot( umap2Dinit, dimRed = "umap", name = "initCellPlot", orientation = "landscape", plotDPI = plotDPI)

umapCondRes <- calcUmapGeneSpace( seuratAll, seurWT = umap2Dinit, experimentType = "allCondWT", Dim = 2, minDist = 4, mySpread = 9, myNeighbors = 25L, myMetric = "cosine")
makeDimPlot( umapCondRes$All, dimRed = "umap", name = "initCellPlotMut", orientation = "landscape", plotDPI = plotDPI)


#now validate clusters
valCutoff 		<- 0.92
valCutoffIdentName	<- paste0( "cutoffIdent", valCutoff)

valUmap <- ValidateClusters( bestUmap, pc.use = 1:8, top.genes = 4, min.connectivity = 0.005, acc.cutoff = 0.75)

for (i in seq(0.76, valCutoff, 0.01)){ cat(i, "\n")
   valUmap <- ValidateClusters( valUmap, pc.use = 1:8, top.genes = 4, min.connectivity = 0.005, acc.cutoff = i)
}

valUmap <- BuildClusterTree( valUmap, pcs.use = 1:8, do.reorder = TRUE, reorder.numeric = TRUE, do.plot = TRUE)
levels(valUmap@ident) <- names( getFinalClusterTypes( valUmap))
valUmap <- BuildClusterTree( valUmap, pcs.use = 1:8, do.reorder = TRUE, reorder.numeric = FALSE, do.plot = TRUE)

png( file.path( clusterTreeDir, "validatedClusterTree.png"))
	PlotClusterTree( valUmap) 
dev.off()

valUmap 	<- StashIdent( valUmap, save.name = valCutoffIdentName)

#make cluster PCA plot
makeInitCellTypePCAPlots( valUmap, nComps = 5, plotDPI = plotDPI, name = "clusterPCAPlot")

#Plot special dotPlots for control cell types
makeDotPlot( valUmap, balanced = TRUE, nLines = length( levels( valUmap@ident)), plotDPI = plotDPI, orientation = "landscape", name = "valClusterDotPlot")

#make Cluster HeatMap

source("R/drawClusterHeatMap.r")
if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "clusterHeatMap.png"),
	height = heatMapHeight + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawClusterHeatMap( valUmap, heatMapHeight, heatMapWidth))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "clusterHeatMap.png"),
	height = heatMapHeight + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawClusterHeatMap( valUmap, heatMapHeight, heatMapWidth))
dev.off()}

umap2Dclust <- umap2Dinit
umap2Dclust@meta.data[, valCutoffIdentName] 	<- valUmap@meta.data[, valCutoffIdentName]
umap2Dclust					<- SetAllIdent( umap2Dclust, id = valCutoffIdentName)

#make cluster plot without pseudotime trajectories
makeDimPlot( umap2Dclust, dimRed = "umap", name = "clusterPlot", col = setClusterColors( umap2Dclust), orientation = "landscape", plotDPI = plotDPI)


#Now we have clustering and can make feature plots accompanied with clusters
source("R/makeFeaturePlots.r")
makeFeaturePlots( umap2Dclust, minCutoff = 3, "umap", plotDPI = plotDPI, name = "initFeaturePlots", orientation = "portrait")

source("R/createSlingShotObjects.r")
if (bestDim == 2){
	seur2D	<- valUmap
	slingUmapObjs <- createSlingShotObjects( 
		seurHD = valUmap, seur2D = seur2D, dimRed2D = "umap", 
		genesUseHD = allGenes, 
		startClust = "eHMP", endClust = c("I","M","X"), 
		distFun = cosineClusterDist)
}else{
	seur2D	<- umap2Dclust
	slingUmapObjs <- createSlingShotObjects( 
		seurHD = valUmap, seur2D = seur2D, dimRed2D = "umap", 
		genesUseHD = allGenes, 
		startClust = "eHMP", endClust = c("I","M","X"), 
		distFun = cosineClusterDist)
}

#we also need pseudotime trajectories
source("R/plot2DAllCurves_ggplot.r") #this file contains function "plot2DAllCurves" which uses ggplot rather than lines
plot2DAllCurves( seur2D, slingUmapObjs, dimRed2D = "umap", lineageToDrawEnds = c("M", "I"), cellsKeepThresh = 0.95, plotDPI = plotDPI, name = "slingClustersUmap")

source("R/calcTSNEGeneSpace.r")
valUmap   <- calcTSNEGeneSpace( valUmap );
save( file = file.path( clusterDataDir, "tsneValUmap.rObj"), valUmap)
#load( file = file.path( clusterDataDir, "tsneValUmap.rObj"))


slingTSNE <- createSlingShotObjects( seurHD = valUmap, seur2D = valUmap, dimRed2D = "tsne", 
	genesUseHD = allGenes,
	startClust = "eHMP", endClust = c("I","M","X"), 
	distFun = cosineClusterDist)
plot2DAllCurves( valUmap, slingTSNE, dimRed2D = "tsne", plotDPI = plotDPI, name = "slingShotClusterTSNEPlot")

#now study control cell type clusters. We start from melanocytes

Mcells 				<- WhichCells( valUmap, ident = "M")
seuratM				<- SubsetData( valUmap, ident.use = "M") 
cntrlM				<- grep( "M", Mcells, value = TRUE)
regM				<- setdiff( Mcells, cntrlM) 					
melIndexLine			<- c( rep( "contM", length( cntrlM)), rep( "regM", length( regM)))
names( melIndexLine)		<- c( cntrlM, regM)
seuratM@meta.data$Mcells 	<- melIndexLine[ rownames( seuratM@meta.data)]
seuratM <- SetAllIdent( seuratM, id = "Mcells")

makeDotPlot( seuratM, balanced = TRUE, nLines = length( levels( seuratM@ident))+2, plotDPI = plotDPI, orientation = "landscape", name = "melanoDotPlot")


#melanocyte subtypes
seuratMSubTypes <- FindClusters( seuratM, dims.use = 1:8, resolution = 1.0)
makeDotPlot( seuratMSubTypes, balanced = TRUE, nLines = length( levels( seuratMSubTypes@ident))+2, plotDPI = plotDPI, orientation = "landscape", name = "melSubTypesDotPlot")

#and melanocyte heatmap
if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "cntrMeloHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawBiclustHeatMap( seuratM, heatMapHeight, heatMapWidth, showCellNames = TRUE))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "cntrMeloHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawBiclustHeatMap( seuratM, heatMapHeight, heatMapWidth, showCellNames = TRUE))
dev.off()}


#For iridophores. 
#Plot special dotPlots for iridophores cell types

Icells 				<- WhichCells( valUmap, ident = "I")
seuratI				<- SubsetData( valUmap, ident.use = "I") 
cntrlI				<- grep( "I", Icells, value = TRUE)
regI				<- setdiff( Icells, cntrlI) 					
irdIndexLine			<- c( rep( "contI", length( cntrlI)), rep( "regI", length( regI)))
names( irdIndexLine)		<- c( cntrlI, regI)
seuratI@meta.data$Icells 	<- irdIndexLine[ rownames( seuratI@meta.data)]
seuratI <- SetAllIdent( seuratI, id = "Icells")

makeDotPlot( seuratI, balanced = TRUE, nLines = length( levels( seuratI@ident))+2, plotDPI = plotDPI, orientation = "landscape", name = "iridoDotPlot")

#For both. 
#Plot special dotPlots for both melanocytes and iridophores cell types

pigmCells 			<- WhichCells( valUmap, ident = c("I", "M"))
seuratPigm			<- SubsetData( valUmap, ident.use = c("I", "M")) 
cntrlI				<- grep( "I", pigmCells, value = TRUE)
regI				<- setdiff( WhichCells( seuratPigm, ident = "I"), cntrlI) 					
irdIndexLine			<- c( rep( "contI", length( cntrlI)), rep( "regI", length( regI)))
names( irdIndexLine)		<- c( cntrlI, regI)
cntrlM				<- grep( "M", pigmCells, value = TRUE)
regM				<- setdiff( WhichCells( seuratPigm, ident = "M"), cntrlM) 					
melIndexLine			<- c( rep( "contM", length( cntrlM)), rep( "regM", length( regM)))
names( melIndexLine)		<- c( cntrlM, regM)
pigmIndexLine			<- c( melIndexLine, irdIndexLine)
seuratPigm@meta.data$pigmCells 	<- pigmIndexLine[ rownames( seuratPigm@meta.data)]
seuratPigm 			<- SetAllIdent( seuratPigm, id = "pigmCells")
seuratPigm@ident 		<- ordered( seuratPigm@ident, levels = c("regM", "contM", "regI", "contI"))

makeDotPlot( seuratPigm, balanced = TRUE, nLines = length( levels( seuratPigm@ident))+1, plotDPI = plotDPI, orientation = "landscape", name = "pigmDotPlot")

#and the iridophore heatmap

if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "cntrIridoHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawBiclustHeatMap( seuratI, heatMapHeight, heatMapWidth, showCellNames = TRUE))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "cntrIridoHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawBiclustHeatMap( seuratI, heatMapHeight, heatMapWidth, showCellNames = TRUE))
dev.off()}

# Now we plot pseudotime heatmaps

source("R/getTargetCurve.r")
source("R/drawTargetHeatMapCells.r")

if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "iridoPseudoCellsHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawTargetHeatMapCells( valUmap, getTargetCurve( slingUmapObjs, target = "I"), heatMapHeight = heatMapHeight, heatMapWidth = heatMapWidth, reorderClusters = TRUE))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "iridoPseudoCellsHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawTargetHeatMapCells( valUmap, getTargetCurve( slingUmapObjs, target = "I"), heatMapHeight = heatMapHeight, heatMapWidth = heatMapWidth, reorderClusters = TRUE))
dev.off()}

#and now melanocytes
if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "melanoPseudoCellsHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawTargetHeatMapCells( valUmap, getTargetCurve( slingUmapObjs, target = "M"), heatMapHeight = heatMapHeight, heatMapWidth = heatMapWidth, reorderClusters = TRUE))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "melanoPseudoCellsHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawTargetHeatMapCells( valUmap, getTargetCurve( slingUmapObjs, target = "M"), heatMapHeight = heatMapHeight, heatMapWidth = heatMapWidth, reorderClusters = TRUE))
dev.off()}

#and for smoothed values, first iridophores
source("R/drawTargetHeatMapCurve.r")

if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "iridoTargetHeatMapSmooth.png"),
	height = heatMapHeight + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawTargetHeatMapCurve( valUmap, getTargetCurve( slingUmapObjs, target = "I"), heatMapHeight, heatMapWidth))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "iridoTargetHeatMapSmooth.png"),
	height = heatMapHeight + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2
)
	draw( drawTargetHeatMapCurve( valUmap,  getTargetCurve( slingUmapObjs, target = "I"), heatMapHeight, heatMapWidth))
dev.off()}

#and melanocytes

if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "melanoTargetHeatMapSmooth.png"),
	height = heatMapHeight + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawTargetHeatMapCurve( valUmap, getTargetCurve( slingUmapObjs, target = "M"), heatMapHeight, heatMapWidth))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "melanoTargetHeatMapSmooth.png"),
	height = heatMapHeight + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2
)
	draw( drawTargetHeatMapCurve( valUmap,  getTargetCurve( slingUmapObjs, target = "M"), heatMapHeight, heatMapWidth))
dev.off()}


#vlnPlots
source("R/makeVlnPlots.r")
makeVlnPlots( valUmap, plotDPI = plotDPI) 


#Now consider mutant cells
seuratMut 	<- SubsetData( seuratAll, ident.use = "sox10-")
seuratMut@project.name <- "sox10_mutants"

#Heatmap of mutant cells with ltk ordering
source("R/drawGeneSortHeatMap.r")
 
#geneSortList <- c("ltk", "foxo1a")

if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "ltkSortedMutantHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawGeneSortHeatMap( seuratMut,  heatMapHeight, heatMapWidth, showCellNames = FALSE))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "ltkSortedmutantHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawGeneSortHeatMap( seuratMut,  heatMapHeight, heatMapWidth, showCellNames = FALSE))
dev.off()}

	

#First make dotplot with mutants

makeDotPlot( seuratAll, balanced = TRUE, nLines = length( levels( seuratAll@ident)), plotDPI = plotDPI, orientation = "landscape", name = "WT_and_MutdotPlot")

#Now cluster mutant cells with best parameters

paramStr	<- grep( "Best", colnames( bestUmap@meta.data), value = TRUE)
paramVector	<- strsplit( paramStr, "_")[[1]] 

#bestDim		<- as.numeric( strsplit( grep( "D", paramVector, value = TRUE), "D")[[1]][2])
bestDim		<- 2   #there is a but in the parameter optimizer
bestMinDist	<- as.numeric( strsplit( grep( "minDist", paramVector, value = TRUE), "minDist")[[1]][2])
bestR		<- as.numeric( strsplit( grep( "R", paramVector, value = TRUE), "R")[[1]][2])

seuratMut		<- calcUmapGeneSpace( seuratMut, Dim = bestDim, myNeighbors = 25L, mySpread = Spread,  
				minDist = bestMinDist,  UMAPRandSeed = 42L, experimentType <- "allCells")$All

seuratMut 		<- FindClusters( seuratMut, reduction.type = 'umap', 
					dims.use = 1:bestDim, 
					k.param = 15, print.output = FALSE, force.recalc = TRUE, save.SNN = TRUE,
					resolution = bestR) 


makeDotPlot( seuratMut, balanced = TRUE, nLines = length( levels( seuratAll@ident)), plotDPI = plotDPI, orientation = "landscape", name = "mutDotPlot")

#and the heatmap

if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "clusterMutHeatMap.png"),
	height = heatMapHeight + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawClusterHeatMap( seuratMut, heatMapHeight, heatMapWidth))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "clusterMutHeatMap.png"),
	height = heatMapHeight + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawClusterHeatMap( seuratMut, heatMapHeight, heatMapWidth))
dev.off()}

seuratCombined			<- MergeSeurat( bestUmap, seuratMut, do.normalize = FALSE, do.scale = TRUE, names.field = NULL)
seuratCombined 			<- RunPCA( seuratCombined, pc.genes = rownames( seuratCombined@data))
seuratCombined			<- BuildSNN( seuratCombined, dims.use = 1:20)
seuratCombined@project.name	<- "allCells"

#restore cluster idents

seuratCombined@meta.data$newId <- "newID"
seuratCombined@meta.data[ names( seuratMut@ident) , "newId"] <- paste0( "mut_", as.character(seuratMut@ident))
seuratCombined@meta.data[ names( bestUmap@ident)  , "newId"] <- as.character( bestUmap@ident)
seuratCombined			<- SetAllIdent( seuratCombined, id = "newId")


seuratCombined			<- BuildClusterTree( seuratCombined, pcs.use = 1:10, do.reorder = TRUE)  


plotClusterTree( seuratCombined, plotDPI = plotDPI, treeName = "combinedClusterTree")


makeDimPlot( seuratMut, dimRed = "umap", col = setClusterColors( seuratMut), orientation = "landscape", plotDPI = plotDPI, name = "mutDimPlot")

makeDotPlot( seuratCombined, balanced = TRUE, nLines = length( levels( seuratCombined@ident)), plotDPI = plotDPI, orientation = "landscape", name = "combDotPlot")

#now validate combined clusters
#validation of mutant clusters only results in merging clusters 0 -- 2 and 1 -- 3

valCutoff 		<- 0.86


valCombined <- ValidateClusters( seuratCombined, pc.use = 1:8, top.genes = 4, min.connectivity = 0.005, acc.cutoff = 0.75)

for (i in seq(0.76, valCutoff, 0.01)){ cat(i, "\n")
   valCombined <- ValidateClusters( valCombined, pc.use = 1:8, top.genes = 4, min.connectivity = 0.005, acc.cutoff = i)
}

valCombined <- BuildClusterTree( valCombined, pcs.use = 1:8, do.reorder = TRUE, reorder.numeric = TRUE, do.plot = TRUE)

levels(valCombined@ident) <- names( getFinalClusterTypes( valCombined))
valCombined <- BuildClusterTree( valCombined, pcs.use = 1:8, do.reorder = TRUE, reorder.numeric = FALSE, do.plot = TRUE)

plotClusterTree( valCombined, plotDPI = plotDPI, treeName = "validatedCombinedClusterTree")



valCombined 	<- StashIdent( valCombined, save.name = valCutoffIdentName)

makeDotPlot( valCombined, balanced = TRUE, nLines = length( levels( valCombined@ident)), plotDPI = plotDPI, orientation = "landscape", name = "valCombDotPlot")

source("R/taqmanAnalysis.r")
