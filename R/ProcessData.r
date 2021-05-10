require( tidyr )
require( ape)
require( ComplexHeatmap)

source("R/ReadSourceFiles.r")
source("R/writeDupFile.r")

plotDPI		<- 100

workDir <- getwd()

rawPath 	 <- file.path( workDir, "SourceData")
resDir		 <- file.path( workDir, "Res")

dir.create( resDir, showWarnings = FALSE)

initialTablesPath 	<- file.path( resDir, "InitialTables")
dir.create( initialTablesPath, showWarnings = FALSE)

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
source("R/qualityControl.r")

#imputation
source("R/makeScTables.r")


#start seurat analysis
plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

source("R/seuratNorm.r")
seuratWT 	<- seuratNorm("WT")
seuratAll	<- seuratNorm("allCells")

seuratAll <- SetAllIdent( seuratAll, id = "genCellTypeIdent")
seuratWT <- SetAllIdent( seuratWT, id = "genCellTypeIdent")

#gene heatmaps


heatMapDir <- file.path( plotDir, "heatMaps")
dir.create( heatMapDir, showWarnings = FALSE)

heatMapDir100dpi	<- file.path( heatMapDir, "100dpi")
dir.create( heatMapDir100dpi, showWarnings = FALSE)

heatMapDir600dpi	<- file.path( heatMapDir, "600dpi")
dir.create( heatMapDir600dpi, showWarnings = FALSE)

#gene covariance
source( "R/drawGeneCovarHeatMap.r")

heatMapHeigth 	<- 10
heatMapWidth 	<- 10
Margin		<- 2

if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "geneCovarHeatMap.png"),
	height = heatMapHeigth + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawGeneCovarHeatMap( seuratWT, heatMapHeigth, heatMapWidth))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "geneCovarHeatMap.png"),
	height = heatMapHeigth + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawGeneCovarHeatMap( seuratWT, heatMapHeigth, heatMapWidth))
dev.off()}

#gene biclusters
source( "R/drawBiclustHeatMap.r")

heatMapHeigth 	<- 7
heatMapWidth 	<- 10
Margin		<- 2

if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "biclustWTHeatMap.png"),
	height = heatMapHeigth + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawBiclustHeatMap( seuratWT, heatMapHeigth, heatMapWidth))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "biclustWTHeatMap.png"),
	height = heatMapHeigth + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawBiclustHeatMap( seuratWT, heatMapHeigth, heatMapWidth))
dev.off()}


dotPlotDir <- file.path( plotDir, "dotPlots")
dir.create( dotPlotDir, showWarnings = FALSE)

source("R/makeDotPlot.r")

nClust = length( levels( seuratWT@ident))
makeDotPlot( seuratWT, balanced = TRUE, nLines = length( levels( seuratWT@ident)), plotDPI = plotDPI, orientation = "landscape", name = "WTdotPlot")


#init cell type distributions, also prepares seuratWT and seuratAll objects
#source("R/makeGeneSpacePlots.r")  

#PCA analysis
source("R/makeInitCellTypePCAPlots.r")

PCAPlotDir	<- file.path( plotDir, "PCAPlots")
dir.create( PCAPlotDir, showWarnings = FALSE)

makeInitCellTypePCAPlots( seuratWT,  nComps = 5, plotDPI = 100, name = "WT_PCAcomps")
makeInitCellTypePCAPlots( seuratAll, nComps = 5, plotDPI = 100, name = "All_PCAcomps")

#find best parameters for UMAP clustering
source("R/findBestUmapClusters.r")
Spread		<- 10
bestUmap 	<- findBestUmapClusters( seuratWT, Spread)
bestDim		<- dim(bestUmap@dr$umap@cell.embeddings)[2]


source("R/getFinalClusterTypes.r")
levels(bestUmap@ident) <- names(getFinalClusterTypes( bestUmap))

#build coarse grain cluster tree
bestUmap <- BuildClusterTree( bestUmap, pcs.use = 1:7, do.reorder = TRUE, reorder.numeric = TRUE, do.plot = FALSE)

levels(bestUmap@ident) <- names(getFinalClusterTypes( bestUmap))

bestUmap <- BuildClusterTree( bestUmap, pcs.use = 1:7, do.reorder = FALSE, reorder.numeric = FALSE, do.plot = FALSE)


source("R/createSlingShotObject.r")
source("R/cosineClusterDist.r")

allGenes <- rownames( bestUmap@data)
sling <- createSlingShotObject( bestUmap, genes.use = allGenes, endClust = c("I","M","X"))

clusterPlotDir		<- file.path( plotDir, "clusterPlots")
dir.create( clusterPlotDir, showWarnings = FALSE)

clusterTreeDir		<- file.path( plotDir, "clusterTrees")
dir.create( clusterTreeDir, showWarnings = FALSE)

clusterTreeDir100dpi	<- file.path( clusterTreeDir, "100dpi")
dir.create( clusterTreeDir100dpi, showWarnings = FALSE)

clusterTreeDir600dpi	<- file.path( clusterTreeDir, "600dpi")
dir.create( clusterTreeDir600dpi, showWarnings = FALSE)

source("R/plotClusterTree.r")

plotClusterTree( bestUmap, plotDPI = plotDPI, name = "coarseGrainClusterTree")

makeDotPlot( bestUmap, balanced = TRUE, nLines = length( levels( bestUmap@ident)), plotDPI = plotDPI, orientation = "landscape", name = "coarseGrainDotPlot")

clusterDataDir		<- file.path(resDir, "clusterData")
dir.create( clusterDataDir, showWarnings = FALSE)

save( bestUmap, file = file.path( clusterDataDir, "bestUmap.rObj"))
bestUmap <- StashIdent( bestUmap, save.name = "bestClustersIdent")

source("R/make2Dmap.r")
umap2Dinit <- make2Dmap( seuratWT)

source("R/setClusterColors.r")

#save( file = file.path( clusterDataDir, "visualisation2Dumap.rObj"), umap2Dinit)
load( file = file.path( clusterDataDir, "visualisation2Dumap.rObj"))
load( file = file.path( clusterDataDir, "bestUmap.rObj"))

source("R/makeDimPlot.r")
makeDimPlot( umap2Dinit, dimRed = "umap", col = setClusterColors( umap2Dinit), orientation = "landscape", plotDPI = plotDPI)


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

#cluster HeatMap

source("R/drawClusterHeatMap.r")
if (plotDPI == 600) {
png( file = file.path( heatMapDir600dpi, "clusterHeatMap.png"),
	height = heatMapHeigth + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawClusterHeatMap( valUmap, heatMapHeigth, heatMapWidth))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( heatMapDir100dpi, "clusterHeatMap.png"),
	height = heatMapHeigth + Margin,
	width =  heatMapWidth + 2*Margin,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawClusterHeatMap( valUmap, heatMapHeigth, heatMapWidth))
dev.off()}


umap2Dclust <- umap2Dinit
umap2Dclust@meta.data[, valCutoffIdentName] 	<- valUmap@meta.data[, valCutoffIdentName]
umap2Dclust					<- SetAllIdent( umap2Dclust, id = valCutoffIdentName)

#Now we have clustering and can make feature plots accompanied with clusters
source("R/makeFeaturePlots.r")
makeFeaturePlots( umap2Dclust, minCutoff = 3, "umap", plotDPI = plotDPI, name = "initFeaturePlots", orientation = "portrait")


#we also need pseudotime trajectories
source("R/plot2DAllCurves_ggplot.r") #this file contains function "plot2DAllCurves" which uses ggplot rather than lines
source("R/cosineClusterDist.r")

plot2DAllCurves( umap2Dclust, valUmap, dims = 1:8, genes.use = allGenes, plotDPI = plotDPI, name = "slingClustersUmap")

source("R/calcTSNEGeneSpace.r")
tsne2Dinit <- calcTSNEGeneSpace( seuratWT, perplexity, )

plot2DAllCurves( tsne2Dinit, valUmap, genes.use = allGenes, dimRed2D = "tsne", plotDPI = plotDPI, name = "slingClustersTSNE")


makeInitCellTypePCAPlots( valUmap, nComps = 5)

clDotPlot 	<- dotPlotBalanced( valUmap, genes.plot = rev(rownames( valUmap@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))
	clDotPlot		<- clDotPlot +
		theme(
			legend.key.size = unit( 1.2, "cm"),			
			legend.position="right",
			legend.title = element_text( size = 24),
			legend.text = element_text( size = 18), 
 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)

png( file.path( dotPlotDir, "clDotPlot.png"), width = 1600, height = 600)
	( clDotPlot)
dev.off()

#Plot special dotPlots for control cell types

Mcells 				<- WhichCells( valUmap, ident = "M")
seuratM				<- SubsetData( valUmap, ident.use = "M") 
cntrlM				<- grep( "M", Mcells, value = TRUE)
regM				<- setdiff( Mcells, cntrlM) 					
melIndexLine			<- c( rep( "contM", length( cntrlM)), rep( "regM", length( regM)))
names( melIndexLine)		<- c( cntrlM, regM)
seuratM@meta.data$Mcells 	<- melIndexLine[ rownames( seuratM@meta.data)]
seuratM <- SetAllIdent( seuratM, id = "Mcells")

mDotPlot 	<- dotPlotBalanced( seuratM, genes.plot = rev(rownames( seuratM@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))
	mDotPlot		<- mDotPlot +
		theme(
			legend.key.size = unit( 1.2, "cm"),			
			legend.position="right",
			legend.title = element_text( size = 24),
			legend.text = element_text( size = 18), 
 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)


png( file.path( dotPlotDir, "MelComp.png"), width = 1600, height = 600)
	( mDotPlot)
dev.off()

#melanocyte subtypes
seuratMSubTypes <- FindClusters( seuratM, dims.use = 1:8, resolution = 1.0)

mstDotPlot 	<- dotPlotBalanced( seuratMSubTypes, genes.plot = rev(rownames( seuratM@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))
	mstDotPlot		<- mstDotPlot +
		theme(
			legend.key.size = unit( 1.2, "cm"),			
			legend.position="right",
			legend.title = element_text( size = 24),
			legend.text = element_text( size = 18), 
 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)


png( file.path( dotPlotDir, "MelSubTypesComp.png"), width = 1600, height = 600)
	( mstDotPlot)
dev.off()



#Plot special dotPlots for iridophores cell types

Icells 				<- WhichCells( valUmap, ident = "I")
seuratI				<- SubsetData( valUmap, ident.use = "I") 
cntrlI				<- grep( "I", Icells, value = TRUE)
regI				<- setdiff( Icells, cntrlI) 					
irdIndexLine			<- c( rep( "contI", length( cntrlI)), rep( "regI", length( regI)))
names( irdIndexLine)		<- c( cntrlI, regI)
seuratI@meta.data$Icells 	<- irdIndexLine[ rownames( seuratI@meta.data)]
seuratI <- SetAllIdent( seuratI, id = "Icells")

iDotPlot 	<- dotPlotBalanced( seuratI, genes.plot = rev(rownames( seuratI@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))
	iDotPlot		<- iDotPlot +
		theme(
			legend.key.size = unit( 1.2, "cm"),			
			legend.position="right",
			legend.title = element_text( size = 24),
			legend.text = element_text( size = 18), 
 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)


png( file.path( dotPlotDir, "IrdComp.png"), width = 1600, height = 600)
	( iDotPlot)
dev.off()


mlGeneSort <- rownames(AverageExpression(seuratM)[order(-AverageExpression(seuratM)$contM),])

png( file.path( heatMapDir, "MlHeatMap.png"), width = 1000, height = 1600)
	( DoHeatmap( seuratM, genes.use = mlGeneSort, use.scaled = FALSE) + theme(
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"))
)
dev.off()

#For iridophores

Icells 				<- WhichCells( valUmap, ident = "I")
seuratI				<- SubsetData( valUmap, ident.use = "I") 
cntrlI				<- grep( "I", Icells, value = TRUE)
regI				<- setdiff( Icells, cntrlI) 					
iphIndexLine			<- c( rep( "contI", length( cntrlI)), rep( "regI", length( regI)))
names( iphIndexLine)		<- c( cntrlI, regI)
seuratI@meta.data$Icells 	<- iphIndexLine[ rownames( seuratI@meta.data)]
seuratI <- SetAllIdent( seuratI, id = "Icells")

png( file.path( dotPlotDir, "IphComp.png"), width = 1600, height = 600)
	 barplot(as.matrix(t(AverageExpression( seuratI))), beside = TRUE, las = 2, legend.text = TRUE)
dev.off()

IpDotPlot 	<- dotPlotBalanced( seuratI, genes.plot = rev(rownames( seuratI@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))
	clDotPlot		<- clDotPlot +
		theme(
			legend.key.size = unit( 1.2, "cm"),			
			legend.position="right",
			legend.title = element_text( size = 24),
			legend.text = element_text( size = 18), 
 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)

png( file.path( dotPlotDir, "IphDotPlot.png"), width = 1600, height = 600)
	( IpDotPlot)
dev.off()

ipGeneSort <- rownames(AverageExpression(seuratI)[ order(-AverageExpression(seuratI)$contI),])

png( file.path( heatMapDir, "IpHeatMap.png"), width = 1200, height = 1600)

	( DoHeatmap( seuratI, genes.use = ipGeneSort, use.scaled = FALSE) + theme(
			legend.key.size = unit( 1.2, "cm"),
			legend.position = "right",
			legend.title = element_text( size = 24),
			legend.text = element_text( size = 18), 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"))
)
dev.off()


source("R/getTargetCurve.r")
source("R/drawHeatMap.r")


png( file.path( heatMapDir, "I_heatMap.png"), width = 600, height = 800)
	(drawHeatMap( umap2Dclust, getTargetCurve( valUmap, target = "I", genes.use = allGenes), do.print = TRUE))  
dev.off()

png( file.path( heatMapDir, "M_heatMap.png"), width = 600, height = 800)
	(drawHeatMap( umap2Dclust, getTargetCurve( valUmap, target = "M", genes.use = allGenes), do.print = TRUE))  
dev.off()

png( file.path( heatMapDir, "X_heatMap.png"), width = 600, height = 800)
	(drawHeatMap( umap2Dclust, getTargetCurve( bestUmap, target = "X", genes.use = allGenes), do.print = TRUE))  
dev.off()

#vlnPlots

vlnPlotDir <- file.path( plotDir, "vlnPlots")
dir.create( vlnPlotDir, showWarnings = FALSE)

source("R/makeNanoStringVlnPlots.r")
makeNanoStringVlnPlots( valUmap) 


#Now consider mutant cells


seuratMut 	<- SubsetData( seuratAll, ident.use = "sox10-")
	

#First make dotplot with mutants

allDotPlot 	<- dotPlotBalanced(seuratAll, genes.plot = rev(rownames(seuratAll@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))
	allDotPlot		<- allDotPlot +
		theme(
			legend.key.size = unit( 1.2, "cm"),			
			legend.position="right",
			legend.title = element_text( size = 24),
			legend.text = element_text( size = 18), 
 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)

png( file.path( dotPlotDir, "allDotPlot.png"), width = 1600, height = 600)
	(allDotPlot)
dev.off()


paramStr	<- tail( colnames( bestUmap@meta.data), 1)

bestDim		<- 2
bestMinDist	<- 3.5
bestR		<- 0.7

seuratMut		<- calcUmapGeneSpace( seuratMut, Dim = bestDim, myNeighbors = 25L, mySpread = Spread,  
				minDist = bestMinDist,  UMAPRandSeed = 42L, experimentType <- "allCells")$All

seuratMut 		<- FindClusters( seuratMut, reduction.type = 'umap', 
					dims.use = 1:bestDim, 
					k.param = 5, print.output = FALSE, force.recalc = TRUE, save.SNN = TRUE,
					resolution = bestR) 




mutDotPlot 	<- dotPlotBalanced(seuratMut, genes.plot = rev(rownames(seuratMut@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))
	mutDotPlot		<- mutDotPlot +
		theme(
			legend.key.size = unit( 1.2, "cm"),			
			legend.position="right",
			legend.title = element_text( size = 24),
			legend.text = element_text( size = 18), 
 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)


png( file.path( dotPlotDir, "mutDotPlot.png"), width = 1600, height = 600)
	(mutDotPlot)
dev.off()

seuratCombined		<- MergeSeurat( bestUmap, seuratMut, do.normalize = FALSE, do.scale = TRUE, names.field = NULL)
seuratCombined 		<- RunPCA( seuratCombined, pc.genes = rownames( seuratCombined@data))
seuratCombined		<- BuildSNN( seuratCombined, dims.use = 1:20)

#restore cluster idents

seuratCombined@meta.data$newId <- "newID"
seuratCombined@meta.data[ names( seuratMut@ident) , "newId"] <- paste0( "mut_", as.character(seuratMut@ident))
seuratCombined@meta.data[ names( bestUmap@ident)  , "newId"] <- as.character( bestUmap@ident)
seuratCombined			<- SetAllIdent( seuratCombined, id = "newId")


seuratCombined			<- BuildClusterTree( seuratCombined, pcs.use = 1:10, do.reorder = TRUE)  
png( file.path( clusterTreeDir, "combinedClusterTree.png"))
	PlotClusterTree( seuratCombined) 
dev.off()


png( file.path( clusterPlotDir, "mutClusterUmap.png"), width = 800, height = 600)
	mutClustUmapPlot <- DimPlot( seuratMut, reduction.use = "umap", cols.use = setClusterColors( seuratMut), pt.size = 2) + 
		theme( axis.text.x = element_text( size = 20), axis.text.y = element_text(size = 20),
		       axis.title.x = element_text( size = 20, margin = margin( t = 5, r = 0, b = 0, l = 0)), axis.title.y = element_text( size = 20))
	(mutClustUmapPlot)
dev.off()

combDotPlot 	<- dotPlotBalanced(seuratCombined, genes.plot = rev(rownames(seuratCombined@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))
	combDotPlot		<- combDotPlot +
		theme(
			legend.key.size = unit( 1.2, "cm"),			
			legend.position="right",
			legend.title = element_text( size = 24),
			legend.text = element_text( size = 18), 
 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)

png( file.path( dotPlotDir, "combDotPlot.png"), width = 1600, height = 600)
	(combDotPlot)
dev.off()

valCutoff 		<- 0.92

valCutoffIdentName	<- paste0( "cutoffIdent", valCutoff)

valCombined <- ValidateClusters( seuratCombined, pc.use = 1:8, top.genes = 4, min.connectivity = 0.005, acc.cutoff = 0.75)

for (i in seq(0.76, valCutoff, 0.01)){ cat(i, "\n")
   valCombined <- ValidateClusters( valCombined, pc.use = 1:8, top.genes = 4, min.connectivity = 0.005, acc.cutoff = i)
}

valCombined <- BuildClusterTree( valCombined, pcs.use = 1:8, do.reorder = TRUE, reorder.numeric = TRUE, do.plot = TRUE)
levels(valCombined@ident) <- names( getFinalClusterTypes( valCombined))
valCombined <- BuildClusterTree( valCombined, pcs.use = 1:8, do.reorder = TRUE, reorder.numeric = FALSE, do.plot = TRUE)

png( file.path( clusterTreeDir, "validatedCombinedClusterTree.png"))
	PlotClusterTree( valCombined) 
dev.off()

valCombined 	<- StashIdent( valCombined, save.name = valCutoffIdentName)

valCombDotPlot 	<- dotPlotBalanced( valCombined, genes.plot = rev(rownames( valCombined@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))
	valCombDotPlot		<- valCombDotPlot +
		theme(
			legend.key.size = unit( 1.2, "cm"),			
			legend.position="right",
			legend.title = element_text( size = 24),
			legend.text = element_text( size = 18), 
 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)

png( file.path( dotPlotDir, "valCombDotPlot.png"), width = 1600, height = 600)
	(valCombDotPlot)
dev.off()


