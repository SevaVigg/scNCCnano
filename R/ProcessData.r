source("R/ReadSourceFiles.r")
source("R/writeDupFile.r")

workDir <- getwd()

rawPath 	 <- file.path( workDir, "SourceData")
resDir		 <- file.path( workDir, "Res")

dir.create( resDir, showWarnings = FALSE)

initialTablesPath 	<- file.path( resDir, "InitialTables")
dir.create( initialTablesPath, showWarnings = FALSE)

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

PCADir		<- file.path( plotDir, "PCA")
dir.create( PCADir, showWarnings = FALSE)

source("R/seuratNorm.r")
seuratWT 	<- seuratNorm("WT")
seuratAll	<- seuratNorm("allCells")

seuratAll <- SetAllIdent( seuratAll, id = "genCellTypeIdent")
seuratWT <- SetAllIdent( seuratWT, id = "genCellTypeIdent")

dotPlotDir <- file.path( plotDir, "dotPlots")
dir.create( dotPlotDir, showWarnings = FALSE)


WTdotPlot 	<- DotPlot(seuratWT, genes.plot = rev(rownames(seuratWT@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))
	WTdotPlot		<- WTdotPlot +
		theme(
			legend.position="none", 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)

png( file.path( dotPlotDir, "WTdotPlot.png"), width = 1600, height = 600)
	(WTdotPlot)
dev.off()


#init cell type distributions, also prepares seuratWT and seuratAll objects
#source("R/makeGeneSpacePlots.r")  

#PCA analysis
source("R/makeInitCellTypePCAPlots.r")

png( file.path( PCADir, "WT_PCAcomps.png"), width = 1536, height = 2048)
	(makeInitCellTypePCAPlots( seuratWT, nComps = 5))
dev.off()
png( file.path( PCADir, "All_PCAcomps.png"), width = 1536, height = 2048)
	(makeInitCellTypePCAPlots( seuratAll, nComps = 5))
dev.off()

#find best parameters for UMAP clustering
source("R/findBestUmapClusters.r")
bestUmap <- findBestUmapClusters( seuratWT)


source("R/getFinalClusterTypes.r")
levels(bestUmap@ident) <- names(getFinalClusterTypes( bestUmap))

cgDotPlot 	<- DotPlot( bestUmap, genes.plot = rev(rownames( bestUmap@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))
	cgDotPlot		<- cgDotPlot +
		theme(
			legend.position="none", 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)

png( file.path( dotPlotDir, "coarseGrainDotPlot.png"), width = 1600, height = 600)
	( cgDotPlot)
dev.off()

#build coarse grain cluster tree
bestUmap <- BuildClusterTree( bestUmap, pcs.use = 1:7, do.reorder = TRUE, reorder.numeric = TRUE, do.plot = FALSE)
levels(bestUmap@ident) <- names(getFinalClusterTypes( bestUmap))
bestUmap <- BuildClusterTree( bestUmap, pcs.use = 1:7, do.reorder = FALSE, reorder.numeric = FALSE, do.plot = FALSE)

clusterDataDir		<- file.path(resDir, "clusterData")
dir.create( clusterDataDir, showWarnings = FALSE)

clusterPlotDir		<- file.path( plotDir, "clusterPlots")
dir.create( clusterPlotDir, showWarnings = FALSE)


png( file.path( clusterPlotDir, "coarseGrainClusterTree.png"))
	PlotClusterTree( bestUmap) 
dev.off()


bestUmap <- BuildClusterTree( bestUmap, pcs.use = 1:7, do.reorder = FALSE, reorder.numeric = FALSE, do.plot = FALSE)


clusterDataDir		<- file.path(resDir, "clusterData")
dir.create( clusterDataDir, showWarnings = FALSE)

clusterPlotDir		<- file.path( plotDir, "clusterPlots")
dir.create( clusterPlotDir, showWarnings = FALSE)


png( file.path( clusterPlotDir, "coarseGrainClusterTree.png"))
	PlotClusterTree( bestUmap) 
dev.off()


save( bestUmap, file = file.path( clusterDataDir, "bestUmap.rObj"))
bestUmap <- StashIdent( bestUmap, save.name = "bestClustersIdent")

source("R/make2Dmap.r")
umap2Dinit <- make2Dmap( seuratWT)

#save( file = file.path( clusterDataDir, "visualisation2Dumap.rObj"), umap2Dinit)
load( file = file.path( clusterDataDir, "visualisation2Dumap.rObj"))
load( file = file.path( clusterDataDir, "bestUmap.rObj"))

source("R/makeFeaturePlots.r")
makeFeaturePlots( umap2Dinit, minCutoff = 3, "umap")

png( file.path( clusterPlotDir, "InitCellTypeUmap.png"), width = 800, height = 600)
	InitCellTypeUmapPlot <- DimPlot( umap2Dinit, reduction.use = "umap", cols.use = setClusterColors( umap2Dinit), pt.size = 2) + 
		theme( axis.text.x = element_text( size = 20), axis.text.y = element_text(size = 20),
		       axis.title.x = element_text( size = 20, margin = margin( t = 5, r = 0, b = 0, l = 0)), axis.title.y = element_text( size = 20))
	(InitCellTypeUmapPlot)
dev.off()


source("R/setClusterColors.r")

valCutoff 		<- 0.89

valCutoffIdentName	<- paste0( "cutoffIdent", valCutoff)


valUmap <- ValidateClusters( bestUmap, pc.use = 1:7, top.genes = 4, min.connectivity = 0, acc.cutoff = 0.7)
for (i in seq(0.71, valCutoff, 0.01)){ cat(i, "\n")
  valUmap <- ValidateClusters( valUmap, pc.use = 1:7, top.genes = 4, min.connectivity = 0, acc.cutoff = i)
}

levels(valUmap@ident) <- names( getFinalClusterTypes( valUmap))
valUmap <- BuildClusterTree( valUmap, pcs.use = 1:8, do.reorder = FALSE, reorder.numeric = FALSE, do.plot = TRUE)

png( file.path( clusterPlotDir, "validatedClusterTree.png"))
	PlotClusterTree( valUmap) 
dev.off()

valUmap 	<- StashIdent( valUmap, save.name = valCutoffIdentName)

umap2Dclust <- umap2Dinit
umap2Dclust@meta.data[, valCutoffIdentName] 	<- valUmap@meta.data[, valCutoffIdentName]
umap2Dclust					<- SetAllIdent( umap2Dclust, id = valCutoffIdentName)

#we also need pseudotime trajectories
source("R/plot2DAllCurves.r")

png( file.path( clusterPlotDir, "InitClusterUmap_Common.png"), width = 800, height = 600)
			InitCellTypeUmapPlot <- DimPlot( umap2Dclust, reduction.use = "umap", cols.use = setClusterColors( umap2Dclust), pt.size = 2) + 
		theme( axis.text.x = element_text( size = 20), axis.text.y = element_text(size = 20),
		       axis.title.x = element_text( size = 20, margin = margin( t = 5, r = 0, b = 0, l = 0)), axis.title.y = element_text( size = 20))
	(InitCellTypeUmapPlot)
	plot2DAllCurves( umap2Dclust, bestUmap@ident, dimRed = "umap")
dev.off()

png( file.path( PCADir, "All_PCAcomps.png"), width = 1536, height = 2048)
	(makeInitCellTypePCAPlots( valUmap, nComps = 5))
dev.off()

clDotPlot 	<- DotPlot( valUmap, genes.plot = rev(rownames( valUmap@data)), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE, cols.use = c("cyan", "red"))
	clDotPlot		<- clDotPlot +
		theme(
			legend.position="none", 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)

png( file.path( dotPlotDir, "clDotPlot.png"), width = 1600, height = 600)
	( clDotPlot)
dev.off()

source("R/createSlingShotObject.r")
source("R/getTargetCurve.r")
source("R/drawHeatMap.r")

heatMapDir <- file.path( plotDir, "heatMaps")
dir.create( heatMapDir, showWarnings = FALSE)

slingWTcg	<- createSlingShotObject( bestUmap, "umap")
slingWT2D	<- createSlingShotObject( umap2Dclust, "umap")

png( file.path( heatMapDir, "I_heatMap.png"), width = 600, height = 800)
	(drawHeatMap( umap2Dclust, getTargetCurve( bestUmap, target = "I"), do.print = TRUE))  
dev.off()
