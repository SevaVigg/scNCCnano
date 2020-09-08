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

seuratWT 	<- seuratNorm("WT")
seuratAll	<- seuratNorm("allCells")

seuratAll <- SetAllIdent( seuratAll, id = "genCellTypeIdent")
seuratWT <- SetAllIdent( seuratWT, id = "genCellTypeIdent")

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

clusterDataDir		<- file.path(resDir, "clusterData")
dir.create( clusterDataDir, showWarnings = FALSE)

save( bestUmap, file = file.path( clusterDataDir, "bestUmap.rObj"))
bestUmap <- StashIdent( bestUmap, save.name = "bestClustersIdent")

#source("R/make2Dmap.r")
#umap2D <- make2Dmap( seuratWT)

load( file = file.path( clusterDataDir, "visualisation2Dumap"))
load( file = file.path( clusterDataDir, "bestUmap.rObj"))

source("R/makeFeaturePlots.r")
makeFeaturePlots( umap2D, minCutoff = 3, "umap")

valCutoff 		<- 0.85

valCutoffIdentName	<- paste0( "cutoffIdent", valCutoff)
umapVal <- ValidateClusters( bestUmap, pc.use = 1:6, top.genes = 3, min.connectivity = 0, acc.cutoff = valCutoff)
umapVal <- BuildClusterTree( umapVal, pcs.use = 1:6, do.reorder = TRUE, reorder.numeric = TRUE)

source("R/getClusterTypes.r")
levels(umapVal@ident) <- names( getClusterTypes( umapVal))

umapVal 				<- StashIdent( umapVal, save.name = valCutoffIdentName)
umap2D@meta.data[,valCutoffIdentName] 	<- umapVal@meta.data[, valCutoffIdentName]
umap2D					<- SetAllIdent( umap2D, id = valCutoffIdentName)
umapVal <- BuildClusterTree( umapVal, pcs.use = 1:6)



