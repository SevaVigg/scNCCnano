source("R/ReadSourceFiles.r")
source("R/writeDupFile.r")

workDir <- getwd()

rawPath 	 <- file.path( workDir, "SourceData")
resPath		 <- file.path( workDir, "Res")

dir.create( resPath, showWarnings = FALSE)

initialTablesPath 	<- file.path( resPath, "initialTables")
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

#Write deduplicated results

write.csv( CellTable$Genes,  file = file.path( initialTablesPath, "expressionTableDedup.csv"   ) )
write.csv( CellTable$Cells,  file = file.path( initialTablesPath, "cellDescripitonsDedup.csv"  ) )
write.csv( CellTable$Probes, file = file.path( initialTablesPath, "ProbesDescripitonsDedup.csv") )



