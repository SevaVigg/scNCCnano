#This is a minimal version of seurat preparing

require(Seurat)
require(methods)

experimentType	<- "allCells"	# may be ("allCells", "WT", and "WT_Sox10")	 

resDir		<- file.path(getwd(), "Res")
scTablesDir	<- file.path( resDir, "scTables")

dataDir		<- file.path( scTablesDir, experimentType)	#this is the data dir showing which source file is to use
logExps	<- read.table( file = file.path( dataDir, "logExpTableDedupQCimp.csv"  ), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
Cells	<- read.table( file = file.path( dataDir, "cellDescriptionsDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )

#reorder cells
cells_ind 	<- order(as.numeric(Cells["hpf",]))			# order with hpf increasing
logExps		<- logExps[, cells_ind]
Cells		<- Cells[, cells_ind]

#rename cell types, prepare the annotated cell table
types 		<- unique(paste0(Cells["hpf",], "_", Cells["CellType",]))
hpf_CellType	<- t(data.frame(hpf_CellType = paste0(Cells["hpf",], "_", Cells["CellType",]), row.names = colnames(Cells)))
Cells		<- rbind(Cells, hpf_CellType)

newTypes   	<- c("18", "21", "24", "Tl", "30", "mitfa-", "sox10-", "36", "48", "I", "M", "60", "72")
names(newTypes)	<- types

#seurat scales the data by mean and sd. Let us scale the data with 

ipmc    	<- CreateSeuratObject( raw.data = as.matrix(logExps) )
ipmc@project.name <- experimentType
ipmc    	<- AddMetaData( object = ipmc, t(Cells), col.name = rownames(Cells) )

newTypeDF	<- data.frame( newType = character(ncol(Cells)), row.names = colnames(Cells) )
cellNamesDF	<- data.frame( cellNames = colnames(Cells), row.names = colnames(Cells))

ipmc		<- AddMetaData( object = ipmc, newTypeDF, col.name = "newType")
ipmc		<- AddMetaData( object = ipmc, cellNamesDF, col.name = "cellNames")

levels(ipmc@ident)    	<- newTypes

ipmc@ident		<- as.factor(unlist( lapply( ipmc@meta.data[ , "hpf_CellType"], function(cell) newTypes[as.character(cell)]) ))
names(ipmc@ident) 	<- names(ipmc@ident) <- colnames(ipmc@data)
ipmc@meta.data$newType 	<- ipmc@ident

ipmc			<- ScaleData(ipmc, do.scale = TRUE, do.center = TRUE)
ipmc			<- StashIdent(object = ipmc, save.name = "originalCellTypes")




