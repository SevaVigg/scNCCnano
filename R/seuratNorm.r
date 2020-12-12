seuratNorm 	<- function( experimentType){

#This is a minimal version of seurat preparing
#it requires scTables prepared by makeScTables.r
#experimentType can be "AllCells", "WT", "WT_Sox10" )

if (!require("Seurat")){
BiocManager::install("Seurat")
library(Seurat)}

if (!require("methods")){
BiocManager::install("methods")
library(methods)}

if(!require("e1071")){
  install.packages("e1071")
}


#experimentType	<- "AllCells"	# may be ("AllCells", "WT", and "WT_Sox10")	 

resDir		<- file.path(getwd(), "Res")
scTablesDir	<- file.path( resDir, "scTables")

dataDir		<- file.path( scTablesDir, experimentType)	#this is the data dir showing which source file is to use
logExps		<- read.table( file = file.path( dataDir, "logExpTableDedupQCimp.csv"  ), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
geneExpsImp	<- read.table( file = file.path( dataDir, "geneExpTableDedupQCimp.csv"  ), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
Cells		<- read.table( file = file.path( dataDir, "cellDescriptionsDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )

#reorder cells
#cells_ind 	<- order(as.numeric(Cells["CellType",]))			# order with CellType increasing
#logExps		<- logExps[, cells_ind]
#geneExpsImp	<- geneExpsImp[, cells_ind]
#Cells		<- Cells[, cells_ind]

celltype 	<- unlist(lapply( Cells, function(x) if (x[6] == "general") return( x[3]) else return( x[6] ))) 

#seurat scales the data by mean and sd. Let us scale the data with 

ipmc    	<- CreateSeuratObject( raw.data = as.matrix(logExps), is.expr = log10(5))
ipmc@project.name <- experimentType
ipmc    	<- AddMetaData( object = ipmc, t(Cells), col.name = rownames(Cells) )

newTypeDF	<- data.frame( newType = character(ncol(Cells)), row.names = colnames(Cells) )
cellNamesDF	<- data.frame( cellNames = colnames(Cells), row.names = colnames(Cells))

ipmc		<- AddMetaData( object = ipmc, newTypeDF, col.name = "newType")
ipmc		<- AddMetaData( object = ipmc, cellNamesDF, col.name = "cellNames")

ipmc@ident	<- as.factor( celltype)

ipmc		<- ScaleData(ipmc, do.scale = TRUE, do.center = TRUE)

ipmc 		<- RunPCA( ipmc, pc.genes = rownames(ipmc@data), weight.by.var = FALSE, do.print = FALSE)

ipmc		<- StashIdent(object = ipmc, save.name = "originalCellTypes")

genCellTypeIdentDF	<- data.frame( genCellTypeIdent =sapply( ipmc@meta.data$originalCellTypes, function(x) if(grepl("^[0-9]", x)) return("R") else return(x)), row.names = colnames( Cells)) 
ipmc		<- AddMetaData( object = ipmc, genCellTypeIdentDF, col.name = "genCellTypeIdent")

#now reorder levels so that I and M go togeather in the plot

ipmc@misc	<- "allGenes"
return( ipmc)
}



