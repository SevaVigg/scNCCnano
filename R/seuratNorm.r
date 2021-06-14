seuratNorm 	<- function( experimentType){

# This snippet prepares a minimal version of seurat object preparing
# it requires scTables prepared by makeScTables.r
# experimentType can be "AllCells", "WT", "WT_Sox10", "Taqman" )
# written by Vsevolod J. Makeev 2017 - 2021

if (!require("Seurat")){
BiocManager::install("Seurat")
library(Seurat)}

if (!require("methods")){
BiocManager::install("methods")
library(methods)}

if(!require("e1071")){
  install.packages("e1071")
}


#experimentType	<- "AllCells"	# may be ("AllCells", "WT", and "taqman")	 

resDir		<- file.path(getwd(), "Res")
scTablesDir	<- file.path( resDir, "scTables")

dataDir		<- file.path( scTablesDir, experimentType)	#this is the data dir showing which source file is to use

if( experimentType == "taqman"){ 
	logExps		<- t(read.table( file = file.path( dataDir, "taqmanLogExps.csv"  ), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE ))
	Cells		<- as.data.frame(colnames( logExps))
	rownames(Cells)	<- colnames( logExps)
	celltype	<- rep( "regular", nrow(Cells))
	names(celltype) <- colnames( logExps)	
	Cells$CellType	<- celltype
	Cells		<- t(Cells)
	isExpressionThreshold <- 1
}else{
	logExps		<- read.table( file = file.path( dataDir, "logExpTableDedupQCimp.csv"  ), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
	geneExpsImp	<- read.table( file = file.path( dataDir, "geneExpTableDedupQCimp.csv"  ), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
	Cells		<- read.table( file = file.path( dataDir, "cellDescriptionsDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
	celltype 	<- unlist(lapply( Cells, function(x) if (x[6] == "regular") return( x[3]) else return( x[6] ))) 
	isExpressionThreshold <- log10(5)
}

#seurat scales the data by mean and sd. Let us scale the data with 

logExps		<- apply( logExps, 1:2, function(x) ifelse( x >= isExpressionThreshold, x, 0))
keep		<- which( apply( logExps, 2, sum) > 0 )
logExps		<- logExps[ , keep]
Cells		<- Cells[ , keep]
celltype	<- celltype[ keep]  

ipmc    	<- CreateSeuratObject( raw.data = as.matrix(logExps), is.expr = isExpressionThreshold)
ipmc@project.name <- experimentType
ipmc    	<- AddMetaData( object = ipmc, t(Cells), col.name = rownames(Cells) )

newTypeDF	<- data.frame( newType = character(ncol(Cells)), row.names = colnames(Cells) )
cellNamesDF	<- data.frame( cellNames = colnames(Cells), row.names = colnames(Cells))

ipmc		<- AddMetaData( object = ipmc, newTypeDF, col.name = "newType")
ipmc		<- AddMetaData( object = ipmc, cellNamesDF, col.name = "cellNames")

ipmc@ident	<- as.factor( celltype)

ipmc		<- ScaleData(ipmc, do.scale = TRUE, do.center = TRUE)

ipmc 		<- RunPCA( ipmc, pc.genes = rownames(ipmc@data), pcs.compute = 10, weight.by.var = FALSE, do.print = FALSE)

ipmc		<- StashIdent(object = ipmc, save.name = "originalCellTypes")

genCellTypeIdentDF	<- data.frame( genCellTypeIdent =sapply( ipmc@meta.data$originalCellTypes, function(x) if(grepl("^[0-9]", x)) return("R") else return(x)), row.names = colnames( Cells)) 
ipmc		<- AddMetaData( object = ipmc, genCellTypeIdentDF, col.name = "genCellTypeIdent")

#a legacy note

ipmc@misc	<- "allGenes"
return( ipmc)
}



