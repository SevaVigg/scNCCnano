#This is a minimal version of seurat preparing
#it requires scTables prepared by makeScTables.r
#this script also sets up one of the cell sets from "AllCells", "WT", "WT_Sox10" )

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
cells_ind 	<- order(as.numeric(Cells["hpf",]))			# order with hpf increasing
logExps		<- logExps[, cells_ind]
geneExpsImp	<- geneExpsImp[, cells_ind]
Cells		<- Cells[, cells_ind]

# 
require(gtools)
dens		<- density( t(as.matrix(logExps)))
expThreshold	<- optimize(approxfun(dens$x,dens$y),interval=c(5,14))$minimum 

#logExps		<- apply( logExps, c(1, 2), function(x) if (x < 5.) 0 else x) 
#logExps		<- apply( logExps, c(1, 2), function(x) if (x > 15) 15 else x) 

#rename cell types, prepare the annotated cell table
celltype 	<- unlist(lapply( Cells, function(x) if (x[6] == "general") return( x[3]) else return( x[6] ))) 

#seurat scales the data by mean and sd. Let us scale the data with 

ipmc    	<- CreateSeuratObject( raw.data = as.matrix(logExps))
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

genCellTypeIdentDF	<- data.frame( genCellTypeIdent =sapply( ipmc@meta.data$originalCellTypes, function(x) if(grepl("^[0-9]", x)) return("G") else return(x)), row.names = colnames( Cells)) 
ipmc		<- AddMetaData( object = ipmc, genCellTypeIdentDF, col.name = "genCellTypeIdent")
ipmc@misc	<- "allGenes"




