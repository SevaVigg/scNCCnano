#This is a minimal version of SingleCellExperiment preparing
#it requires scTables prepared by makeScTables.r
#this script also sets up one of the cell sets from "allCells", "WT", "WT_Sox10" )

if (!require("SC3")){
BiocManager::install("SC3")
library(SC3)}

if (!require("SingleCellExperiment"))
BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)

resDir		<- file.path(getwd(), "Res")
scTablesDir	<- file.path( resDir, "scTables")

experimentType	<- "allCells"	# may be ("allCells", "WT", and "WT_Sox10")	

dataDir		<- file.path( scTablesDir, experimentType)	#this is the data dir showing which source file is to use
logExps		<- read.table( file = file.path( dataDir, "logExpTableDedupQCimp.csv"  ), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
geneExpsImp	<- read.table( file = file.path( dataDir, "geneExpTableDedupQCimp.csv"  ), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
Cells		<- read.table( file = file.path( dataDir, "cellDescriptionsDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )

#reorder cells
cells_ind 	<- order(as.numeric(Cells["hpf",]))			# order with hpf increasing
logExps		<- logExps[, cells_ind]
geneExpsImp	<- geneExpsImp[, cells_ind]
Cells		<- Cells[, cells_ind]

#rename cell types, prepare the annotated cell table
celltype 	<- unlist(lapply( Cells, function(x) if (x[6] == "regular") return( x[3]) else return( x[6] ))) 

#seurat scales the data by mean and sd. Let us scale the data with 

ipmc_sce    	<- SingleCellExperiment(
    			assays = list(
        			counts = as.matrix(geneExpsImp),
        			logcounts = as.matrix(logExps)
    				)) 
rowData(ipmc_sce)$feature_symbol <- rownames(ipmc_sce)

ipmc_sce 	<- sc3(ipmc_sce, ks = 5:20, gene_filter = FALSE, biology = TRUE)

