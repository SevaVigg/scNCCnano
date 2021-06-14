# This script prepares files, which can be read from Seurat
# It also runs the procedure that performes imputation 
# Preparing clean  data for WT and WT+sox10 mutants, which it places into different folders, so seuratNorm.r can read the correct source files
# This script is the next in the pipeline after the qualityControl.r
#
# written by Vsevolod J. Makeev 2017 - 2021

source("R/imputeDropouts.r")

resDir		<- file.path( getwd(), "Res")
scTablesDir	<- file.path( resDir, "scTables")
dir.create( scTablesDir, showWarnings = FALSE)

dirQCres <- file.path( resDir, "QualityControl", "ResQC") 

genesQC	<- read.table( file = file.path( dirQCres, "New_NormalizedExTable.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
cellsQC	<- read.table( file = file.path( dirQCres, "New_cellDescripitonsDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )


exps		<- genesQC
allGenes	<- rownames(genesQC)

allCellsDir 	<- file.path( scTablesDir, "allCells")
dir.create( allCellsDir , showWarnings = FALSE)

#remove values, that are too close to zero
noiseTol	<- 30
exps		<- apply( exps, c(1,2), function(x) if(x>noiseTol) x else 0)

#remove cells with all zero values
exps 		<- exps[, which( apply( exps, 2, sum) >0)]
logExps 	<- log10(1+exps)

#impute for dropouts
randSeedAll 	<- as.integer(Sys.time())
logExpsImp	<- imputeDropouts(logExps, randSeedAll)

logNoiseTol	<- log10( noiseTol)

#remove imputed values, that are too close to zero
logExpsImp	<- apply( logExpsImp, c(1,2), function(x) if(x>logNoiseTol) x else 0)

#remove cells with all zero values
logExpsImp	<- logExpsImp[, which( apply( logExpsImp, 2, sum) >0)]
cellsQCImp	<- cellsQC[, colnames(logExpsImp)]
genesQCImp	<- 10^logExpsImp - 1

cat( file = file.path( allCellsDir, "ImputeSeedAll.txt"), randSeedAll)
write.table( logExps, file = file.path( allCellsDir, "logExpTableDedupQC.csv"), sep = "\t")
write.table( logExpsImp, file = file.path( allCellsDir, "logExpTableDedupQCimp.csv"), sep = "\t")
write.table( cellsQCImp, file = file.path( allCellsDir, "cellDescriptionsDedupQC.csv"), sep = "\t")
write.table( genesQCImp, file = file.path( allCellsDir, "geneExpTableDedupQCimp.csv"), sep = "\t")

WTDir		<- file.path( scTablesDir, "WT")
dir.create(WTDir, showWarnings = FALSE)	

randSeedWT 	<- as.integer(Sys.time())
logExpsWT	<- logExps[, setdiff( seq_along(colnames(logExps)), grep("(sox10|mitfa)", colnames(logExps)))]
logExpsWTImp	<- imputeDropouts(logExpsWT, randSeedWT)

logExpsWTImp	<- apply( logExpsWTImp, c(1,2), function(x) if(x>logNoiseTol) x else 0)
logExpsWTImp	<- logExpsWTImp[, which( apply( logExpsWTImp, 2, sum) >0)]
cellsQCWTImp	<- cellsQC[, colnames(logExpsWTImp)]
genesQCWTImp	<- 10^logExpsWTImp - 1

cat( file = file.path( WTDir, "ImputeSeedWT.txt"), randSeedWT)
write.table( logExpsWT, file = file.path( WTDir, "logExpTableDedupQC.csv"), sep = "\t")
write.table( logExpsWTImp, file = file.path( WTDir, "logExpTableDedupQCimp.csv"), sep = "\t")
write.table( cellsQCWTImp, file = file.path( WTDir, "cellDescriptionsDedupQC.csv"), sep = "\t")
write.table( genesQCWTImp, file = file.path( WTDir, "geneExpTableDedupQCimp.csv"), sep = "\t")

cat( "End of imputation\n")

