#This script removes "non-interesting" genes from the dataset
#Preparing clean  data for WT, WTandSox10, WTandAllmutants
#This script is the next in the pipeline after the qualityControl.r
#this script also imputes for dropouts

source("R/imputeDropouts.r")

resDir		<- file.path( getwd(), "Res")
scTablesDir	<- file.path( resDir, "scTables")
dir.create( scTablesDir, showWarnings = FALSE)

dirQCres <- file.path( resDir, "QualityControl", "ResQC") 

genesQC	<- read.table( file = file.path( dirQCres, "NormalizedExTable.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
cellsQC	<- read.table( file = file.path( dirQCres, "cellDescripitonsDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )

allGenes	<- rownames(genesQC)
logExps 	<- log2(1+genesQC)

#remove genes with too low values
genesToSubset	<- c("csf1r", "sox5", "dpf3", "ets1a", "fgfr3_v2", "mycl1a", "smad9", "pax3_v2", "hbp1")
logExps		<- logExps[ !(rownames(logExps) %in% genesToSubset), ]
genesQC		<- genesQC[ !(rownames(genesQC) %in% genesToSubset), ]

allCellsDir 	<- file.path( scTablesDir, "allCells")
dir.create( allCellsDir , showWarnings = FALSE)

randSeedAll 	<- as.integer(Sys.time())
logExpsImp	<- imputeDropouts(logExps, randSeedAll)

#noiseTol	<- log2(19)
noiseTol	<- 0

#remove values, that are too close to zero
logExpsImp	<- apply( logExpsImp, c(1,2), function(x) if(x>noiseTol) x else 0)
logExpsImp	<- logExpsImp[, which( apply( logExpsImp, 2, sum) >0)]
cellsQCImp	<- cellsQC[, colnames(logExpsImp)]
genesQCImp	<- genesQC[, colnames(logExpsImp)]

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

logExpsWTImp	<- apply( logExpsWTImp, c(1,2), function(x) if(x>noiseTol) x else 0)
logExpsWTImp	<- logExpsWTImp[, which( apply( logExpsWTImp, 2, sum) >0)]
cellsQCWTImp	<- cellsQC[, colnames(logExpsWTImp)]
genesQCWTImp	<- genesQC[, colnames(logExpsWTImp)]

cat( file = file.path( WTDir, "ImputeSeedWT.txt"), randSeedWT)
write.table( logExpsWT, file = file.path( WTDir, "logExpTableDedupQC.csv"), sep = "\t")
write.table( logExpsWTImp, file = file.path( WTDir, "logExpTableDedupQCimp.csv"), sep = "\t")
write.table( cellsQCWTImp, file = file.path( WTDir, "cellDescriptionsDedupQC.csv"), sep = "\t")
write.table( genesQCWTImp, file = file.path( WTDir, "geneExpTableDedupQCimp.csv"), sep = "\t")

WT_Sox10Dir	<- file.path( scTablesDir, "WT_Sox10")
dir.create( WT_Sox10Dir, showWarnings = FALSE)

randSeedWTsox10 	<- as.integer(Sys.time())
logExpsWTsox10		<- logExps[, setdiff( seq_along(colnames(logExps)), grep("(mitfa)", colnames(logExps)))]
logExpsWTsox10Imp	<- imputeDropouts(logExpsWTsox10, randSeedWTsox10)

logExpsWTsox10Imp	<- apply( logExpsWTsox10Imp, c(1,2), function(x) if(x>noiseTol) x else 0)
logExpsWTsox10Imp	<- logExpsWTsox10Imp[, which( apply( logExpsWTsox10Imp, 2, sum) >0)]
cellsWTsox10Imp		<- cellsQC[ , colnames(logExpsWTsox10Imp)]
genesQCWTsox10Imp	<- genesQC[, colnames(logExpsWTImp)]

cat( file = file.path( WT_Sox10Dir, "ImputeSeedWTsox10.txt"), randSeedWTsox10)
write.table( logExpsWTsox10,   file = file.path( WT_Sox10Dir, "logExpTableDedupQC.csv"), sep = "\t")
write.table( logExpsWTsox10Imp,   file = file.path( WT_Sox10Dir, "logExpTableDedupQCimp.csv"), sep = "\t")
write.table( cellsWTsox10Imp, file = file.path( WT_Sox10Dir, "cellDescriptionsDedupQC.csv"), sep = "\t")
write.table( genesQCWTsox10Imp, file = file.path( WT_Sox10Dir, "geneExpTableDedupQCimp.csv"), sep = "\t")

