#Preparing clean  data for WT, WTandSox10, WTandAllmutants

source("R/imputeDropouts.r")

resDir		<- file.path( getwd(), "Res")
scTablesDir	<- file.path( resDir, "scTables")
dir.create( scTablesDir, showWarnings = FALSE)

dirQCres <- file.path( resDir, "QualityControl", "ResQC") 

genesQC	<- read.table( file = file.path( dirQCres, "NormalizedExTable.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
cellsQC	<- read.table( file = file.path( dirQCres, "cellDescripitonsDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )

allGenes	<- rownames(genesQC)
logExps 	<- log2(1+genesQC)

allCellsDir 	<- file.path( scTablesDir, "allCells")
dir.create( allCellsDir , showWarnings = FALSE)

randSeedAll 	<- as.integer(Sys.time())
logExpsImp	<- imputeDropouts(logExps, randSeedAll)

cat( file = file.path( allCellsDir, "ImputeSeedAll.txt"), randSeedAll)
write.table( logExps, file = file.path( allCellsDir, "logExpTableDedupQC.csv"), sep = "\t")
write.table( logExpsImp, file = file.path( allCellsDir, "logExpTableDedupQCimp.csv"), sep = "\t")
write.table( cellsQC, file = file.path( allCellsDir, "cellDescriptionsDedupQC.csv"), sep = "\t")

WTDir		<- file.path( scTablesDir, "WT")
dir.create(WTDir, showWarnings = FALSE)	

randSeedWT 	<- as.integer(Sys.time())
logExpsWT	<- logExps[, setdiff( seq_along(colnames(logExps)), grep("(sox10|mitfa)", colnames(logExps)))]
logExpsWTImp	<- imputeDropouts(logExpsWT, randSeedWT)
cellsDescWT	<- cellsQC[, setdiff( seq_along(colnames(cellsQC)), grep("(sox10|mitfa)", colnames(cellsQC)))]

cat( file = file.path( WTDir, "ImputeSeedWT.txt"), randSeedWT)
write.table( logExpsWT, file = file.path( WTDir, "logExpTableDedupQC.csv"), sep = "\t")
write.table( logExpsWTImp, file = file.path( WTDir, "logExpTableDedupQCimp.csv"), sep = "\t")
write.table( cellsQC, file = file.path( WTDir, "cellDescriptionsDedupQC.csv"), sep = "\t")

WT_Sox10Dir	<- file.path( scTablesDir, "WT_Sox10")
dir.create( WT_Sox10Dir, showWarnings = FALSE)

randSeedWTsox10 	<- as.integer(Sys.time())
logExpsWTsox10		<- logExps[, setdiff( seq_along(colnames(logExps)), grep("(mitfa)", colnames(logExps)))]
logExpsWTsox10Imp	<- imputeDropouts(logExpsWTsox10, randSeedWTsox10)
cellsDescWTsox10	<- cellsQC[, setdiff( seq_along(colnames(cellsQC)), grep("(mitfa)", colnames(cellsQC)))]

cat( file = file.path( WT_Sox10Dir, "ImputeSeedWTsox10.txt"), randSeedWTsox10)
write.table( logExpsWTsox10,   file = file.path( WT_Sox10Dir, "logExpTableDedupQC.csv"), sep = "\t")
write.table( logExpsWTsox10Imp,   file = file.path( WT_Sox10Dir, "logExpTableDedupQCimp.csv"), sep = "\t")
write.table( cellsDescWTsox10, file = file.path( WT_Sox10Dir, "cellDescriptionsDedupQC.csv"), sep = "\t")

