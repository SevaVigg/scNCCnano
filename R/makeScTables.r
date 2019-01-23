#Preparing clean  data for WT, WTandSox10, WTandAllmutants

resDir		<- file.path( getwd(), "Res")
scTablesDir	<- file.path( resDir, "scTables")
dir.create( scTablesDir, showWarnings = FALSE)

dirQCres <- file.path( resDir, "QualityControl", "ResQC") 

genesQC	<- read.table( file = file.path( dirQCres, "NormalizedExTable.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
cellsQC	<- read.table( file = file.path( dirQCres, "cellDescripitonsDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )

rownames(genesQC)[which(rownames(genesQC)=="Kanamycin Pos")] <- "Kanamycin_Pos"

allGenes	<- rownames(genesQC)
logExps 	<- log2(1+genesQC)

allCellsDir 	<- file.path( scTablesDir, "allCells")
dir.create( allCellsDir , showWarnings = FALSE)
write.table( logExps, file = file.path( allCellsDir, "logExpTableDedupQC_all.csv"), sep = "\t")
write.table( cellsQC, file = file.path( allCellsDir, "cellDescriptionsDedupQC_all.csv"), sep = "\t")

WTDir		<- file.path( scTablesDir, "WT")
dir.create(WTDir, showWarnings = FALSE)	
logExpsWT	<- logExps[, setdiff( seq_along(colnames(logExps)), grep("(sox10|mitfa)", colnames(logExps)))]
cellsDescWT	<- cellsQC[, setdiff( seq_along(colnames(cellsQC)), grep("(sox10|mitfa)", colnames(cellsQC)))]
write.table( logExpsWT, file = file.path( WTDir, "logExpTableDedupQC_wt.csv"), sep = "\t")
write.table( cellsQC, file = file.path( WTDir, "cellDescriptionsDedupQC_wt.csv"), sep = "\t")

WT_Sox10Dir	<- file.path( scTablesDir, "WT_Sox10")
dir.create( WT_Sox10Dir, showWarnings = FALSE)
logExpsWTsox10		<- logExps[, setdiff( seq_along(colnames(logExps)), grep("(mitfa)", colnames(logExps)))]
cellsDescWTsox10	<- cellsQC[, setdiff( seq_along(colnames(cellsQC)), grep("(mitfa)", colnames(cellsQC)))]
write.table( logExpsWTsox10,   file = file.path( WT_Sox10Dir, "logExpTableDedupQC_wtSox10.csv"), sep = "\t")
write.table( cellsDescWTsox10, file = file.path( WT_Sox10Dir, "cellDescriptionsDedupQC_wtSox10.csv"), sep = "\t")

