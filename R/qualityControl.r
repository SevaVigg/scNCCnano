if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}

if (!require("vsn")){
BiocManager::install("vsn", version = "3.8")
library(vsn)}

if(!require("NanoStringNorm")){
  install.packages("NanoStringNorm", dependencies = TRUE)
}

require(preprocessCore)
require(data.table)

resDir		<- file.path( getwd(), "Res")
inTablesDir	<- file.path( resDir, "InitialTables")
qualDir 	<- file.path( resDir, "QualityControl")
plotDir 	<- file.path( qualDir, "Plot")
resQCDir	<- file.path( qualDir, "ResQC")

dir.create( qualDir, showWarnings = FALSE)
dir.create( plotDir, showWarnings = FALSE)
dir.create( resQCDir, showWarnings = FALSE)

qualContLogFile	<- file(paste0(resQCDir, .Platform$file.sep, "qualContLog.txt"), open = "w")

Genes 	<- read.csv( file = file.path( inTablesDir, "expressionTableDedup.csv" ) , stringsAsFactors = FALSE)
Cells	<- read.csv( file = file.path( inTablesDir, "cellDescripitonsDedup.csv") , stringsAsFactors = FALSE)
Probes	<- read.csv( file = file.path( inTablesDir, "ProbesDescripitonsDedup.csv") , stringsAsFactors = FALSE)

geneNames <- Genes[,1]
cellNames <- Cells[,1]

rownames(Cells)	<- cellNames
rownames(Genes)	<- geneNames
Cells[,1]	<- NULL
Genes[,1]	<- NULL

# We have a number of x18 cells with missing mitfa values; these values are updated as medians
# X18_mitfa 	<- Genes["mitfa", grep("X18", colnames(Genes))]
# X18_mitfa_med	<- median(as.numeric(X18_mitfa["mitfa", !is.na(X18_mitfa[1,])]))
# Genes["mitfa", is.na(Genes["mitfa",])] <- X18_mitfa_med
# cat(file = qualContLogFile, ncol(X18_mitfa), " cells with missing mitfa values updated as medians\n")

#remove mitfa gene

colnames(Genes) <- colnames(Cells) <- paste0(Cells["hpf",],"_", Cells["CellType",])

geneMatrix	<- as.matrix(Genes)

normGeneMatrix	<- normalize.quantiles(geneMatrix)

rownames(normGeneMatrix) <- Probes[, "Gene.Name"]
colnames(normGeneMatrix) <- Cells["FileName", ]

batches		<- unique(unlist(Cells["FileName",]))
batchDates	<- unlist(lapply(batches, function(x) Cells["dateEx", grep(x, Cells["FileName",], fixed = TRUE)][1]))
batchTypes	<- unlist(lapply(batches, function(x) Cells["CellType", grep(x, Cells["FileName",], fixed = TRUE)][1]))
batchHpf	<- unlist(lapply(batches, function(x) Cells["hpf", grep(x, Cells["FileName",], fixed = TRUE)][1]))

batchVals	<- lapply(batches, function(x) geneMatrix[, which(Cells["FileName", ]==x)])
batchNormVals 	<- lapply(batches, function(x) normGeneMatrix[, which(Cells["FileName", ]==x)])

names(batchNormVals) <- as.character(batches)
names(batchVals)     <- as.character(batches)

batchMedVals	<- unlist(lapply(batchVals, median))
batchMedNormVals<- unlist(lapply(batchNormVals, median))

qualMatrix	<- as.data.frame(rbind(batchHpf, batchDates, batchTypes, batchMedVals, batchMedNormVals), stringsAsFactors = FALSE)
names(qualMatrix)	<- batches
qualMatrix	<- qualMatrix[, order(as.numeric(qualMatrix["batchMedNormVals",]))]

png(paste0(plotDir, .Platform$file.sep,"batchMedians.png"))
	plot(as.numeric(qualMatrix["batchMedNormVals", ]), main = "Medians of batches", xlab = "Batches", ylab = "CountMedian")
	abline(h = 40, col = "red")
	abline(h = min(as.numeric(qualMatrix["batchMedNormVals", qualMatrix["batchTypes",]=="M"])), col = "blue")
dev.off()

batchVecs <- lapply(batchVals, function(x) log(as.vector(x)))

png(paste0(plotDir, .Platform$file.sep, "LogNotNormedBatchBoxPlot.png"), width = 960)
	boxplot(batchVecs[colnames(qualMatrix)], 
		main = "Not normalized, not filtered batches (log-scale)", 
		las =2 , 
		cex.axis = 0.7, 
		boxwex = 1.2, 
		at = seq(1, 2*length(batches), 2))
dev.off()

batchProbl 	<- qualMatrix[, which(as.numeric(qualMatrix["batchMedNormVals",])<40)]	#threshold to keep M (MedNormVal ~ 60)
cellsProbl_Ind	<- which(Cells["FileName",] %in% colnames(batchProbl))


Genes_f		<- Genes[,-cellsProbl_Ind]
Cells_f		<- Cells[,-cellsProbl_Ind]

cat(file = qualContLogFile, ncol(Genes_f), " cells remaining after removing batches with low medians\n")

posSpikes	<- c(128, 32, 8, 2, 0.5, 0.125)						# Pos probes in fM
posProbes	<- c("POS_A", "POS_B", "POS_C", "POS_D", "POS_E", "POS_F")

coefs		<- apply( log2(Genes_f[posProbes,]), 2, function(x) lm( x~log2(posSpikes))$coefficients[2] )

sortLogCoefs	<- sort(coefs)

png(paste0(plotDir, .Platform$file.sep, "LogSortPosCoefs.png"))
	plot(sortLogCoefs, main = "Positive quality control coef values")
	abline( h = 1,   col = "green")
	abline( h = 1.1, col = "red")
	abline( h = 0.9, col = "red")
dev.off()

keepCoefs <- which(coefs > 0.9 & coefs < 1.1)

Genes_ff	<- Genes_f[,keepCoefs]
Cells_ff	<- Cells_f[,keepCoefs]

cat(file = qualContLogFile, ncol(Genes_ff), " cells remaining after removing cells with poor positive control values\n")

KanamycinPos_Distr		<- log( as.numeric(Genes_ff["Kanamycin Pos", ]))
names(KanamycinPos_Distr)	<- colnames(Genes_ff)

KanamycinPosTopQuant	<- 0.995
KanamycinPosBotQuant	<- 0.005

png( file.path( plotDir, "Kanamycin_Distr_plot.png"))
	plot(sort( KanamycinPos_Distr))
	abline( h = quantile( KanamycinPos_Distr, KanamycinPosTopQuant), col = "red")
	abline( h = quantile( KanamycinPos_Distr, KanamycinPosBotQuant), col = "red")
	abline( h = min( KanamycinPos_Distr[ grep("M", names( KanamycinPos_Distr))]), col = "black")
	abline( h = max( KanamycinPos_Distr[ grep("M", names( KanamycinPos_Distr))]), col = "black")
	abline( h = max( KanamycinPos_Distr[ grep("I", names( KanamycinPos_Distr))]), col = "cyan")
	abline( h = min( KanamycinPos_Distr[ grep("I", names( KanamycinPos_Distr))]), col = "cyan")
dev.off()


KanamycinPos_keep	<- which( KanamycinPos_Distr > quantile( KanamycinPos_Distr, KanamycinPosBotQuant)
				& KanamycinPos_Distr < quantile( KanamycinPos_Distr, KanamycinPosTopQuant))

Genes_fff	<- Genes_ff[,KanamycinPos_keep]
Cells_fff	<- Cells_ff[,KanamycinPos_keep]

cat(file = qualContLogFile, ncol(Genes_fff), " cells remaining after removing cells with Kanamycin Pos dropouts\n")



rpl13_Distr	<- log(as.numeric(Genes_fff["rpl13", ]))
names(rpl13_Distr) <- colnames(Genes_fff)

rpl13TopQuant	<- 0.97
rpl13BotQuant	<- 0.03

rpl13_keep	<- which( rpl13_Distr > quantile(rpl13_Distr, rpl13BotQuant) & rpl13_Distr < quantile(rpl13_Distr, rpl13TopQuant))

png( file.path( plotDir, "rpl13_Distr_plot.png"))
	plot(sort(rpl13_Distr))
	abline( h = quantile( rpl13_Distr, rpl13TopQuant), col = "red")
	abline( h = quantile( rpl13_Distr, rpl13BotQuant), col = "red")
	abline( h = min( rpl13_Distr[ grep("M", names(rpl13_Distr))]), col = "black")
	abline( h = max( rpl13_Distr[ grep("M", names(rpl13_Distr))]), col = "black")
	abline( h = max( rpl13_Distr[ grep("I", names(rpl13_Distr))]), col = "cyan")
	abline( h = min( rpl13_Distr[ grep("I", names(rpl13_Distr))]), col = "cyan")
dev.off()

Genes_ffff	<- Genes_fff[,rpl13_keep]
Cells_ffff	<- Cells_fff[,rpl13_keep]

cat(file = qualContLogFile, ncol(Genes_ffff), " cells remaining after removing cells with rpl13 dropouts\n")

normIteration <- 0
repeat{									       #iterations over the background level								      
normIteration <- normIteration + 1
cat("Starting itration ", normIteration, "\n")

NanoTable	<- cbind(Probes[,"Class.Name"], Probes[,"Gene.Name"], Probes[,"Accession.."], Genes_ffff, stringsAsFactors = FALSE)
colnames(NanoTable)[1:3] <- c("Code.Class", "Name", "Accession")

NanoTable[c("Kanamycin Pos", "rpl13"), "Code.Class"] <- "Housekeeping"
NanoTableNormed <- NanoStringNorm(x = NanoTable, CodeCount = "sum", Background = "mean", SampleContent = "housekeeping.geo.mean")

nf	 	<- NanoTableNormed$sample.summary.stats.norm$pos.norm.factor	#recommended values for normfactors are between 0.3 and 3
bg	 	<- NanoTableNormed$sample.summary.stats.norm$background.level	#recommended values for background are within 3 sd from the mean
sm	 	<- NanoTableNormed$sample.summary.stats.norm$Sample.Missing	#recommended values for the number of missing samples is less than 0.9
rc	 	<- NanoTableNormed$sample.summary.stats.norm$rna.content	#recommended values for RNA content are within 3 sd from the mean
nf_drop		<- which( nf < 0.3 | nf > 3)
sm_drop		<- which( sm > 0.9 )
bg_drop		<- which( abs( bg - mean(bg)) > 3*sd(bg))
rc_drop		<- which( abs( rc - mean(rc)) > 3*sd(rc))

cat( file = qualContLogFile, "Iteration ", normIteration, length(nf_drop), " cells with inadequate norm factors\n")
cat( file = qualContLogFile, "Iteration ", normIteration, length(bg_drop), " cells with adequate backgrounds\n")
cat( file = qualContLogFile, "Iteration ", normIteration, length(sm_drop), " cells with less than 0.9 missing samples\n")
cat( file = qualContLogFile, "Iteration ", normIteration, length(rc_drop), " cells with not too large RNA content\n")

"%u%"		<- union							#union() supports only two arguments
norm_drop	<- nf_drop %u% sm_drop %u% bg_drop %u% rc_drop

if( length(norm_drop) == 0){break}

Genes_ffff	<- Genes_ffff[, -norm_drop]
Cells_ffff	<- Cells_ffff[, -norm_drop]

}										# repeat

cat( file = qualContLogFile, ncol(Genes_ffff), " cells remaining after removing cells with poor normalization statistics\n")

NanoTable	<- cbind(Probes[,"Class.Name"], Probes[,"Gene.Name"], Probes[,"Accession.."], Genes_ffff, stringsAsFactors = FALSE)
colnames(NanoTable)[1:3] <- c("Code.Class", "Name", "Accession")
NanoTable[c("Kanamycin Pos", "rpl13"), "Code.Class"] <- "Housekeeping"
Genes_n		 <- NanoStringNorm(x = NanoTable, CodeCount = "sum", Background = "mean", SampleContent = "housekeeping.geo.mean", return.matrix.of.endogenous.probes = TRUE)

geneMatrix_n	<- as.matrix(Genes_n)

batches_n	<- unique(unlist(Cells_ffff["batch",]))
batchDates_n	<- unlist(lapply(batches_n, function(x) Cells_ffff["dateEx", grep(x, Cells_ffff["batch",])][1]))
batchTypes_n	<- unlist(lapply(batches_n, function(x) Cells_ffff["CellType", grep(x, Cells_ffff["batch",])][1]))
batchHpf_n	<- unlist(lapply(batches_n, function(x) Cells_ffff["hpf", grep(x, Cells_ffff["batch",])][1]))

batchVals_n	<- lapply(batches_n, function(x) geneMatrix_n[, which(Cells_ffff["batch", ]==x)])

names(batchVals_n) <- as.character(batches_n)

batchMedVals_n	<- unlist(lapply(batchVals_n, median))

qualMatrix_n	<- as.data.frame(rbind(batchHpf_n, batchDates_n, batchTypes_n, batchMedVals_n), stringsAsFactors = FALSE)
names(qualMatrix_n)	<- batches_n
qualMatrix_n	<- qualMatrix_n[, order(as.numeric(qualMatrix_n["batchMedVals_n",]))]

png(paste0(plotDir, .Platform$file.sep,"NormalizedBatchMedians.png"))
	plot(as.numeric(qualMatrix_n["batchMedVals_n", ]), main = "Medians of normalized batches", xlab = "Batches", ylab = "CountMedian")
dev.off()

batchVecs_n <- lapply(batchVals_n, function(x) log(as.vector(x)+1))

png(paste0(plotDir, .Platform$file.sep, "LogNormalizedBatchBoxPlot.png"), width = 960)
	boxplot(batchVecs_n[colnames(qualMatrix_n)], 
		main = "Normalized and filtered batches (log-scale)", 
		las =2 , 
		cex.axis = 0.7, 
		boxwex = 1.2, 
		at = seq(1, 2*length(batches_n), 2))
dev.off()


write.table(Genes_n, file = file.path(resQCDir, "NormalizedExTable.csv"), sep = "\t")



cells_dt 	<- as.data.table( t(Cells_ffff))
cells_tbl	<- cells_dt[, .N, .(CellType,hpf)]

write.table(cells_tbl, file=qualContLogFile, sep = "\t")

close(qualContLogFile)	

write.table( Genes_ffff, file = file.path( resQCDir, "expressionTableDedupQC.csv"), sep = "\t" )
write.table( Cells_ffff, file = file.path( resQCDir, "cellDescripitonsDedupQC.csv"),sep = "\t" )









 


		
 

