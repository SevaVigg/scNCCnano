if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}

if (!require("vsn")){
BiocManager::install("vsn", version = "3.8")
library(vsn)}

require(gtools)
require(dplyr)
require(tidyr)
require(ggplot2)

if(!require("NanoStringNorm")){
  install.packages("NanoStringNorm", dependencies = TRUE)
}

require( NanoStringNorm)
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

qualDF		<- data.frame( t(Genes), t(Cells))

#First we filter cells with probes

#now check negative probes
negProbes	<- grep("NEG", colnames( qualDF), value= TRUE)
negProbThres	<- quantile( unlist( qualDF[ ,negProbes]), 0.97)
negProbTable	<- apply( qualDF[ ,negProbes], c(1,2), function(x) if( x > negProbThres) 1 else 0)
qualDF$negProbSum	<- apply( negProbTable, 1, sum)

negProbDF 	<- gather( qualDF[ , negProbes], negProbe, negProbVal)
negThreshold 	<- quantile( negProbDF$negProbVal, 0.97)

ggplot(data = negProbDF, aes( x = negProbVal)) +
	geom_vline( xintercept = negThreshold, col = "red", linetype = "longdash") +
	geom_histogram( fill = "blue", binwidth = 2) +
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 0, 
		hjust = 1, 
		vjust = 0.5, 
		size = 12),
		axis.text.y = element_text( size = 12)) +
	scale_y_continuous( name = "Counts", breaks = seq(1000, 5000, by = 1000), expand = c(0,0)) +
	scale_x_continuous( name = "Negative probe value", expand = c(0,0) ) +
	ggtitle( "Histogram of negative probe expression, \naggregation by all genes and negative probes") 

#now check positive probes
posSpikes		<- c(128, 32, 8, 2, 0.5, 0.125)						# Pos probes in fM
posProbes		<- grep("POS", rownames(Genes), value= TRUE)
qualDF$posProbCoef 	<- apply( log2(Genes[ posProbes, ]), 2, function(x) lm( x~log2(posSpikes))$coefficients[2] )
leftPosProbThrsld	<- quantile( qualDF$posProbCoef, 0.01)
rightPosProbThrsld	<- quantile( qualDF$posProbCoef, 0.99)

ggplot( data = qualDF[ , "posProbCoef", drop = FALSE], aes( x = posProbCoef)) +
	geom_histogram( fill = "blue", binwidth = 0.01) +
	geom_vline( xintercept = leftPosProbThrsld, col = "red", linetype = "longdash") +
	geom_vline( xintercept = rightPosProbThrsld, col = "red", linetype = "longdash") +
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 0, 
		hjust = 1, 
		vjust = 0.5, 
		size = 12),
		axis.text.y = element_text( size = 12)) +
	scale_y_continuous( name = "Counts", expand = c(0,0)) +
	scale_x_continuous( name = "Positive probe log(value) regression coefficient", expand = c(0,0), limits = c(0, 1.5) ) +
	ggtitle( "Histogram of positive probe log-expression, \nregression coefficient for spiked controls") 

#remove cells with 3 or more counting negative probes

negProbCellIndex <- qualDF$negProbSum > 2
gualDF		 <- qualDF[ -negProbCellIndex, ]

posProbCellIndex <- qualDF$posProbCoef < leftPosProbThrsld | qualDF$posProbCoef > rightPosProbThrsld
qualDF		 <- qualDF[ -posProbCellIndex, ]

#remove genes with too low values
#now we remove cells with a very low expression for all genes.
testGenes 	<- setdiff( geneNames, c( "Kanamycin Pos", "rpl13", grep("(NEG_|POS_)", geneNames, value = TRUE)))
log2Exps	<- log2( Genes[ testGenes, ])	
dens		<- density( t( log2Exps ))
expThreshold	<- optimize(approxfun(dens$x,dens$y),interval=c(5,15))$minimum
geneThreshold	<- 2/3*expThreshold
poorCells 	<- which(sapply(log2Exps, function(x) sum(x > geneThreshold) < 3))  #cells with poor values for all genes but houskeeping
qualDF		<- qualDF[ -poorCells, ]

#make Probe value plot and Probe value 
geneAverageData 		<- data.frame( gene <- factor(rownames(Genes), levels = rownames(Genes)), avExp = apply( qualDF[ , rownames(Genes) ], 2, mean))
geneAverageData$top10median	<- apply( qualDF[ , rownames(Genes) ] , 2, function(x) median(tail(sort(unlist(x)), 10))) 

ggplot( data = geneAverageData) + 
	aes( x = gene, y = avExp) +
	geom_col(fill = "blue") +
	geom_hline( yintercept = geneThreshold, color = "red") + 
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 90, 
		hjust = 1, 
		vjust = 0.5, 
		size = 7)) +
	scale_x_discrete( name = "Gene" ) +
	scale_y_continuous( name = "Expression", expand = c(0,0)) +
	ggtitle( "Average gene expression over all cells in the dataset")

ggplot( data = geneAverageData) + 
	aes( x = gene, y = log10(avExp)) +
	geom_col(fill = "blue") +
	geom_hline( yintercept = log(geneThreshold), color = "red") + 
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 90, 
		hjust = 1, 
		vjust = 0.5, 
		size = 7)) +
	scale_x_discrete( name = "Gene" ) +
	scale_y_continuous( name = "Expression", expand = c(0,0)) +
	ggtitle( "Average gene expression over all cells in the dataset")

ggplot( data = geneAverageData) + 
	aes( x = gene, y = log10(avExp)) +
	geom_col(fill = "blue") +
	geom_hline( yintercept = log10(geneThreshold), color = "red") + 
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 90, 
		hjust = 1, 
		vjust = 0.5, 
		size = 7)) +
	scale_x_discrete( name = "Gene" ) +
	scale_y_continuous( name = "Expression", expand = c(0,0)) +
	ggtitle( "Average gene expression over all cells in the dataset, log-scale")

ggplot( data = geneAverageData) + 
	aes( x = gene, y = log10( top10median)) +
	geom_col(fill = "blue") +
	geom_hline( yintercept = log10(geneThreshold), color = "red") + 
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 90, 
		hjust = 1, 
		vjust = 0.5, 
		size = 7)) +
	scale_x_discrete( name = "Gene" ) +
	scale_y_continuous( name = "Expression", expand = c(0,0)) +
	ggtitle( "Median expression over top 10 cells, \nexpressing the gene in the dataset, log-scale")

poorGenes 	<- rownames(log2Exps)[which( apply(log2Exps, 1, function(x) sum(x > geneThreshold) < 5))]  #genes poorly expressed in all cell types
genes2exclude	<- character(0)
#genes2exclude	<- c("csf1r", "sox5", "dpf3", "ets1a", "fgfr3_v2", "mycl1a", "smad9", "pax3_v2", "hbp1")
poorGenes	<- c(poorGenes, genes2exclude)
goodGenes	<- setdiff( geneNames, poorGenes)
Genes_p		<- Genes_p[ goodGenes, ]
Probes		<- Probes[ Probes[ ,"Gene.Name"] %in% goodGenes, ]

geneMatrix	<- as.matrix(Genes_p)
normGeneMatrix	<- normalize.quantiles(geneMatrix)

rownames(normGeneMatrix) <- rownames(Genes_p)
colnames(normGeneMatrix) <- Cells_p["FileName", ]

qualDF		<- data.frame( t(Genes_p), t(Cells_p))		#this data frame contains cells and cell-related information. Cells are in rows

batchList	<- split( qualDF[, goodGenes], list(qualDF$FileName, qualDF$hpf), drop = TRUE)
batchMedian	<- sapply( batchList, function(x) median(unlist(x)))

qualNormDF	<- data.frame( t(normGeneMatrix), t(Cells_p))
batchNormList	<- split( qualNormDF[ , goodGenes], list(qualNormDF$FileName, qualDF$hpf), drop = TRUE)

batchDF		<- data.frame( median = batchMedian, batch = names(batchMedian)) #this data frame contains batches and batch-related information
batchDF$normMedian 	<- sapply( batchNormList, function(x) median(unlist(x)))
batchDF$fileName	<- sapply(strsplit(as.character(batchDF$batch), "[.]"), function(b) paste(head( b, length(b)-1), collapse = "."))
batchDF$cellType	<- sapply(strsplit(as.character(batchDF$batch), "[.]"), function(b){
					fn <- paste(head( b, length(b)-1), collapse = ".") 
					return( unique( qualDF[qualDF$FileName == fn, "CellType"]))})
batchDF$hpf	<- sapply(strsplit(as.character(batchDF$batch), "[.]"), function(b){
					hp <- paste(tail( b, 1), collapse = ".") 
					return( unique( qualDF[qualDF$hpf == hp, "hpf"]))})

batchDF$batch 	<- factor(batchDF$batch, levels = batchDF$batch) #ggplot requires an ordered factor as x, otherwise it reorders columns
batchDF		<- batchDF[order(batchDF$cellType),]		 #but the ordering we need is according the cell type

normMedianThreshold <- min(batchDF$normMedian[ batchDF$cellType=="Tl"]) #keeping batches that retain all Tails samples

batchMedianPlot <- ggplot( data = batchDF) +
	aes( x = batch, y = normMedian) + 
	geom_col( fill = "blue") +
	geom_hline( yintercept = normMedianThreshold, color = "magenta") +
	theme( 	panel.background = element_rect( fill = "gray80"),
		title 		= element_text( size = 17),
		axis.title = element_text( size = 15),
		axis.text.x = element_text(
			angle = 90,
			hjust = 1,
			vjust = 0.5,
			size = 11), 
		axis.text.y = element_text( size = 13)
		) +
	scale_x_discrete( name = "batch", labels = paste(batchDF$cellType, batchDF$hpf, sep = ".")) +
	scale_y_continuous( name = "Median expression") +
	ggtitle("Median gene expression over all genes and all cells in the batch. \nQuantile normalization" )

png(paste0(plotDir, .Platform$file.sep,"batchMedians.png"), width = 960)
	(batchMedianPlot)
dev.off()

batchVecs <- lapply(batchList, function(x) log(unlist(x)))	# legacy variable, probably not needed

#png(paste0(plotDir, .Platform$file.sep, "LogNotNormedBatchBoxPlot.png"), width = 960)
#	boxplot(batchVecs[colnames(qualMatrix)], 
#		main = "Not normalized, not filtered batches (log-scale)", 
#		las =2 , 
#		cex.axis = 0.7, 
#		boxwex = 1.2, 
#		at = seq(1, 2*length(batches), 2))
#dev.off()

batchProbl 	<- batchDF[ batchDF$normMedian < normMedianThreshold, ]	#Problematic batches; threshold to keep Tl (MedNormVal ~ 90)
cellsProbl_Ind	<- which( qualDF$FileName %in% batchProbl$fileName & qualDF$hpf %in% batchProbl$hpf)

qualDF_f	<- qualDF[ -cellsProbl_Ind, ]
Genes_f		<- t(qualDF_f[ , goodGenes])

cat(file = qualContLogFile, nrow( qualDF_f), " cells remaining after removing batches with low medians\n")

ggplot( 


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
cat("Starting iteration ", normIteration, "\n")

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









 


		
 

