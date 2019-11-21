if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}

if (!require("vsn")){
BiocManager::install("vsn", version = "3.8")
library(vsn)}

require(gtools)
require(dplyr)
require(tidyr)
require(ggplot2)
require(grid)

if(!require("gridGraphics")){
install.packages("gridGraphics")
library("gridGraphics")}

require(cowplot)

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

#First we filter cells according to control probes

#now check negative probes
negProbes	<- grep("NEG", colnames( qualDF), value= TRUE)
negProbThres	<- quantile( unlist( qualDF[ ,negProbes]), 0.97)
negProbTable	<- apply( qualDF[ ,negProbes], c(1,2), function(x) if( x > negProbThres) 1 else 0)
qualDF$negProbSum	<- apply( negProbTable, 1, sum)

negProbDF 	<- gather( qualDF[ , negProbes], negProbe, negProbVal)
negThreshold 	<- quantile( negProbDF$negProbVal, 0.97)

negProbCellIndex	<- which(qualDF$negProbSum < 3)
nCellNeg		<- length(negProbCellIndex)
nCellNegGrob 		<- grobTree( textGrob( paste0(nCellNeg, " cells remaining"), x=0.5,  y=0.95, hjust=0,
  				gp=gpar(col="black", fontsize=13, fontface="bold")))

negProbDistrPlot <- ggplot(data = negProbDF, aes( x = negProbVal)) +
	geom_vline( xintercept = negThreshold, col = "red", linetype = "longdash") +
	geom_histogram( fill = "blue", binwidth = 2) +
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 0, 
		hjust = 1, 
		vjust = 0.5, 
		size = 12),
		axis.text.y = element_text( size = 12),
		plot.title  = element_text( hjust = 0)) +
	scale_y_continuous( name = "Counts", breaks = seq(1000, 5000, by = 1000), expand = c(0,0)) +
	scale_x_continuous( name = "Negative probe value", expand = c(0,0) ) +
	annotation_custom( nCellNegGrob) +
	ggtitle( "Negative probes, \nall neg probes in a cell aggregated") 

qualDF		 	<- qualDF[ negProbCellIndex, ]

#now check positive probes
posSpikes		<- c(128, 32, 8, 2, 0.5, 0.125)						# Pos probes in fM
posProbes		<- grep("POS", colnames( qualDF), value= TRUE)
qualDF$posProbCoef 	<- apply( log2( qualDF[ , posProbes ]), 1, function(x) lm( x~log2(posSpikes))$coefficients[2] )
leftPosProbThrsld	<- quantile( qualDF$posProbCoef, 0.01)
rightPosProbThrsld	<- quantile( qualDF$posProbCoef, 0.99)

posProbCellIndex 	<- which(qualDF$posProbCoef > leftPosProbThrsld & qualDF$posProbCoef < rightPosProbThrsld)
nCellPos		<- length( posProbCellIndex)
nCellPosGrob 		<- grobTree( textGrob( paste0(nCellPos, " cells remaining"), x=0.1,  y=0.95, hjust=0,
  				gp=gpar(col="black", fontsize=13, fontface="bold")))

posProbDistrPlot <- ggplot( data = qualDF[ , "posProbCoef", drop = FALSE], aes( x = posProbCoef)) +
	geom_histogram( fill = "blue", binwidth = 0.01) +
	geom_vline( xintercept = leftPosProbThrsld, col = "red", linetype = "longdash") +
	geom_vline( xintercept = rightPosProbThrsld, col = "red", linetype = "longdash") +
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 0, 
		hjust = 1, 
		vjust = 0.5, 
		size = 12),
		axis.text.y = element_text( size = 12),
		plot.title  = element_text( hjust = 0)) +
	scale_y_continuous( name = "Counts", expand = c(0,0)) +
	scale_x_discrete( name = "log10 pos regression coef.", expand = c(0,0), limits = c(0, 1.5) ) +
	annotation_custom( nCellPosGrob) +
	ggtitle( "Positive probes, \nlog2 expression regression coefficient") 
#function(){ plot(dens, main = "Combined expression density of exogenous probes", xlab = "Expression", ylab = "Density")
#			abline( v = geneThreshold, col = "red", lty = 5)
#			text( 5.3, 0.62, labels = paste0(nCells, " cells remaining"))
#			}	


qualDF		 <- qualDF[ posProbCellIndex, ]

testGenes 	<- setdiff( geneNames, c( "Kanamycin Pos", "rpl13", grep("(NEG_|POS_)", geneNames, value = TRUE)))
log10Exps	<- log10( Genes[ testGenes, ])	

dens		<- density( t( log10Exps )) 
expMinimal	<- optimize(approxfun(dens$x,dens$y),interval=c(2,4))$minimum
geneThreshold	<- expMinimal
poorCells 	<- which(sapply(log10Exps, function(x) sum(x > geneThreshold) < 3))  #cells with poor values for all genes but houskeeping
qualDF		<- qualDF[ -poorCells, ]
nCells		<- nrow( qualDF)
nCellsExpGrob 	<- grobTree( textGrob( paste0(nCells, " cells remaining"), x=0.5,  y=0.95, hjust=0,
  				gp=gpar(col="black", fontsize=13, fontface="bold")))


densPlot	<- ggplot( data = data.frame( dens$x, dens$y), mapping = aes( dens$x, dens$y)) + 
			geom_line() +
			geom_vline( xintercept = geneThreshold, color = "red") +
			scale_x_continuous( name = "Expression" ) +
			scale_y_continuous( name = "cell density") +
			annotation_custom( nCellsExpGrob) +
			ggtitle( "Density of gene expression \nover exogenous probes")


#remove genes with very low expression values in all cells
geneAverageData 		<- data.frame( gene <- factor(rownames(Genes), levels = rownames(Genes)), avExp = apply( qualDF[ , rownames(Genes) ], 2, mean))
geneAverageData$top5		<- apply( qualDF[ , rownames(Genes) ] , 2, function(x) head(tail(sort(unlist(x)), 5), 1))	# smallest of the top five

poorGenes 	<- rownames(log10Exps)[which( apply(log10Exps, 1, function(x) sum(x > geneThreshold) < 6))]  #genes poorly expressed in all cell types
genes2exclude	<- character(0)
#genes2exclude	<- c("csf1r", "sox5", "dpf3", "ets1a", "fgfr3_v2", "mycl1a", "smad9", "pax3_v2", "hbp1")
poorGenes	<- c(poorGenes, genes2exclude)
goodGenes	<- setdiff( geneNames, poorGenes)
poorGenesI	<- which( colnames( qualDF) %in% poorGenes)

geneAverageData$filterColor	<- sapply( rownames(Genes), function(x) if( x %in% poorGenes) "red" else "blue")

geneTop5CellPlot <- ggplot( data = geneAverageData) + 
	aes( x = gene, y = log10( top5 )) +
	geom_col(fill = geneAverageData$filterColor) +
	geom_hline( yintercept = geneThreshold, color = "red") + 
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 90, 
		hjust = 1, 
		vjust = 0.5, 
		size = 7),
		plot.title  = element_text( hjust = 0)) +
	scale_x_discrete( name = "Gene" ) +
	scale_y_continuous( name = "Expression", expand = c(0,0)) +
	ggtitle( "5th rank statistics of gene expression")

#now we remove cells with a very low expression for all genes.


#make Probe value plot and Probe value 
geneAverageDistributionPlot <- ggplot( data = geneAverageData) + 
	aes( x = gene, y = avExp) +
	geom_col( fill = geneAverageData$filterColor) +
	geom_hline( yintercept = 10^geneThreshold, color = "red") + 
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 90, 
		hjust = 1, 
		vjust = 0.5, 
		size = 7),
		plot.title  = element_text( hjust = 0)) +
	scale_x_discrete( name = "Gene" ) +
	scale_y_continuous( name = "Expression", expand = c(0,0)) +
	ggtitle( "Expression average over cells")

geneLogAverageDistributionPlot <- ggplot( data = geneAverageData) + 
	aes( x = gene, y = log10(avExp)) +
	geom_col(fill = geneAverageData$filterColor) +
	geom_hline( yintercept = geneThreshold, color = "red") + 
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 90, 
		hjust = 1, 
		vjust = 0.5, 
		size = 7),
		plot.title  = element_text( hjust = 0)) +
	scale_x_discrete( name = "Gene" ) +
	scale_y_continuous( name = "Log 10 expression", expand = c(0,0)) +
	ggtitle( "log10 expression average over cells")

geneTop10CellMedianPlot <- ggplot( data = geneAverageData) + 
	aes( x = gene, y = log10( top10median)) +
	geom_col(fill = geneAverageData$filterColor) +
	geom_hline( yintercept = geneThreshold, color = "red") + 
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 90, 
		hjust = 1, 
		vjust = 0.5, 
		size = 7),
		plot.title  = element_text( hjust = 0)) +
	scale_x_discrete( name = "Gene" ) +
	scale_y_continuous( name = "Expression", expand = c(0,0)) +
	ggtitle( "Median log10 expression in top 10 cells")

#Now exclude genes that are found in small number of cells
qualDF		<- qualDF[ , -poorGenesI]
Probes		<- Probes[ Probes[ ,"Gene.Name"] %in% goodGenes, ]

geneMatrix	<- t(as.matrix( qualDF[ , goodGenes]))
normGeneMatrix	<- normalize.quantiles(geneMatrix)

rownames(normGeneMatrix) <- rownames( geneMatrix)
colnames(normGeneMatrix) <- colnames( geneMatrix)

#qualDF		<- data.frame( t(Genes_p), t(Cells_p))		#this data frame contains cells and cell-related information. Cells are in rows

batchList	<- split( qualDF[, goodGenes], list(qualDF$FileName, qualDF$hpf), drop = TRUE)
batchMedian	<- sapply( batchList, function(x) median(unlist(x)))

qualNormDF	<- qualDF; qualNormDF[ , goodGenes] <- t(normGeneMatrix)
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

normMedianThreshold <- 90 #keeping batches that retain all control samples

#batchMedianPlot <- ggplot( data = batchDF) +
#	aes( x = batch, y = normMedian) + 
#	geom_col( fill = "blue") +
#	geom_hline( yintercept = normMedianThreshold, color = "magenta") +
#	theme( 	panel.background = element_rect( fill = "gray80"),
#		title 		= element_text( size = 17),
#		axis.title = element_text( size = 15),
#		axis.text.x = element_text(
#			angle = 90,
#			hjust = 1,
#			vjust = 0.5,
#			size = 11), 
#		axis.text.y = element_text( size = 13),
#		plot.title  = element_text( hjust = 0)
#		) +
#	scale_x_discrete( name = "batch", labels = paste(batchDF$cellType, batchDF$hpf, sep = ".")) +
#	scale_y_continuous( name = "Median expression") +
#	ggtitle("Median expression in batches, \ngene and cells aggregated, \nquantile normalization" )

#png(paste0(plotDir, .Platform$file.sep,"batchMedians.png"), width = 960)
#	(batchMedianPlot)
#dev.off()

#batchVecs <- lapply(batchList, function(x) log(unlist(x)))	# legacy variable, probably not needed

batchGeneEx 		<- gather( qualDF[c(goodGenes, "hpf", "FileName", "CellType")], key = "Gene", value = "Expression", goodGenes)  

#we will use CellType.hpf levels to lable columns
cellHpfList		<- split( batchGeneEx, list(batchGeneEx$FileName, batchGeneEx$hpf), drop = TRUE)
batchGeneEx	 	<- do.call( rbind, lapply( cellHpfList, 
				function(x) data.frame(x, 
					cellHpf = paste(x$CellType, x$hpf, sep = "."),
					fileHpf = paste(x$FileName, x$hpf, sep = ".") )))
fileHpfLevels		<- levels( batchGeneEx$fileHpf)   
fileHpfLabels		<- sapply(fileHpfLevels, function(x) unique(batchGeneEx[ batchGeneEx$fileHpf == x, "cellHpf"]))
fileHpfLabels		<- factor( fileHpfLabels, levels = sort(levels( fileHpfLabels)))
fileHpfLabels		<- sort(fileHpfLabels)
#fileHpfLabels		<- fileHpfLabels[ sort(names(fileHpfLabels))]
batchGeneEx$fileHpf	<- factor( batchGeneEx$fileHpf, levels = names(fileHpfLabels))
fileHpfTags		<- as.character( fileHpfLabels)
 
batchBoxPlotNotNormalized <- ggplot( data = batchGeneEx, aes( x = fileHpf, y = log10(Expression))) +
	geom_boxplot( fill = "blue") +
	theme( 	panel.background = element_rect( fill = "gray80"),
	title 		= element_text( size = 17),
	axis.title = element_text( size = 15),
	axis.text.x = element_text(
		angle = 90,
		hjust = 1,
		vjust = 0.5,
		size = 11), 
	axis.text.y = element_text( size = 13),
	plot.title  = element_text( hjust = 0)
	) +
	geom_hline( yintercept = log10(normMedianThreshold), col = "red", linetype = "longdash") +
	scale_x_discrete( name = "batch", labels = fileHpfTags) +
	scale_y_continuous( name = "log10 expression", expand = c(0,0)) +
	ggtitle( "log-expression in batches, \nnot normalized")

batchQuantNormGeneEx <- gather( qualNormDF[c(goodGenes, "hpf", "FileName", "CellType")], key = "batchFile", value = "Expression", goodGenes)  

cellHpfList		<- split( batchQuantNormGeneEx, list(batchQuantNormGeneEx$FileName, batchQuantNormGeneEx$hpf), drop = TRUE)
batchQuantNormGeneEx	 	<- do.call( rbind, lapply( cellHpfList, 
				function(x) data.frame(x, 
					cellHpf = paste(x$CellType, x$hpf, sep = "."),
					fileHpf = paste(x$FileName, x$hpf, sep = ".") )))
fileHpfLevels		<- levels( batchQuantNormGeneEx$fileHpf)   
fileHpfLabels		<- sapply(fileHpfLevels, function(x) unique(batchQuantNormGeneEx[ batchQuantNormGeneEx$fileHpf == x, "cellHpf"]))
fileHpfLabels		<- factor( fileHpfLabels, levels = sort(levels( fileHpfLabels)))
fileHpfLabels		<- sort(fileHpfLabels)
#fileHpfLabels		<- fileHpfLabels[ sort(names(fileHpfLabels))]
batchQuantNormGeneEx$fileHpf	<- factor( batchQuantNormGeneEx$fileHpf, levels = names(fileHpfLabels))
fileHpfTags		<- as.character( fileHpfLabels)

batchBoxPlotQuantileNormalized <- ggplot( data = batchQuantNormGeneEx, aes( x = fileHpf, y = log10(Expression))) +
	geom_boxplot( fill = "blue") +
	theme( 	panel.background = element_rect( fill = "gray80"),
	title 		= element_text( size = 17),
	axis.title = element_text( size = 15),
	axis.text.x = element_text(
		angle = 90,
		hjust = 1,
		vjust = 0.5,
		size = 11), 
	axis.text.y = element_text( size = 13),
	plot.title  = element_text( hjust = 0)
	) +
	geom_hline( yintercept = log10(normMedianThreshold), col = "red", linetype = "longdash") +
	scale_x_discrete( name = "batch", 
			labels = fileHpfTags) +
	scale_y_continuous( name = "log10 expression", expand = c(0,0)) +
	ggtitle( "log expression in batches, \nquantile normalized")

batchAgrMedian <- aggregate( list( Expression = batchQuantNormGeneEx$Expression), by = list( batchMedian = batchQuantNormGeneEx$fileHpf), FUN = median)
batchAgrMedian$hpfTag 	<- fileHpfTags 

batchMedianPlotQuantileNormalized <- ggplot( data = batchAgrMedian, aes( x = batchMedian, y = Expression)) +
	geom_col( fill = "blue") +
	theme( 	panel.background = element_rect( fill = "gray80"),
	title 		= element_text( size = 17),
	axis.title = element_text( size = 15),
	axis.text.x = element_text(
		angle = 90,
		hjust = 1,
		vjust = 0.5,
		size = 11), 
	axis.text.y = element_text( size = 13),
	plot.title  = element_text( hjust = 0)
	) +
	geom_hline( yintercept = normMedianThreshold, col = "red", linetype = "longdash") +
	scale_x_discrete( name = "batch", 
			labels = fileHpfTags) +
	scale_y_continuous( name = "Median expression", expand = c(0,0)) +
	ggtitle( "Median expression in batches, \nquantile normalized")

batchProbl 	<- batchDF[ batchDF$normMedian < normMedianThreshold, ]	#Problematic batches; threshold to keep Tl (MedNormVal ~ 90)
cellsProbl_Ind	<- which( qualDF$FileName %in% batchProbl$fileName & qualDF$hpf %in% batchProbl$hpf)

qualDF_f	<- qualDF[ -cellsProbl_Ind, ]
Genes_f		<- t(qualDF_f[ , goodGenes])

cat(file = qualContLogFile, nrow( qualDF_f), " cells remaining after removing batches with low medians\n")

#
#png(paste0(plotDir, .Platform$file.sep, "LogSortPosCoefs.png"))
#	plot(sortLogCoefs, main = "Positive quality control coef values")
#	abline( h = 1,   col = "green")
#	abline( h = 1.1, col = "red")
#	abline( h = 0.9, col = "red")
#dev.off()

#keepCoefs <- which(coefs > 0.9 & coefs < 1.1)
#
#Genes_ff	<- Genes_f[,keepCoefs]
#Cells_ff	<- Cells_f[,keepCoefs]
#
#cat(file = qualContLogFile, ncol(Genes_ff), " cells remaining after removing cells with poor positive control values\n")
#

KanamycinPosTopQuant	<- quantile( log10(qualDF_f$Kanamycin.Pos), 0.995)
KanamycinPosBotQuant	<- quantile( log10(qualDF_f$Kanamycin.Pos), 0.005)

kanamycinExpressionDistributionPlot <- ggplot( data = qualDF_f, aes( x = log10(Kanamycin.Pos))) +
	geom_histogram( fill = "blue", binwidth = 0.05) +
	geom_vline( xintercept = KanamycinPosTopQuant, col = "red", linetype = "longdash") +
	geom_vline( xintercept = KanamycinPosBotQuant, col = "red", linetype = "longdash") +
	geom_vline( xintercept = min( log10( qualDF_f$Kanamycin.Pos[ grep("M", rownames(qualDF_f))])), col = "black", linetype = "longdash") +
	geom_vline( xintercept = max( log10( qualDF_f$Kanamycin.Pos[ grep("M", rownames(qualDF_f))])), col = "black", linetype = "longdash") +
	geom_vline( xintercept = min( log10( qualDF_f$Kanamycin.Pos[ grep("I", rownames(qualDF_f))])), col = "cyan", linetype = "longdash") +
	geom_vline( xintercept = max( log10( qualDF_f$Kanamycin.Pos[ grep("I", rownames(qualDF_f))])), col = "cyan", linetype = "longdash") +
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 0, 
		hjust = 1, 
		vjust = 0.5, 
		size = 12),
		axis.text.y = element_text( size = 12)) +
	scale_y_continuous( name = "Counts", expand = c(0,0)) +
	scale_x_discrete( name = "log10 Kanamycin expression", expand = c(0,0), limits = c(0, 7) ) +
	ggtitle( "Kanamycin expresspion") 


#
#png( file.path( plotDir, "Kanamycin_Distr_plot.png"))
#	plot(sort( KanamycinPos_Distr))
#	abline( h = quantile( KanamycinPos_Distr, KanamycinPosTopQuant), col = "red")
#	abline( h = quantile( KanamycinPos_Distr, KanamycinPosBotQuant), col = "red")
#	abline( h = min( KanamycinPos_Distr[ grep("M", names( KanamycinPos_Distr))]), col = "black")
#	abline( h = max( KanamycinPos_Distr[ grep("M", names( KanamycinPos_Distr))]), col = "black")
#	abline( h = max( KanamycinPos_Distr[ grep("I", names( KanamycinPos_Distr))]), col = "cyan")
#	abline( h = min( KanamycinPos_Distr[ grep("I", names( KanamycinPos_Distr))]), col = "cyan")
#dev.off()


KanamycinPos_keep	<- which( log10(qualDF_f$Kanamycin.Pos) > KanamycinPosBotQuant  
				& log10(qualDF_f$Kanamycin.Pos) < KanamycinPosTopQuant)

qualDF_ff	<- qualDF_f[ KanamycinPos_keep, ]
Genes_ff	<- t(qualDF_f[ KanamycinPos_keep, goodGenes])

cat(file = qualContLogFile, nrow( qualDF_ff), " cells remaining after removing cells with Kanamycin Pos dropouts\n")

rpl13TopQuant	<- quantile( log10(qualDF_ff$rpl13), 0.97)
rpl13BotQuant	<- quantile( log10(qualDF_ff$rpl13), 0.03)

rpl13ExpressionDistributionPlot <- ggplot( data = qualDF_ff, aes( x = log10(rpl13))) +
	geom_histogram( fill = "blue", binwidth = 0.05) +
	geom_vline( xintercept = rpl13TopQuant, col = "red", linetype = "longdash") +
	geom_vline( xintercept = rpl13BotQuant, col = "red", linetype = "longdash") +
	geom_vline( xintercept = min( log10( qualDF_ff$rpl13[ grep("M", rownames(qualDF_ff))])), col = "black", linetype = "longdash") +
	geom_vline( xintercept = max( log10( qualDF_ff$rpl13[ grep("M", rownames(qualDF_ff))])), col = "black", linetype = "longdash") +
	geom_vline( xintercept = min( log10( qualDF_ff$rpl13[ grep("I", rownames(qualDF_ff))])), col = "cyan", linetype = "longdash") +
	geom_vline( xintercept = max( log10( qualDF_ff$rpl13[ grep("I", rownames(qualDF_ff))])), col = "cyan", linetype = "longdash") +
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 12), 
		axis.text.x = element_text(angle = 0, 
		hjust = 1, 
		vjust = 0.5, 
		size = 12),
		plot.title  = element_text( hjust = 0),
		axis.text.y = element_text( size = 12)) +
	scale_y_continuous( name = "Counts", expand = c(0,0)) +
	scale_x_discrete( name = "log10 Rpl13 Expression", expand = c(0,0), limits = c(0, 7) ) +
	ggtitle( "Rpl13 expression") 


#
#png( file.path( plotDir, "Kanamycin_Distr_plot.png"))
#	plot(sort( rpl13_Distr))
#	abline( h = quantile( rpl13_Distr, rpl13TopQuant), col = "red")
#	abline( h = quantile( rpl13_Distr, rpl13BotQuant), col = "red")
#	abline( h = min( rpl13_Distr[ grep("M", names( rpl13_Distr))]), col = "black")
#	abline( h = max( rpl13_Distr[ grep("M", names( rpl13_Distr))]), col = "black")
#	abline( h = max( rpl13_Distr[ grep("I", names( rpl13_Distr))]), col = "cyan")
#	abline( h = min( rpl13_Distr[ grep("I", names( rpl13_Distr))]), col = "cyan")
#dev.off()


rpl13_keep	<- which( log10(qualDF_ff$rpl13) > rpl13BotQuant  
				& log10(qualDF_ff$rpl13) < rpl13TopQuant)

qualDF_fff	<- qualDF_ff[ rpl13_keep, ]

#qualDF_fff[ , "Kanamycin.Pos"] <- qualDF_fff[ 1, "Kanamycin.Pos"]
#qualDF_fff[ , "rpl13"]	       <- qualDF_fff[ 1, "rpl13"]
Genes_fff	<- t(qualDF_ff[ rpl13_keep, goodGenes])



#rpl13_Distr	<- log(as.numeric(Genes_fff["rpl13", ]))
#names(rpl13_Distr) <- colnames(Genes_fff)
#
#rpl13TopQuant	<- 0.97
#rpl13BotQuant	<- 0.03
#
#rpl13_keep	<- which( rpl13_Distr > quantile(rpl13_Distr, rpl13BotQuant) & rpl13_Distr < quantile(rpl13_Distr, rpl13TopQuant))
#
#png( file.path( plotDir, "rpl13_Distr_plot.png"))
#	plot(sort(rpl13_Distr))
#	abline( h = quantile( rpl13_Distr, rpl13TopQuant), col = "red")
#	abline( h = quantile( rpl13_Distr, rpl13BotQuant), col = "red")
#	abline( h = min( rpl13_Distr[ grep("M", names(rpl13_Distr))]), col = "black")
#	abline( h = max( rpl13_Distr[ grep("M", names(rpl13_Distr))]), col = "black")
#	abline( h = max( rpl13_Distr[ grep("I", names(rpl13_Distr))]), col = "cyan")
#	abline( h = min( rpl13_Distr[ grep("I", names(rpl13_Distr))]), col = "cyan")
#dev.off()

#Genes_ffff	<- Genes_fff[,rpl13_keep]
#Cells_ffff	<- Cells_fff[,rpl13_keep]

cat(file = qualContLogFile, ncol(Genes_fff), " cells remaining after removing cells with rpl13 dropouts\n")

normIteration <- 0
repeat{									       #iterations over the background level								      
normIteration <- normIteration + 1
cat("Starting iteration ", normIteration, "\n")

#NanoTable	<- cbind(Probes[,"Class.Name"], Probes[,"Gene.Name"], Probes[,"Accession.."], Genes_fff, stringsAsFactors = FALSE)
NanoTable	<- data.frame( Probes[, c("Class.Name", "Gene.Name", "Accession..")], t(qualDF_fff[ , goodGenes]), stringsAsFactors = FALSE)

colnames(NanoTable)[1:3] <- c("Code.Class", "Name", "Accession")
rownames(NanoTable) 	<- NanoTable$Name

myCodeCount 	<- "geo.mean"
myBackground	<- "mean"
mySampleContent	<- "housekeeping.geo.mean"

NanoTable[c("Kanamycin.Pos", "rpl13"), "Code.Class"] <- "Housekeeping"
NanoTableNormed <- NanoStringNorm(x = NanoTable, CodeCount = myCodeCount, Background = myBackground, SampleContent = mySampleContent)

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

cat( "Cells to drop ->" , length(norm_drop), "\n")

if( length(norm_drop) == 0){break}

#Genes_ffff	<- Genes_ffff[, -norm_drop]
#Cells_ffff	<- Cells_ffff[, -norm_drop]

qualDF_fff	<- qualDF_fff[ -norm_drop , ]
Genes_fff	<- t(qualDF_fff[ , goodGenes])

}										# repeat

cat( file = qualContLogFile, ncol(Genes_fff), " cells remaining after removing cells with poor normalization statistics\n")

#we assigned the normalized value as a fresh round of normalization only on cells which survived the previous round

NanoTable	<- data.frame( Probes[, c("Class.Name", "Gene.Name", "Accession..")], t(qualDF_fff[ , goodGenes]), stringsAsFactors = FALSE)

colnames(NanoTable)[1:3] <- c("Code.Class", "Name", "Accession")
rownames(NanoTable) 	<- NanoTable$Name

NanoTable[c("Kanamycin.Pos", "rpl13"), "Code.Class"] <- "Housekeeping"
Genes_n <- NanoStringNorm(x = NanoTable, CodeCount = myCodeCount, Background = myBackground, SampleContent = mySampleContent, return.matrix.of.endogenous.probes = TRUE)
Genes_n	<- Genes_n + 1

goodCells	<- colnames( Genes_n)
signGenes	<- rownames( Genes_n)

qualNanoNormDF 	<- data.frame( t(Genes_n), qualDF_fff[ goodCells, c("num", "batch", "hpf", "dateEx", "Desk", "CellType", "FileName", "negProbSum", "posProbCoef")])  
batchNanoNormGeneEx <- gather( qualNanoNormDF[  , c( signGenes, "hpf", "FileName", "CellType")], key = "Gene", value = "Expression", signGenes)  

cellHpfList		<- split( batchNanoNormGeneEx, list(batchNanoNormGeneEx$FileName, batchNanoNormGeneEx$hpf), drop = TRUE)
batchNanoNormGeneEx	 	<- do.call( rbind, lapply( cellHpfList, 
				function(x) data.frame(x, 
					cellHpf = paste(x$CellType, x$hpf, sep = "."),
					fileHpf = paste(x$FileName, x$hpf, sep = ".") )))
fileHpfLevels		<- levels( batchNanoNormGeneEx$fileHpf)   
fileHpfLabels		<- sapply(fileHpfLevels, function(x) unique(batchNanoNormGeneEx[ batchNanoNormGeneEx$fileHpf == x, "cellHpf"]))
fileHpfLabels		<- factor( fileHpfLabels, levels = sort(levels( fileHpfLabels)))
fileHpfLabels		<- sort(fileHpfLabels)
#fileHpfLabels		<- fileHpfLabels[ sort(names(fileHpfLabels))]
batchNanoNormGeneEx$fileHpf	<- factor( batchNanoNormGeneEx$fileHpf, levels = names(fileHpfLabels))
fileHpfTags		<- as.character( fileHpfLabels)

boxPlotNanoStringNormalized <- ggplot( data = batchNanoNormGeneEx, aes( x = fileHpf, y = log10(Expression))) +
	geom_boxplot( fill = "blue") +
	theme( 	panel.background = element_rect( fill = "gray80"),
	title 		= element_text( size = 17),
	axis.title = element_text( size = 15),
	axis.text.x = element_text(
		angle = 90,
		hjust = 1,
		vjust = 0.5,
		size = 11), 
	axis.text.y = element_text( size = 13), 
	plot.title  = element_text( hjust = 0)
	) +
	geom_hline( yintercept = log10(normMedianThreshold), col = "red", linetype = "longdash") +
	scale_x_discrete( name = "batch", 
		labels = fileHpfTags) +
	scale_y_continuous( name = "Gene expression", expand = c(0,0)) +
	ggtitle( "log-expression in batches,\nNanoString normalization")


write.table(Genes_n, file = file.path(resQCDir, "NormalizedExTable.csv"), sep = "\t")

cells_tbl <- table(qualNanoNormDF$CellType)

write.table(cells_tbl, file=qualContLogFile, sep = "\t")

close(qualContLogFile)	

write.table( Genes_n, file = file.path( resQCDir, "New_expressionTableDedupQC.csv"), sep = "\t" )
write.table( qualDF_fff[ , c("num", "batch", "hpf", "dateEx", "Desk", "CellType", "FileName", "negProbSum", "posProbCoef")] , file = file.path( resQCDir, "New_cellDescripitonsDedupQC.csv"),sep = "\t" )

#and now make the final figure

firstLine	<- plot_grid(
			negProbDistrPlot 			+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			posProbDistrPlot 			+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			geneTop5CellPlot 			+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
#			+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			nrow = 1,
			align = "h",
			labels = c("A", "B", "C"),
			label_size = 25,
			rel_widths = c( 1, 1)
			)

secondLine	<- plot_grid( 
			densPlot				+ theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
			kanamycinExpressionDistributionPlot 	+ theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
			rpl13ExpressionDistributionPlot		+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			geneAverageDistributionPlot 		+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			nrow = 1,
			align = "h",
			labels = c("D", "E", "F", "G"),
			label_size = 25,
			rel_widths = c( 1, 1)
			)

thirdLine	<- plot_grid(
			batchMedianPlotQuantileNormalized 	+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			batchBoxPlotNotNormalized		+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			nrow = 1,
			align = "h",
			labels = c("H", "I"),
			label_size = 25,
			rel_widths = c( 1, 1)
			)

forthLine	<- plot_grid(
			batchBoxPlotQuantileNormalized		+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			boxPlotNanoStringNormalized		+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			nrow = 1,
			align = "h",
			labels = c("J", "K"),
			label_size = 25,
			rel_widths = c( 1, 1)
			)



qualityControlPlot 	<- plot_grid( 
				firstLine  + theme( plot.margin = unit( c( 0.3, 0, 0.3, 1), "inches")),  
				secondLine + theme( plot.margin = unit( c( 0.3, 0, 0.3, 1), "inches")),  
				thirdLine + theme( plot.margin = unit( c( 0.3, 0, 0.3, 1), "inches")),  
				forthLine + theme( plot.margin = unit( c( 0.3, 0, 0., 1), "inches")),  
					align = "v",
					labels = '',
					rel_heights = c(1, 1),
					ncol = 1) +
			   theme( plot.margin = unit( c( 1, 1, 1, 1), "inches"))

png( file.path( plotDir, paste0("qualityControlMainPlot", ".png")), width = 1536, height = 2048)
	plot( qualityControlPlot)
dev.off()









 


		
 

