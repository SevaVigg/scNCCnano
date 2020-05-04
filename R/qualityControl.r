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

#remove mitfa 

qualDF		<- data.frame( t(Genes), t(Cells))
qualDF$Cell	<- gsub( ".", "_", rownames(qualDF), fixed = TRUE)

#Remove mitfa-/- mutants

qualDF		<- qualDF[ -grep( "mitfa", qualDF$Cell), ]

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
nCellNegGrob 		<- grobTree( textGrob( paste0(nCellNeg, " cells remaining"), x=0.55,  y=0.95, hjust=0,
  				gp=gpar(col="black", fontsize = 26, fontface="bold")))

negProbDistrPlot <- ggplot(data = negProbDF, aes( x = negProbVal)) +
	geom_vline( xintercept = negThreshold, col = "red", linetype = "solid") +
	geom_histogram( fill = "deepskyblue2", binwidth = 2) +
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.title 	= element_text( size = 30),
		axis.text 	= element_text( size = 24), 
		axis.text.x 	= element_text( angle = 0, hjust = 1, vjust = 0.5),
		axis.text.y	= element_text( margin = margin( l = 10)),
		plot.title  	= element_text( size = 34, hjust = 0)) +
	scale_y_continuous( name = "#Probes", breaks = seq(1000, 5000, by = 1000), expand = c(0,0)) +
	scale_x_continuous( name = "Negative probe counts", expand = c(0,0) ) +
	annotation_custom( nCellNegGrob) +
	ggtitle( "Negative probes") 

qualDF		 	<- qualDF[ negProbCellIndex, ]

#now check positive probes
posSpikes		<- c(128, 32, 8, 2, 0.5, 0.125)						# Pos probes in fM
posProbes		<- grep("POS", colnames( qualDF), value= TRUE)
qualDF$posProbCoef 	<- apply( log2( qualDF[ , posProbes ]), 1, function(x) lm( x~log2(posSpikes))$coefficients[2] )
leftPosProbThrsld	<- quantile( qualDF$posProbCoef, 0.01)
rightPosProbThrsld	<- quantile( qualDF$posProbCoef, 0.99)

posProbCellIndex 	<- which(qualDF$posProbCoef > leftPosProbThrsld & qualDF$posProbCoef < rightPosProbThrsld)
nCellPos		<- length( posProbCellIndex)
nCellPosGrob 		<- grobTree( textGrob( paste0(nCellPos, " cells remaining"), x=0.05,  y=0.95, hjust=0,
  				gp=gpar(col="black", fontsize = 26, fontface="bold")))

posProbDistrPlot <- ggplot( data = qualDF[ , "posProbCoef", drop = FALSE], aes( x = posProbCoef)) +
	geom_histogram( aes(y = ..count..), fill = "deepskyblue2", binwidth = 0.01) +
	geom_vline( xintercept = leftPosProbThrsld, col = "red", linetype = "solid") +
	geom_vline( xintercept = rightPosProbThrsld, col = "red", linetype = "solid") +
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.title 	= element_text( size = 30),  
		axis.text 	= element_text( size = 24), 
		axis.text.x 	= element_text( angle = 0, hjust = 1, vjust = 0.5),
		axis.text.y	= element_text( margin = margin( l = 10)), 
		plot.title  	= element_text( size = 34, hjust = 0)) +
	scale_y_continuous( name = "Counts", expand = c(0,0)) +
	scale_x_continuous( name = "regression coef. of log positive probe counts", expand = c(0,0), breaks = seq(0, 1.25, 0.25)) +
	coord_cartesian(xlim = c(0, 1.25)) +  
	annotation_custom( nCellPosGrob) +
	ggtitle( "Positive probes") 
				

qualDF		 <- qualDF[ posProbCellIndex, ]

#We check agreement between external (Kanamycin) and internal (rpl13) controls. We do this by making a linear model and looking at residuals

lmKanRpl		<- lm( qualDF$Kanamycin.Pos ~ qualDF$rpl13 )
qualDF$KanRplRes	<- lmKanRpl$residuals
leftKanRplThrsld	<- quantile( qualDF$KanRplRes, 0.025)
rightKanRplThrsld	<- quantile( qualDF$KanRplRes, 0.975)

kanRplIndex 		<- which(qualDF$KanRplRes > leftKanRplThrsld & qualDF$KanRplRes < rightKanRplThrsld)
nCellKanRpl		<- length( kanRplIndex)
nCellKanRplGrob 	<- grobTree( textGrob( paste0(nCellKanRpl, " cells remaining"), x=0.05,  y=0.95, hjust=0,
  				gp=gpar(col="black", fontsize = 26, fontface="bold")))

kanRplDistrPlot <- ggplot( data = qualDF[ , "KanRplRes", drop = FALSE], aes( x = KanRplRes)) +
	geom_histogram( aes(y = ..count..), fill = "deepskyblue2", bins = 50) +
	geom_vline( xintercept = leftKanRplThrsld, col = "red", linetype = "solid") +
	geom_vline( xintercept = rightKanRplThrsld, col = "red", linetype = "solid") +
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.title 	= element_text( size = 30),  
		axis.text 	= element_text( size = 24), 
		axis.text.x 	= element_text( angle = 0, hjust = 1, vjust = 0.5),
		axis.text.y	= element_text( margin = margin( l = 10)), 
		plot.title  	= element_text( size = 34, hjust = 0)) +
	scale_y_continuous( name = "Counts", expand = c(0,0)) +
	scale_x_continuous( name = "regression coef. of Kanamycin and rpl13 probes", expand = c(0,0) ) +
	coord_cartesian(xlim = c( min(qualDF$KanRplRes), max(qualDF$KanRplRes))) +  
	annotation_custom( nCellKanRplGrob) +
	ggtitle( "Kanamycin rpl13 agreement") 

qualDF	<- qualDF[ kanRplIndex, ]
	
exoGenes 	<- setdiff( geneNames, c( "Kanamycin.Pos", "rpl13", grep("(NEG_|POS_)", geneNames, value = TRUE)))
log10Exps	<- log10( qualDF[ ,exoGenes])	

dens		<- density(as.matrix(log10Exps)) 
expLogMinimal	<- optimize(approxfun(dens$x,dens$y),interval=c(2,4))$minimum
geneLogThreshold <- expLogMinimal
genesPerCell	<- 3	 
poorCellsDens 	<- which( apply(log10Exps, 1, function(x) sum(x > geneLogThreshold) <= genesPerCell ))  #cells with poor values for all genes but houskeeping
qualDF		<- qualDF[ -poorCellsDens, ]
nCells		<- nrow( qualDF)
nCellsExpGrob 	<- grobTree( textGrob( paste0("At least ", genesPerCell, " genes per cell \nabove the threshold\n", nCells, " cells remaining"), x=0.45,  y=0.9, hjust=0,
  				gp=gpar(col="black", hjust = "bottom", fontsize = 26, fontface="bold")))


densPlot	<- ggplot( data = data.frame( dens$x, dens$y), mapping = aes( dens$x, dens$y)) + 
			geom_line() +
			geom_vline( xintercept = geneLogThreshold, color = "red") +
			scale_x_continuous( name = "log10 probe counts" ) +
			scale_y_continuous( name = "Probe counts density") +
			annotation_custom( nCellsExpGrob) +
			theme(	plot.title 	= element_text( size = 34, hjust = 0),
				axis.title	= element_text( size = 30), 
				axis.text	= element_text( size = 24),
				axis.text.y	= element_text( margin = margin( l = 10))
				) +
			ggtitle( "Exogenous probes")


#remove genes with very low expression values in all cells
geneAverageData 		<- data.frame( gene <- factor(rownames(Genes), levels = rownames(Genes)), avExp = apply( qualDF[ , rownames(Genes) ], 2, mean))
geneAverageData$top5		<- apply( qualDF[ , rownames(Genes) ] , 2, function(x) head(tail(sort(unlist(x)), 5), 1))	# smallest of the top five

cellsPerGene	<- 5 
poorGenes 	<- colnames(log10Exps)[which( apply(log10Exps, 2, function(x) sum(x > expLogMinimal) <= cellsPerGene))]  #genes poorly expressed in all cell types
genes2exclude	<- character(0)
#genes2exclude	<- c("csf1r", "sox5", "dpf3", "ets1a", "fgfr3_v2", "mycl1a", "smad9", "pax3_v2", "hbp1")
poorGenes	<- c(poorGenes, genes2exclude)
goodGenes	<- setdiff( geneNames, poorGenes)
exoGenes	<- setdiff( exoGenes, poorGenes)

geneAverageData$filterColor	<- sapply( rownames(Genes), function(x) if( x %in% poorGenes) "red" else "deepskyblue2")


nGenesGrob 			<- grobTree( textGrob( paste0( "Gene scoring in at least ", cellsPerGene, " cells\n", length( exoGenes), " genes remainng"), 
				x=0.05,  y=0.9, hjust=0, gp=gpar(col="black", fontsize = 26, fontface="bold")))


geneTop5CellPlot <- ggplot( data = geneAverageData) + 
	aes( x = gene, y = log10( top5 )) +
	geom_col(fill = geneAverageData$filterColor) +
	geom_hline( yintercept = expLogMinimal, color = "red") + 
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.title	= element_text( size = 30), 
		axis.text 	= element_text( size = 24), 
		axis.text.x 	= element_text( size = 16, angle = 90, hjust = 1, vjust = 0.5),
		axis.text.y	= element_text( margin = margin( l = 10)), 
		plot.title  	= element_text( size = 34, hjust = 0)) +
	scale_x_discrete( name = "Probe" ) +
	scale_y_continuous( name = "log10 probe count", expand = c(0,0), limits = c(0, 7)) +
	annotation_custom( nGenesGrob) +
	ggtitle( "top 5th probe count")

#now we remove cells with a very low expression for all genes.


#make Probe value plot and Probe value 
geneAverageDistributionPlot <- ggplot( data = geneAverageData) + 
	aes( x = gene, y = avExp) +
	geom_col( fill = geneAverageData$filterColor) +
	geom_hline( yintercept = 10^geneLogThreshold, color = "red") + 
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 24), 
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
	geom_hline( yintercept = geneLogThreshold, color = "red") + 
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 24), 
		axis.text.x = element_text(angle = 90, 
		hjust = 1, 
		vjust = 0.5, 
		size = 16),
		plot.title  = element_text( hjust = 0)) +
	scale_x_discrete( name = "Gene" ) +
	scale_y_continuous( name = "Log 10 probe counts", expand = c(0,0)) +
	ggtitle( "log10 probe counts average over cells")

geneTop10CellMedianPlot <- ggplot( data = geneAverageData) + 
	aes( x = gene, y = log10( top10median)) +
	geom_col(fill = geneAverageData$filterColor) +
	geom_hline( yintercept = geneLogThreshold, color = "red") + 
	theme(	panel.background = element_rect( fill = "gray80"),
		axis.text = element_text(size = 24), 
		axis.text.x = element_text(angle = 90, 
		hjust = 1, 
		vjust = 0.5, 
		size = 7),
		plot.title  = element_text( hjust = 0)) +
	scale_x_discrete( name = "Gene" ) +
	scale_y_continuous( name = "Expression", expand = c(0,0)) +
	ggtitle( "Median log10 probe counts in top 10 cells")

#Now exclude genes that are found in small number of cells
qualDF		<- qualDF[ , setdiff( colnames(qualDF), poorGenes) ]
Probes		<- Probes[ Probes[ ,"Gene.Name"] %in% goodGenes, ]

#calculate quantile normalization
geneMatrix	<- t(as.matrix( qualDF[ , exoGenes]))
normGeneMatrix	<- normalize.quantiles(geneMatrix)

rownames(normGeneMatrix) <- rownames( geneMatrix)
colnames(normGeneMatrix) <- colnames( geneMatrix)


qualNormDF	<- qualDF; qualNormDF[ , exoGenes] <- t(normGeneMatrix)

#Now we compare batches

batchGeneEx 		<- gather( qualDF[c(exoGenes, "hpf", "FileName", "CellType", "Cell")], key = "Gene", value = "Expression", exoGenes)
batchGeneEx		<- batchGeneEx[ batchGeneEx$Expression > 1, ]
batchGeneEx$cellHpf	<- paste( batchGeneEx$CellType, batchGeneEx$hpf, sep = ".") 
batchGeneEx$fileHpf	<- paste( batchGeneEx$FileName, batchGeneEx$hpf, sep = ".") 
batchGeneEx$fileHpf	<- factor( batchGeneEx$fileHpf, levels = unique( batchGeneEx[ order(batchGeneEx$cellHpf, batchGeneEx$fileHpf), "fileHpf"] ))

batchQuantNormGeneEx 	<- gather( qualNormDF[c(exoGenes, "hpf", "FileName", "CellType", "Cell")], key = "Gene", value = "Expression", exoGenes) 
batchQuantNormGeneEx	<- batchQuantNormGeneEx[ batchQuantNormGeneEx$Expression > 1, ]
batchQuantNormGeneEx$cellHpf	<- paste( batchQuantNormGeneEx$CellType, batchQuantNormGeneEx$hpf, sep = ".") 
batchQuantNormGeneEx$fileHpf	<- paste( batchQuantNormGeneEx$FileName, batchQuantNormGeneEx$hpf, sep = ".")
batchQuantNormGeneEx$fileHpf	<- factor( batchQuantNormGeneEx$fileHpf, 
					levels = unique(batchQuantNormGeneEx[ order(batchQuantNormGeneEx$cellHpf, batchQuantNormGeneEx$fileHpf), "fileHpf"] ))

cellHpfList		<- split( batchGeneEx, list(batchGeneEx$FileName, batchGeneEx$hpf), drop = TRUE)

#batchDF is the data.frame containing information about batches
batchDF			<- aggregate(   list( medianExp = batchGeneEx$Expression), 
					by = list( 	FileName	= batchGeneEx$FileName,
							cellHpf 	= batchGeneEx$cellHpf, 
							fileHpf 	= batchGeneEx$fileHpf, 
							CellType	= batchGeneEx$CellType,
							hpf 		= batchGeneEx$hpf),
					FUN = median)

batchQuantNormMedian		<- aggregate(   list( medianExp = batchQuantNormGeneEx$Expression), 
					by = list( 	FileName	= batchQuantNormGeneEx$FileName,
							cellHpf 	= batchQuantNormGeneEx$cellHpf, 
							fileHpf 	= batchQuantNormGeneEx$fileHpf,
							CellType	= batchQuantNormGeneEx$CellType,
							hpf 		= batchQuantNormGeneEx$hpf), 
					FUN = median)
batchDF$normMedian	<- batchQuantNormMedian$medianExp

batchDF$fileHpf		<- factor(  batchDF$fileHpf, levels = unique( batchDF[ order(batchDF$cellHpf, batchDF$fileHpf), "fileHpf"] ))
batchDF			<- batchDF[ order( batchDF$fileHpf), ]


batchDF$Good		<- batchDF$normMedian >= 10^geneLogThreshold

fileHpfLevels		<- levels( batchGeneEx$fileHpf)   
fileHpfLabels		<- sapply( fileHpfLevels, function(x) unique(batchGeneEx[ batchGeneEx$fileHpf == x, "cellHpf"]))
fileHpfTags		<- as.character( fileHpfLabels)

batchBoxPlotNotNormalized <- ggplot( data = batchGeneEx, aes( x = fileHpf, y = log10(Expression))) +
	geom_boxplot( fill = sapply(batchDF$Good, function(x) if(x) "deepskyblue2" else "red")) +
	theme( 	panel.background = element_rect( fill = "gray80"),
	plot.title 	= element_text( size = 34, hjust = 0),
	axis.title 	= element_text( size = 30),
	axis.text	= element_text( size = 24), 	
	axis.text.x 	= element_text( size = 16, angle = 90, hjust = 1, vjust = 0.5), 
	axis.text.y 	= element_text( margin = margin( l = 10))
	) +
	geom_hline( yintercept = geneLogThreshold, col = "red", linetype = "solid") +
	scale_x_discrete( name = "batch", labels = fileHpfTags) +
	scale_y_continuous( name = "log10 probe counts", expand = c(0,0), limits = c(0, 7)) +
	ggtitle( "log10 probe counts in batches, not normalized")

batchBoxPlotQuantileNormalized <- ggplot( data = batchQuantNormGeneEx, aes( x = fileHpf, y = log10(Expression))) +
	geom_boxplot( fill = sapply(batchDF$Good, function(x) if(x) "deepskyblue2" else "red")) +
	theme( 	panel.background = element_rect( fill = "gray80"),
		plot.title 	= element_text( size = 34, hjust = 0),
		axis.title 	= element_text( size = 30),
		axis.text	= element_text( size = 24),
		axis.text.x 	= element_text( size = 16, angle = 90, hjust = 1, vjust = 0.5), 
		axis.text.y = element_text( margin = margin( l = 10))
	) +
	geom_hline( yintercept = geneLogThreshold, col = "red", linetype = "solid") +
	scale_x_discrete( name = "batch", labels = fileHpfTags) +
	scale_y_continuous( name = "log10 probe counts", expand = c(0,0), limits = c(0,7)) +
	ggtitle( "log10 probe counts in batches, \nquantile normalized")

batchAgr		<- aggregate( list( medianExp = batchQuantNormGeneEx$Expression), by = list( batch = batchQuantNormGeneEx$fileHpf), FUN = quantile, probs = 0.5)
batchAgrQuart		<- aggregate( list( quart3Exp = batchQuantNormGeneEx$Expression), by = list( batch = batchQuantNormGeneEx$fileHpf), FUN = quantile, probs = 0.75)
batchAgr$quart3Exp	<- batchAgrQuart$quart3Exp
batchAgr$hpfTag 	<- fileHpfTags 
cellsBatchGood		<- which( qualDF$FileName %in% batchDF[ batchDF$Good, "FileName"] & qualDF$hpf %in% batchDF[ batchDF$Good, "hpf"])
nCellsBatch		<- length( cellsBatchGood)
nCellsBatchGrob 		<- grobTree( textGrob( paste0(nCellsBatch, " cells remaining"), x=0.05,  y=0.95, hjust=0,
  				gp=gpar(col="black", fontsize = 26, fontface="bold")))

batchQuartPlotQuantileNormalized <- ggplot( data = batchAgr, aes( x = batch, y = medianExp)) +
	geom_col( fill = sapply(batchDF$Good, function(x) if(x) "deepskyblue2" else "red")) +
	theme( 	panel.background = element_rect( fill = "gray80"),
		plot.title 	= element_text( size = 34, hjust = 0),
		axis.title 	= element_text( size = 30),
		axis.text	= element_text( size = 24), 
		axis.text.x 	= element_text( size = 16, angle = 90, hjust = 1, vjust = 0.5), 
		axis.text.y 	= element_text( margin = margin( l= 10))
		) +
	geom_hline( yintercept = 10^geneLogThreshold, col = "red", linetype = "solid") +
	scale_x_discrete( name = "batch", 
			labels = fileHpfTags) +
	scale_y_continuous( name = "Median probe counts", expand = c(0,0), limits = c( 0, (max(batchAgr$medianExp) %/% 100 + 1) * 100)) +
	annotation_custom( nCellsBatchGrob) +
	ggtitle( "Median probe counts in batches, \nquantile normalized")

qualDF_f	<- qualDF[ cellsBatchGood, ]
Genes_f		<- t(qualDF_f[ , goodGenes])


#cat(file = qualContLogFile, nrow( qualDF_f), " cells remaining after removing batches with low medians\n")
#cat(file = qualContLogFile, ncol(Genes_f), " cells remaining after removing cells with rpl13 dropouts\n")


#Now we remove cells with very low Kanamycin counts
KanamycinPosTopQuant	<- quantile( qualDF_f$Kanamycin.Pos, 1)
KanamycinPosBotQuant	<- median( qualDF_f$Kanamycin.Pos)/50
goodCellsKana		<- which( qualDF_f$Kanamycin.Pos >= KanamycinPosBotQuant & qualDF_f$Kanamycin.Pos <= KanamycinPosTopQuant)
nCells			<- length( goodCellsKana)
nCellsKanaGrob 		<- grobTree( textGrob( paste0(nCells, " cells remaining"), x=0.05,  y=0.95, hjust=0,
  				gp=gpar(col="black", fontsize = 26, fontface="bold")))

kanamycinExpressionDistributionPlot <- ggplot( data = qualDF_f, aes( x = log10(Kanamycin.Pos))) +
	geom_histogram( fill = "deepskyblue2", binwidth = 0.05) +
	geom_vline( xintercept = log10(KanamycinPosTopQuant), col = "red", linetype = "solid") +
	geom_vline( xintercept = log10(KanamycinPosBotQuant), col = "red", linetype = "solid") +
	geom_vline( xintercept = min( log10( qualDF_f$Kanamycin.Pos[ grep("M", rownames(qualDF_f))])), col = "black", linetype = "longdash") +
	geom_vline( xintercept = max( log10( qualDF_f$Kanamycin.Pos[ grep("M", rownames(qualDF_f))])), col = "black", linetype = "longdash") +
	geom_vline( xintercept = min( log10( qualDF_f$Kanamycin.Pos[ grep("I", rownames(qualDF_f))])), col = "cyan", linetype = "longdash") +
	geom_vline( xintercept = max( log10( qualDF_f$Kanamycin.Pos[ grep("I", rownames(qualDF_f))])), col = "cyan", linetype = "longdash") +
	theme(	panel.background = element_rect( fill = "gray80"),
		plot.title = element_text( size = 34, hjust = 0), 
		axis.title = element_text( size = 30), 
		axis.text = element_text( size = 24), 
		axis.text.x = element_text( angle = 0, hjust = 1, vjust = 0.5),
		axis.text.y = element_text( margin = margin( l = 10))
		) +
	scale_y_continuous( name = "#cells", expand = c(0,0)) +
	scale_x_continuous( name = "log10 Kanamycin counts") +
	coord_cartesian(xlim = c(0, 7), expand = FALSE) + 
	annotation_custom( nCellsKanaGrob) +
	ggtitle( "Kanamycin counts") 

qualDF_f		<- qualDF_f[ goodCellsKana, ]

rpl13TopQuant	<- quantile( qualDF_f$rpl13, 1)
rpl13BotQuant	<- quantile( qualDF_f$rpl13, 0.05)
goodCellsRpl	<- which( qualDF_f$rpl13 >= rpl13BotQuant & qualDF_f$rpl13 <= rpl13TopQuant)
nCells		<- length( goodCellsRpl)
nCellsRplGrob 	<- grobTree( textGrob( paste0(nCells, " cells remaining"), x=0.05,  y=0.95, hjust=0,
  				gp=gpar(col="black", fontsize = 26, fontface="bold")))


rpl13ExpressionDistributionPlot <- ggplot( data = qualDF_f, aes( x = log10(rpl13))) +
	geom_histogram( fill = "deepskyblue2", binwidth = 0.05) +
	geom_vline( xintercept = log10( rpl13TopQuant), col = "red", linetype = "solid") +
	geom_vline( xintercept = log10( rpl13BotQuant), col = "red", linetype = "solid") +
	geom_vline( xintercept = min( log10( qualDF_f$rpl13[ grep("M", rownames(qualDF_f))])), col = "black", linetype = "longdash") +
	geom_vline( xintercept = max( log10( qualDF_f$rpl13[ grep("M", rownames(qualDF_f))])), col = "black", linetype = "longdash") +
	geom_vline( xintercept = min( log10( qualDF_f$rpl13[ grep("I", rownames(qualDF_f))])), col = "cyan", linetype = "longdash") +
	geom_vline( xintercept = max( log10( qualDF_f$rpl13[ grep("I", rownames(qualDF_f))])), col = "cyan", linetype = "longdash") +
	theme(	panel.background = element_rect( fill = "gray80"),
		plot.title  	= element_text( size = 34, hjust = 0),
		axis.title 	= element_text( size = 30), 
		axis.text 	= element_text( size = 24), 
		axis.text.x 	= element_text(angle = 0, hjust = 1, vjust = 0.5),
		axis.text.y 	= element_text( margin = margin( l = 10))
		) +
	scale_y_continuous( name = "#cells", expand = c(0,0)) +
	scale_x_continuous( name = "log10 rpl13 counts", expand = c(0,0) ) +
	coord_cartesian( xlim = c(0, 7), expand = FALSE) +
	annotation_custom( nCellsRplGrob) +
	ggtitle( "Rpl13 counts")

#We need to separate all NanoString tags into two sets, selecting those with expressions of 1 (possible dropouts) which we do not normalize  


normIteration <- 0
repeat{									       #iterations over the background level								      
nanoGenes	<- t(qualDF_f[ , goodGenes] - 1)  #We need zero counts and the landscape table    

normIteration <- normIteration + 1
cat("Starting iteration ", normIteration, "\n")

#NanoTable	<- cbind(Probes[,"Class.Name"], Probes[,"Gene.Name"], Probes[,"Accession.."], Genes_fff, stringsAsFactors = FALSE)
NanoTable	<- data.frame( Probes[, c("Class.Name", "Gene.Name", "Accession..")], nanoGenes , stringsAsFactors = FALSE)

colnames(NanoTable)[1:3] <- c("Code.Class", "Name", "Accession")
rownames(NanoTable) 	<- NanoTable$Name

myCodeCount 	<- "none"
myBackground	<- "mean.2sd"
mySampleContent	<- "housekeeping.sum"


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

qualDF_f	<- qualDF_f[ -norm_drop , ]
Genes_f		<- t(qualDF_f[ , goodGenes])

}	#repeat iterations

cat( file = qualContLogFile, ncol(Genes_f), " cells remaining after removing cells with poor normalization statistics\n")

#we assigned the normalized value as a fresh round of normalization only on cells which survived the previous round


nanoGenes	<- t(qualDF_f[ , goodGenes] - 1)  #We need true zero counts and the landscape table    
NanoTable	<- data.frame( Probes[, c("Class.Name", "Gene.Name", "Accession..")], nanoGenes, stringsAsFactors = FALSE)

colnames(NanoTable)[1:3] <- c("Code.Class", "Name", "Accession")
rownames(NanoTable) 	<- NanoTable$Name

NanoTable[c("Kanamycin.Pos", "rpl13"), "Code.Class"] <- "Housekeeping"
Genes_n <- NanoStringNorm(x = NanoTable, CodeCount = myCodeCount, Background = myBackground, SampleContent = mySampleContent, return.matrix.of.endogenous.probes = TRUE)

Genes_n <- trunc(Genes_n) 

goodCells	<- colnames( Genes_n)
signGenes	<- rownames( Genes_n)

qualNanoNormDF 	<- data.frame( t(Genes_n), qualDF_f[ goodCells, c("num", "batch", "hpf", "dateEx", "Desk", "CellType", "FileName", "negProbSum", "posProbCoef")])  
batchNanoNormGeneEx <- gather( qualNanoNormDF[  , c( signGenes, "hpf", "FileName", "CellType")], key = "Gene", value = "Expression", signGenes)  

batchNanoNormGeneEx$cellHpf	<- paste( batchNanoNormGeneEx$CellType, batchNanoNormGeneEx$hpf, sep = ".") 
batchNanoNormGeneEx$fileHpf	<- paste( batchNanoNormGeneEx$FileName, batchNanoNormGeneEx$hpf, sep = ".") 
batchNanoNormGeneEx$fileHpf	<- factor( batchNanoNormGeneEx$fileHpf, levels = unique( batchNanoNormGeneEx[ order(batchNanoNormGeneEx$cellHpf, batchNanoNormGeneEx$fileHpf), "fileHpf"] ))

batchDF			<- aggregate(   list( medianExp = batchNanoNormGeneEx$Expression), 
					by = list( 	FileName	= batchNanoNormGeneEx$FileName,
							cellHpf 	= batchNanoNormGeneEx$cellHpf, 
							fileHpf 	= batchNanoNormGeneEx$fileHpf, 
							CellType	= batchNanoNormGeneEx$CellType,
							hpf 		= batchNanoNormGeneEx$hpf),
					FUN = median)

fileHpfLevels		<- levels( batchNanoNormGeneEx$fileHpf)   
fileHpfLabels		<- sapply( fileHpfLevels, function(x) unique(batchNanoNormGeneEx[ batchNanoNormGeneEx$fileHpf == x, "cellHpf"]))
fileHpfTags		<- as.character( fileHpfLabels)

nCellNanoGrob 		<- grobTree( textGrob( paste0( length( goodCells), " cells remaining"), x=0.05,  y=0.95, hjust=0,
  				gp=gpar(col="black", fontsize = 26, fontface="bold")))

batchNanoNormGeneExNonZero <- batchNanoNormGeneEx[ batchNanoNormGeneEx$Expression > 0, ]

boxPlotNanoStringNormalized <- ggplot( data = batchNanoNormGeneExNonZero, aes( x = fileHpf, y = log10(Expression))) +
	geom_boxplot( fill = "deepskyblue2") +
	theme( 	panel.background = element_rect( fill = "gray80"),
		plot.title = element_text( size = 34, hjust = 0),
		axis.title = element_text( size = 30),
		axis.text.x = element_text(
			angle = 90,
			hjust = 1,
			vjust = 0.5,
			size = 16), 
		axis.text.y = element_text( size = 24) 
		) +
	geom_hline( yintercept = geneLogThreshold, col = "red", linetype = "solid") +
	scale_x_discrete( name = "batch", 
			  labels = fileHpfTags) +
	scale_y_continuous( name = "log10 probe counts", expand = c(0,0), limits = c(-0.5,7)) +
	annotation_custom( nCellNanoGrob) +
	ggtitle( "log10 probe counts in batches,\nNanoString normalization")

write.table(Genes_n, file = file.path(resQCDir, "NormalizedExTable.csv"), sep = "\t")

cells_tbl <- table(qualNanoNormDF$CellType)

write.table(cells_tbl, file=qualContLogFile, sep = "\t")

close(qualContLogFile)	

write.table( Genes_n, file = file.path( resQCDir, "New_NormalizedExTable.csv"), sep = "\t" )
write.table( t(qualDF_f[ , c("num", "batch", "hpf", "dateEx", "Desk", "CellType", "FileName", "negProbSum", "posProbCoef")]) , file = file.path( resQCDir, "New_cellDescripitonsDedupQC.csv"),sep = "\t" )

#and now make the final figure. plot_grid opens an empty window, I don't know how to suppress this

firstLine	<- plot_grid(
			negProbDistrPlot 			+ theme( plot.margin = unit( c( 0, 0.4, 0, 0.4), "inches")), 
			posProbDistrPlot 			+ theme( plot.margin = unit( c( 0, 0.4, 0, 0.4), "inches")), 
			densPlot				+ theme( plot.margin = unit( c( 0, 0.4, 0, 0.4), "inches")), 
#			+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			nrow = 1,
			align = "h",
			labels = c("A", "B", "C"),
			label_size = 40,
			rel_widths = c( 1/3, 1/3, 1/3 )
			)

secondLine	<- plot_grid( 
			geneTop5CellPlot	 		+ theme( plot.margin = unit( c( 0, 0.4, 0, 0.4), "inches")), 
			batchBoxPlotNotNormalized		+ theme( plot.margin = unit( c( 0, 0.4, 0, 0.4), "inches")), 
			nrow = 1,
			align = "h",
			labels = c("D", "E"),
			label_size = 40,
			rel_widths = c(1, 1)
			)

thirdLine	<- plot_grid(
			batchQuartPlotQuantileNormalized 	+ theme( plot.margin = unit( c( 0, 0.4, 0, 0.4), "inches")), 
			batchBoxPlotQuantileNormalized		+ theme( plot.margin = unit( c( 0, 0.4, 0, 0.4), "inches")), 
			nrow = 1,
			align = "h",
			labels = c("F", "G"),
			label_size = 40,
			rel_widths = c( 1, 1)
			)

forthLine	<- plot_grid(
			kanamycinExpressionDistributionPlot 	+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			rpl13ExpressionDistributionPlot		+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			boxPlotNanoStringNormalized		+ theme( plot.margin = unit( c( 0, 0.1, 0, 0.1), "inches")), 
			nrow = 1,
			align = "h",
			labels = c("H", "I", "J"),
			label_size = 40,
			rel_widths = c( 0.3, 0.3, 0.4)
			)



qualityControlPlot 	<- plot_grid( 
				firstLine  +  theme( plot.margin = unit( c( 0.4, 0, 0.4, 1), "inches")),  
				secondLine +  theme( plot.margin = unit( c( 0.4, 0, 0.4, 1), "inches")),  
				thirdLine  +  theme( plot.margin = unit( c( 0.4, 0, 0.4, 1), "inches")),  
				forthLine  +  theme( plot.margin = unit( c( 0.4, 0, 0.4, 1), "inches")),  
					align = "v",
					labels = '',
					rel_heights = c(1, 1),
					ncol = 1, 
					label_size = 40) +
			   theme( plot.margin = unit( c( 1, 1, 1, 1), "inches"))

png( file.path( plotDir, paste0("qualityControlMainPlot", ".png")), width = 2480, height = 3506)
	plot( qualityControlPlot)
dev.off()

#Plot separate plots 

#png( file.path( plotDir, paste0("negProbDistrPlot", ".png")), width = 2480, height = 3506)
#	(negProbDistrPlot)
#dev.off()
#
#png( file.path( plotDir, paste0("posProbDistrPlot", ".png")), width = 2480, height = 3506)
#	(posProbDistrPlot)
#dev.off()
#
#png( file.path( plotDir, paste0("densPlot", ".png")), width = 2480, height = 3506)
#	(densPlot)
#dev.off()
#
#png( file.path( plotDir, paste0("geneTop5CellPlot", ".png")), width = 2480, height = 3506)
#	(geneTop5CellPlot)
#dev.off()
#
#png( file.path( plotDir, paste0("batchBoxPlotNotNormalized", ".png")), width = 2480, height = 3506)
#	(batchBoxPlotNotNormalized)
#dev.off()
#
#png( file.path( plotDir, paste0("batchQuartPlotQuantileNormalized", ".png")), width = 2480, height = 3506)
#	(batchQuartPlotQuantileNormalized)
#dev.off()
#
#png( file.path( plotDir, paste0("batchBoxPlotQuantileNormalized", ".png")), width = 2480, height = 3506)
#	(batchBoxPlotQuantileNormalized)
#dev.off()
#
#png( file.path( plotDir, paste0("kanamycinExpressionDistributionPlot", ".png")), width = 2480, height = 3506)
#	(kanamycinExpressionDistributionPlot)
#dev.off()
#
#png( file.path( plotDir, paste0("rpl13ExpressionDistributionPlot", ".png")), width = 2480, height = 3506)
#	(rpl13ExpressionDistributionPlot)
#dev.off()
#
#png( file.path( plotDir, paste0("boxPlotNanoStringNormalized", ".png")), width = 2480, height = 3506)
#	(boxPlotNanoStringNormalized)
#dev.off()
















 


		
 

