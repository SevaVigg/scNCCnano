# Requiers seuratObj created by seuratNorm.r
#
#this snippet is use to make TSNE, UMAP and PCA plots for initial cell types and clusters in the gene space
#it requires to run first seuratNorm.r, which among others sets the appropriate experimentType ( allCells, WT, WT_sox10)
#in seuratObj@project.name
#
# Directory structure :     Res <- Plots <- geneSpace

makeGeneSpacePlotsAll <- function( seuratObj){

source("R/getClusterTypes.r")
source("R/setCellTypeColors.r")
source("R/setClusterColors.r")

library(ape)

if(!require("gridGraphics")){
install.packages("gridGraphics")
library("gridGraphics")}

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypePlotDir <- file.path(plotDir, seuratObj@project.name)
dir.create( experimentTypePlotDir, showWarnings = FALSE)

geneSpacePlotDir <- file.path( experimentTypePlotDir, "geneSpacePlots")
dir.create( geneSpacePlotDir, showWarnings = FALSE)

geneSetPlotDir	 <- file.path( geneSpacePlotDir, "allGenes")
dir.create( geneSetPlotDir, showWarnings = FALSE)

seuratObj		<- StashIdent( seuratObj, save.name = 'hpfIdent')

levels(seuratObj@ident) <- c(levels(seuratObj@ident), "R")
seuratObj@ident[ grep("regular", names(seuratObj@ident))] <- "R"
seuratObj@ident 	<- droplevels(seuratObj@ident)

seuratObj		<- StashIdent( seuratObj, save.name = 'genCellTypeIdent')

TSNESeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( geneSetPlotDir, "TSNESeed.txt"), TSNESeed, "\n")

seuratObj <- RunTSNE( seuratObj, genes.use = rownames(seuratObj@data), seed.use = TSNESeed, 
	theta = 0, eta = 100, max_iter = 300, perplexity = 20, verbose = FALSE)

UMAPSeed <- as.numeric(as.POSIXct(Sys.time()))
UMAPSeed <- 20
	cat( file = file.path( geneSetPlotDir, "UMAPSeed.txt"), TSNESeed, "\n")
seuratObj <- RunUMAP( seuratObj, genes.use = rownames(seuratObj@data), max.dim = 2, seed.use = UMAPSeed, 
	n_neighbors = 20, min_dist = 0.15, metric = "cosine")



plotList <- list()

	tsnePlot 	<- DimPlot( object = seuratObj, reduction.use = 'tsne', cols.use = setCellTypeColors( seuratObj), pt.size = 2, do.return = TRUE)
	tsnePlot 	<- tsnePlot +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			axis.title.x = element_text( margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 0,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		) +
	xlab( label = "tSNE1")+ylab( label = "tSNE2")

	umapPlot	<- DimPlot(object = seuratObj, reduction.use = 'umap',  cols.use = setCellTypeColors( seuratObj), pt.size = 2, do.return = TRUE)
	umapPlot 	<- umapPlot +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			axis.title.x = element_text( margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 0,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		)  


	noiseTol		<- log2(19)
	seuratObjDenoise	<- seuratObj
	seuratObjDenoise@data	<- apply( seuratObjDenoise@data, c(1,2), function(x) if(x>noiseTol) x else 0)

	hpfClusterTreePlot <- function() {	seuratObj <- SetAllIdent( seuratObj, id = 'hpfIdent') 
					   	seuratObj <- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = FALSE)
						par( mar = c(5,5,5,5))
					   	PlotClusterTree( seuratObj, type = "phylogram", cex = 2); nodelabels( text = "  ")}
	genClusterTreePlot <- function() {	seuratObj <- SetAllIdent( seuratObj, id = 'genCellTypeIdent')
					   	seuratObj <- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = FALSE)
						par( mar = c(5,5,5,5))
						PlotClusterTree( seuratObj, type = "phylogram", cex = 2); nodelabels( text = "  ")}
	
	seuratObjDenoise <- SetAllIdent( seuratObj, id = 'hpfIdent')
	dotPlotHpf 	 <- DotPlot(seuratObjDenoise, genes.plot = rownames(seuratObjDenoise@data), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE)
	dotPlotHpf	<- dotPlotHpf +
		theme(
			legend.position="none", 
			axis.text.y = element_text( size = 30),
			axis.text.x = element_text( size = 20, angle = 90),
			axis.title  = element_text( size = 25, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)	 
	seuratObjDenoise <- SetAllIdent( seuratObj, id = 'genCellTypeIdent')
	dotPlotGen 	<- DotPlot(seuratObjDenoise, genes.plot = rownames(seuratObjDenoise@data), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE)
	dotPlotGen	<- dotPlotGen +
		theme(
			legend.position="none", 
			axis.text.y = element_text( size = 30),
			axis.text.x = element_text( size = 20, angle = 90),
			axis.title  = element_text( size = 25, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)	
	
	firstLine	<- plot_grid( 
				tsnePlot + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				umapPlot + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				nrow = 1,
				labels = c("A", "B"),
				label_size = 25)
	secondLine	<- plot_grid( 
				hpfClusterTreePlot,
				genClusterTreePlot, 
				nrow = 1, 
				labels = c("C", "D"),
				label_size = 25,
				rel_widths = c( 1, 1))
	thirdLine	<- plot_grid( 
				dotPlotHpf + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				nrow = 1,
				labels = c("E", "F"),
				label_size = 25)
	forthLine	<- plot_grid( 
				dotPlotGen + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				nrow = 1,
				labels = c("E", "F"),
				label_size = 25)
	 
	gridPlot 	<- plot_grid( 
				firstLine  + theme( plot.margin = unit( c( 0, 0, 1, 1), "inches")),  
				secondLine + theme( plot.margin = unit( c( 0, 0, 1, 1), "inches")),  
				thirdLine  + theme( plot.margin = unit( c( 0, 0, 1, 0), "inches")),  
				forthLine  + theme( plot.margin = unit( c( 0, 0, 1, 0), "inches")),  
					labels = '',
					rel_heights = c(1, 1),
					ncol = 1) +
			   theme( plot.margin = unit( c( 1, 1, 5, 1), "inches"))


png( file.path( geneSetPlotDir, "InitCellTypesPlots.png"), width = 1536, height = 2048)
	plot( gridPlot)
dev.off()

seuratObj	<- SetAllIdent( seuratObj, id = "genCellTypeIdent")

source("R/makeInitCellTypePCAPlots.r")
pcaPlot <- makeInitCellTypePCAPlots(seuratObj, 5)	#get ggplot2 PCA diagrams with cell colors, uses its own directorial structure

png( file.path( geneSpacePlotDir, "PCAPlots.png"), width = 1536, height = 2048)
	plot( pcaPlot)
dev.off()

seuratObj 	<- SetAllIdent(seuratObj, id = "originalCellTypes")

return( seuratObj)

#remove values, that are too close to zero

}
