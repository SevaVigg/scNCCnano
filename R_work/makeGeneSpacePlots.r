# Requiers seuratObj created by seuratNorm.r
#
#this snippet is use to make TSNE, UMAP and PCA plots for initial cell types and clusters in the gene space
#it requires to run first seuratNorm.r, which among others sets the appropriate experimentType ( allCells, WT, WT_sox10)
#in seuratObj@project.name
#
# Directory structure :     Res <- Plots <- geneSpace

source("R/getClusterTypes.r")
source("R/calcTSNEGeneSpace.r")
source("R/setCellTypeColors.r")
source("R/setClusterColors.r")
source("R/calcUmapGeneSpace.r")
source("R/seuratNorm.r")

library(ape)
require(proxy)

if(!require("gridGraphics")){
install.packages("gridGraphics")
library("gridGraphics")}

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

#experimentTypePlotDir <- file.path(plotDir, seuratObj@project.name)
#dir.create( experimentTypePlotDir, showWarnings = FALSE)

geneSpacePlotDir <- file.path( plotDir, "geneSpacePlots")
dir.create( geneSpacePlotDir, showWarnings = FALSE)

seuratWT 	<- seuratNorm("WT")
seuratAll	<- seuratNorm("allCells")

seuratWT@meta.data$genCellTypeIdent 	<- factor( seuratWT@meta.data$genCellTypeIdent, levels = c( "Tl", "G", "I", "M"))
seuratWT <- SetAllIdent( seuratWT, id = "genCellTypeIdent")

initCellTypePlotDir	 <- file.path( geneSpacePlotDir, "initCellTypes")
dir.create( initCellTypePlotDir, showWarnings = FALSE) 

TSNESeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( initCellTypePlotDir, "TSNESeed.txt"), TSNESeed, "\n")

seuratWT <- calcTSNEGeneSpace( seuratWT, TSNESeed, Norm = FALSE)
#seuratWT	<- RunTSNE( seuratWT, genes.use = rownames( seuratWT@data), perplexity = 20, 


UMAPSeed <- as.numeric(as.POSIXct(Sys.time()))
#UMAPSeed <- 42

cat( file = file.path( initCellTypePlotDir, "UMAPSeed.txt"), UMAPSeed, "\n")

#Prepare tsneWT plot

	tsnePlotWT 	<- DimPlot( object = seuratWT, reduction.use = 'tsne', cols.use = setCellTypeColors( seuratWT), pt.size = 2, do.return = TRUE)
	tsnePlotWT 	<- tsnePlotWT +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 24),
			axis.title  = element_text( size = 30),
			axis.title.x = element_text( margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 10,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		) +
	xlab( label = "tSNE1")+ylab( label = "tSNE2")

#reorder values to have a cuter DotPlot

seuratAll@meta.data$genCellTypeIdent 	<- factor( seuratAll@meta.data$genCellTypeIdent, levels = c("M", "I", "sox10-", "G", "Tl"))
seuratAll	<- SetAllIdent( seuratAll, id = 'genCellTypeIdent')

seuratAll <- calcTSNEGeneSpace( seuratAll, TSNESeedAllCells, Norm = TRUE, initSeurObj = seuratWT)

#seuratAll	<- RunTSNE( seuratAll, genes.use = rownames( seuratAll@data), perplexity = 20, 
#seruatObj 	<- RunUMAP( seuratAll, genes.use = rownames( seuratAll@data), n_neighbors = 15L, min.dist = 0.01, metric = "cosine", seed.use = UMAPSeedAllCells)

umapRes 	<- calcUmapGeneSpace( seuratAll, seurWT = seuratWT, UMAPRandSeed = UMAPSeed, Dim = 2, experimentType = 'allCondWT', mySpread = 1, minDist = 0.65)	

	tsnePlotAllCells 	<- DimPlot( object = seuratAll, reduction.use = 'tsne', cols.use = setCellTypeColors( seuratAll), pt.size = 2, do.return = TRUE)
	tsnePlotAllCells 	<- tsnePlotAllCells +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 24),
			axis.title  = element_text( size = 30),
			axis.title.x = element_text( margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 0,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		) +
	xlab( label = "tSNE1")+ylab( label = "tSNE2")

	umapPlotAllWT	<- DimPlot(object = umapRes$WT, reduction.use = 'umap',  cols.use = setCellTypeColors( seuratWT), pt.size = 2, do.return = TRUE)
	umapPlotAllWT 	<- umapPlotAllWT +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 24),
			axis.title  = element_text( size = 30, face = "bold"),
			axis.title.x = element_text( margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 10,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		)

	umapPlotAllCells	<- DimPlot(object = umapRes$All, reduction.use = 'umap',  cols.use = setCellTypeColors( seuratAll), pt.size = 2, do.return = TRUE)
	umapPlotAllCells 	<- umapPlotAllCells +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 24),
			axis.title  = element_text( size = 30, face = "bold"),
			axis.title.x = element_text( margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 0,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		)  


	noiseTol		<- log2(19)
	seuratAllDenoise	<- seuratAll
	seuratAllDenoise@data	<- apply( seuratAllDenoise@data, c(1,2), function(x) if(x>noiseTol) x else 0)

	hpfClusterTreePlot <- function() {	seuratAll <- SetAllIdent( seuratAll, id = 'originalCellTypes') 
					   	seuratAll <- BuildClusterTree( seuratAll, genes.use = rownames(seuratAll@data), do.plot = FALSE, do.reorder = FALSE)
					   	PlotClusterTree( seuratAll, type = "phylogram", cex = 2); nodelabels( text = "  ")}
	genClusterTreePlot <- function() {	seuratAll <- SetAllIdent( seuratAll, id = 'genCellTypeIdent')
					   	seuratAll <- BuildClusterTree( seuratAll, genes.use = rownames(seuratAll@data), do.plot = FALSE, do.reorder = FALSE)
						PlotClusterTree( seuratAll, type = "phylogram", cex = 2); nodelabels( text = "  ")}

	dotPlot 	<- DotPlot(seuratAllDenoise, genes.plot = rownames(seuratAllDenoise@data), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE)
	dotPlot		<- dotPlot +
		theme(
			legend.position="none", 
			axis.text.y = element_text( size = 24),
			axis.text.x = element_text( size = 24, angle = 90),
			axis.title  = element_text( size = 30, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)	 

	
	firstLine	<- plot_grid( 
				tsnePlotWT + theme( plot.margin = unit( c( 0, 0.15, 0, 0.15), "inches")), 
				umapPlotAllWT + theme( plot.margin = unit( c( 0, 0.15, 0, 0.15), "inches")), 
				nrow = 1,
				labels = c("A", "B"),
				label_size = 40)
	secondLine	<- plot_grid( 
				tsnePlotAllCells + theme( plot.margin = unit( c( 0, 0.15, 0, 0.15), "inches")), 
				umapPlotAllCells + theme( plot.margin = unit( c( 0, 0.15, 0, 0.15), "inches")), 
				nrow = 1,
				labels = c("С", "D"),
				label_size = 40)
	thirdLine	<- plot_grid( 
				hpfClusterTreePlot,
				genClusterTreePlot, 
				nrow = 1, 
				labels = c("E", "F"),
				label_size = 40,
				rel_widths = c( 1, 1))
	forthLine	<- plot_grid( 
				dotPlot + theme( plot.margin = unit( c( 0, 0.15, 0, 0.15), "inches")), 
				nrow = 1,
				labels = c("E"),
				label_size = 40)
	 
	gridPlot 	<- plot_grid( 
				firstLine  + theme( plot.margin = unit( c( 1, 0, 1, 0), "inches")),  
				secondLine + theme( plot.margin = unit( c( 1, 0, 1, 0), "inches")),  
				thirdLine  + theme( plot.margin = unit( c( 1, 0, 1, 0), "inches")),  
				forthLine  + theme( plot.margin = unit( c( 1, 0, 3, 0), "inches")),  
					labels = '',
					rel_heights = c(3, 3, 2, 3),
					ncol = 1) +
			   theme( plot.margin = unit( c( 1, 1, 5, 1), "inches"))


png( file.path( initCellTypePlotDir, "InitCellTypesPlots.png"), width = 2480, height = 3506)
	plot( gridPlot)
dev.off()


