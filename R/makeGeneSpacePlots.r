# Requiers seuratObj created by seuratNorm.r
#
#this snippet is use to make TSNE, UMAP and PCA plots for initial cell types and clusters in the gene space
#it requires to run first seuratNorm.r, which among others sets the appropriate experimentType ( allCells, WT, WT_sox10)
#in seuratObj@project.name
#
# Directory structure :     Res <- Plots <- geneSpace

#makeGeneSpacePlots <- function( seuratObj){

source("R/getClusterTypes.r")
source("R/calcTSNEGeneSpace.r")
source("R/setCellTypeColors.r")
source("R/setClusterColors.r")
source("R/calcUMAPGeneSpace.r")

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

experimentType <- "WT"
source("R/seuratNorm.r")

seuratObj <- ipmc

seuratObj <- SetAllIdent( seuratObj, id = "genCellTypeIdent")

geneSetPlotDir	 <- file.path( geneSpacePlotDir, seuratObj@misc)
dir.create( geneSetPlotDir, showWarnings = FALSE) 

TSNESeedWT <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( geneSetPlotDir, "TSNESeedWT.txt"), TSNESeedWT, "\n")

seuratObj <- calcTSNEGeneSpace( seuratObj, TSNESeedWT)
#seuratObj	<- RunTSNE( seuratObj, genes.use = rownames( seuratObj@data), perplexity = 20, 


UMAPSeedWT <- as.numeric(as.POSIXct(Sys.time()))

#UMAPSeed <- 20

cat( file = file.path( geneSetPlotDir, "UMAPSeed.txt"), UMAPSeedWT, "\n")

#seruatObj 	<- RunUMAP( seuratObj, genes.use = rownames( seuratObj@data), n_neighbors = 15L, min.dist = 0.01, metric = "cosine", seed.use = UMAPSeedWT)

seuratObj <- calcUMAPGeneSpace( seuratObj, UMAPSeedWT, 2)						#the last parameters is UMAP dimension


#Prepare allGenes plots

	tsnePlotWT 	<- DimPlot( object = seuratObj, reduction.use = 'tsne', cols.use = setCellTypeColors( seuratObj), pt.size = 2, do.return = TRUE)
	tsnePlotWT 	<- tsnePlotWT +
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

	umapPlotWT	<- DimPlot(object = seuratObj, reduction.use = 'umap',  cols.use = setCellTypeColors( seuratObj), pt.size = 2, do.return = TRUE)
	umapPlotWT 	<- umapPlotWT +
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

experimentType <- "AllCells"
source("R/seuratNorm.r")

seuratObj <- ipmc
seuratObj <- SetAllIdent( seuratObj, id = "genCellTypeIdent")

geneSetPlotDir	 <- file.path( geneSpacePlotDir, seuratObj@misc)
dir.create( geneSetPlotDir, showWarnings = FALSE) 

TSNESeedAllCells <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( geneSetPlotDir, "TSNESeedAllCells.txt"), TSNESeedAllCells, "\n")

seuratObj <- calcTSNEGeneSpace( seuratObj, TSNESeedAllCells)
#seuratObj	<- RunTSNE( seuratObj, genes.use = rownames( seuratObj@data), perplexity = 20, 


UMAPSeedAllCells <- as.numeric(as.POSIXct(Sys.time()))

#UMAPSeedAllCells <- 20

cat( file = file.path( geneSetPlotDir, "UMAPSeedAllCells.txt"), UMAPSeedAllCells, "\n")

#seruatObj 	<- RunUMAP( seuratObj, genes.use = rownames( seuratObj@data), n_neighbors = 15L, min.dist = 0.01, metric = "cosine", seed.use = UMAPSeedAllCells)

seuratObj <- calcUMAPGeneSpace( seuratObj, UMAPSeedAllCells, 2)						#the last parameters is UMAP dimension

	tsnePlotAllCells 	<- DimPlot( object = seuratObj, reduction.use = 'tsne', cols.use = setCellTypeColors( seuratObj), pt.size = 2, do.return = TRUE)
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

	umapPlotAllCells	<- DimPlot(object = seuratObj, reduction.use = 'umap',  cols.use = setCellTypeColors( seuratObj), pt.size = 2, do.return = TRUE)
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
	seuratObjDenoise	<- seuratObj
	seuratObjDenoise@data	<- apply( seuratObjDenoise@data, c(1,2), function(x) if(x>noiseTol) x else 0)

	hpfClusterTreePlot <- function() {	seuratObj <- SetAllIdent( seuratObj, id = 'originalCellTypes') 
					   	seuratObj <- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = FALSE)
					   	PlotClusterTree( seuratObj, type = "phylogram", cex = 2); nodelabels( text = "  ")}
	genClusterTreePlot <- function() {	seuratObj <- SetAllIdent( seuratObj, id = 'genCellTypeIdent')
					   	seuratObj <- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = FALSE)
						PlotClusterTree( seuratObj, type = "phylogram", cex = 2); nodelabels( text = "  ")}

	dotPlot 	<- DotPlot(seuratObjDenoise, genes.plot = rownames(seuratObjDenoise@data), x.lab.rot = TRUE, dot.scale = 10, 
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
				tsnePlotWT + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				umapPlotWT + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				nrow = 1,
				labels = c("A", "B"),
				label_size = 40)
	secondLine	<- plot_grid( 
				tsnePlotAllCells + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				umapPlotAllCells + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
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
				dotPlot + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				nrow = 1,
				labels = c("E"),
				label_size = 40)
	 
	gridPlot 	<- plot_grid( 
				firstLine  + theme( plot.margin = unit( c( 1, 0, 3, 0), "inches")),  
				secondLine + theme( plot.margin = unit( c( 1, 0, 3, 0), "inches")),  
				thirdLine  + theme( plot.margin = unit( c( 1, 0, 3, 0), "inches")),  
				forthLine  + theme( plot.margin = unit( c( 1, 0, 3, 0), "inches")),  
					labels = '',
					rel_heights = c(3, 3, 2, 3),
					ncol = 1) +
			   theme( plot.margin = unit( c( 1, 1, 5, 1), "inches"))


png( file.path( geneSetPlotDir, "InitCellTypesPlots.png"), width = 2480, height = 3506)
	plot( gridPlot)
dev.off()



#levels(seuratObj@ident) <- c(levels(seuratObj@ident), "G")
#seuratObj@ident[ grep("general", names(seuratObj@ident))] <- "G"
#seuratObj@ident <- droplevels(seuratObj@ident)
#source("R/plotInitCellTypePCAs.r")
#png( file.path( PCAPlotDirName, "geneSpacePlotsDir.png"), width = 480, height = 640)
#	plotInitCellTypePCAs( seuratObj, 5)
#dev.off() 

seuratObj <- SetAllIdent(seuratObj, id = "originalCellTypes")

#plotInitCellTypePCAs(seuratObj, 6)	#Plot PCA diagrams with cell colors, uses its own directorial structure

#remove values, that are too close to zero


