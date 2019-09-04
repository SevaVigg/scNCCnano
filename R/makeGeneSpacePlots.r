# Requiers ipmc created by seuratNorm.r
#
#this snippet is use to make TSNE, UMAP and PCA plots for initial cell types and clusters in the gene space
#it requires to run first seuratNorm.r, which among others sets the appropriate experimentType ( allCells, WT, WT_sox10)
#in ipmc@project.name
#
# Directory structure :     Res <- Plots <- geneSpace

source("R/getClusterTypes.r")
source("R/calcTSNEGeneSpace.r")
source("R/setCellTypeColors.r")
source("R/setClusterColors.r")
source("R/calcUMAPGeneSpace.r")

if(!require("gridGraphics")){
install.packages("gridGraphics")
library("gridGraphics")}


resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypePlotDir <- file.path(plotDir, ipmc@project.name)
dir.create( experimentTypePlotDir, showWarnings = FALSE)

geneSpacePlotsDir <- file.path( experimentTypePlotDir, "geneSpacePlots")
dir.create( geneSpacePlotsDir, showWarnings = FALSE)


levels(ipmc@ident) <- c(levels(ipmc@ident), "G")
ipmc@ident[ grep("general", names(ipmc@ident))] <- "G"
ipmc@ident <- droplevels(ipmc@ident)



ipmc <- BuildClusterTree( ipmc, genes.use = rownames(ipmc@data), do.plot = FALSE, do.reorder = FALSE)
#png( file.path( geneSpacePlotsDir, "InitCellTypeTree.png"))
	clusterTreePlot <- function() {PlotClusterTree( ipmc)}
#dev.off()

TSNESeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( geneSpacePlotsDir, "TSNESeed.txt"), TSNESeed, "\n")
#ipmc <- calcTSNEGeneSpace( ipmc, TSNESeed)

UMAPSeed <- as.numeric(as.POSIXct(Sys.time()))
	cat( file = file.path( geneSpacePlotsDir, "UMAPSeed.txt"), TSNESeed, "\n")
#ipmc <- calcUMAPGeneSpace( ipmc, UMAPSeed, 2)						#the last parameters is UMAP dimension


plotList <- list()


	tsnePlot 	<- DimPlot( object = ipmc, reduction.use = 'tsne', cols.use = setCellTypeColors( ipmc), pt.size = 2, do.return = TRUE)
	tsnePlot 	<- tsnePlot +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			panel.background = element_rect(fill = "gray60")
		)

	umapPlot	<- DimPlot(object = ipmc, reduction.use = 'umap',  cols.use = setCellTypeColors( ipmc), pt.size = 2, do.return = TRUE)
	umapPlot 	<- umapPlot +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			panel.background = element_rect(fill = "gray60")
		)


	noiseTol		<- log2(19)
	ipmcDenoise		<- ipmc
	ipmcDenoise@data	<- apply( ipmcDenoise@data, c(1,2), function(x) if(x>noiseTol) x else 0)

	dotPlot 	<- DotPlot(ipmcDenoise, genes.plot = rownames(ipmcDenoise@data), x.lab.rot = TRUE, dot.scale = 5, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE)
	dotPlot		<- dotPlot +
		theme(
			legend.position="none", 
			axis.text.y = element_text( size = 30),
			axis.text.x = element_text( size = 20, angle = 90),
			axis.title  = element_text( size = 25, face = "bold")
			#panel.background = element_rect(fill = "gray90")
		)

	gridPlot 	<- plot_grid( tsnePlot, umapPlot, dotPlot, clusterTreePlot)

png( file.path( geneSpacePlotsDir, "InitCellTypesPlots.png"), width = 1536, height = 2048)
	plot(gridPlot)
dev.off()



#levels(ipmc@ident) <- c(levels(ipmc@ident), "G")
#ipmc@ident[ grep("general", names(ipmc@ident))] <- "G"
#ipmc@ident <- droplevels(ipmc@ident)
#source("R/plotInitCellTypePCAs.r")
#png( file.path( PCAPlotDirName, "geneSpacePlotsDir.png"), width = 480, height = 640)
#	plotInitCellTypePCAs( ipmc, 5)
#dev.off() 

ipmc <- SetAllIdent(ipmc, id = "originalCellTypes")

#plotInitCellTypePCAs(ipmc, 6)	#Plot PCA diagrams with cell colors, uses its own directorial structure

#remove values, that are too close to zero




