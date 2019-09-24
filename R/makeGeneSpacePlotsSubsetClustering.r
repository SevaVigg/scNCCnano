# Requiers seuratObj created by seuratNorm.r
#
#this snippet is use to make TSNE, UMAP and PCA plots for initial cell types and clusters in the gene space
#it requires to run first seuratNorm.r, which among others sets the appropriate experimentType ( allCells, WT, WT_sox10)
#in seuratObj@project.name
#
# Directory structure :     Res <- Plots <- geneSpace
# 
# before clustering for vizualization tSNE plot must be prepared by makeGeneSpacePlotsSubset
# which also prepares hpfIdent and genCellTypeIdent columns in seuratObj@meta.data


#makeGeneSpacePlotsSubsetClustering <- function( seuratObj){

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

geneSetPlotDir	 <- file.path( geneSpacePlotDir, "subsetGenes")
dir.create( geneSetPlotDir, showWarnings = FALSE)


#pca based clusters

seuratObj 	<- FindClusters( seuratObj, reduction.type = 'pca', dims.use = 1:7, k.param = 10, print.output = FALSE, force.recalc = TRUE, resolution = 0.6)
seuratObj 	<- BuildClusterTree( seuratObj, pcs.use = 1:7, do.plot = FALSE, do.reorder = TRUE) 
#This functions sorts clusters accourding to their size, so we need to assign cluster types

clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj		<- StashIdent( seuratObj, save.name = 'pcaClusters')

#now we are going to cluster in HiD UMAP spaceall plots are related to clustering

seuratObj		<- RunUMAP( seuratObj, genes.use = rownames( seuratObj@data), max.dim = 14, reduction.name = 'umap', n_neighbors = 15L, 
				min_dist = 0.1, metric = "cosine", seed.use = 2 )
seuratObj 		<- FindClusters( seuratObj, reduction.type = 'umap', dims.use = 1:14, k.param = 7, print.output = FALSE, force.recalc = TRUE, resolution = 0.666666)
seuratObj 		<- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = TRUE) 
clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj		<- StashIdent( seuratObj, save.name = 'umap14Clusters')

#we have spoiled umap needed for visualization, so we recalculate umap

seuratObj		<- RunUMAP( seuratObj, genes.use = rownames( seuratObj@data), max.dim = 2, reduction.name = 'umap', n_neighbors = 15L, 
				min_dist = 0.15, metric = "cosine", seed.use = 2 )

#now prepare plots

plotList <- list()

# 	1. Initial cell types with tSNE

seuratObj		<- SetAllIdent( seuratObj, id = 'genCellTypeIdent')

	tsnePlotCells 	<- DimPlot( object = seuratObj, reduction.use = 'tsne', cols.use = setCellTypeColors( seuratObj), pt.size = 2, do.return = TRUE)
	tsnePlotCells 	<- tsnePlotCells +
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

# 	2. Clusters with tSNE

seuratObj		<- SetAllIdent( seuratObj, id = 'umap14Clusters')

	tsnePlotClusters <- DimPlot( object = seuratObj, reduction.use = 'tsne', cols.use = setClusterColors( seuratObj), pt.size = 2, do.return = TRUE)
	tsnePlotClusters <- tsnePlotClusters +
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

#	3. Initial cell types with UMAP

seuratObj	<- SetAllIdent( seuratObj, id = 'genCellTypeIdent')

	umapPlotCells	<- DimPlot(object = seuratObj, reduction.use = 'umap',  cols.use = setCellTypeColors( seuratObj), pt.size = 2, do.return = TRUE)
	umapPlotCells 	<- umapPlotCells +
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
#	4. Clusters with UMAP

seuratObj	<- SetAllIdent( seuratObj, id = 'umap14Clusters')

	umapPlotClusters <- DimPlot(object = seuratObj, reduction.use = 'umap',  cols.use = setClusterColors( seuratObj), pt.size = 2, do.return = TRUE)
	umapPlotClusters <- umapPlotClusters +
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
#	5. Cluster tree

	umapClusterTreePlot <- function() {	seuratObj <- SetAllIdent( seuratObj, id = 'umap14Clusters') 
					   	seuratObj <- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = FALSE)
						par( mar = c(5,5,5,5))
					   	PlotClusterTree( seuratObj, type = "phylogram", cex = 2); nodelabels( text = "  ")}
#	6. DotPlot for clusters

	noiseTol		<- log2(19)
	seuratObjDenoise	<- seuratObj
	seuratObjDenoise@data	<- apply( seuratObjDenoise@data, c(1,2), function(x) if(x>noiseTol) x else 0)


	seuratObjDenoise <- SetAllIdent( seuratObj, id = 'umap14Clusters')
	dotPlotRes 	 <- DotPlot(seuratObjDenoise, genes.plot = rownames(seuratObjDenoise@data), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE)
	dotPlotRes	<- dotPlotRes +
		theme(
			legend.position="none", 
			axis.text.y = element_text( size = 30),
			axis.text.x = element_text( size = 20, angle = 90),
			axis.title  = element_text( size = 25, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)	 

#Prepare the final panel of plots

	firstLine	<- plot_grid( 
				tsnePlotCells + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				umapPlotCells + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				nrow = 1,
				labels = c("A", "B"),
				label_size = 25)
	secondLine	<- plot_grid( 
				tsnePlotClusters + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				umapPlotClusters + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				nrow = 1, 
				labels = c("C", "D"),
				label_size = 25,
				rel_widths = c( 1, 1))
	thirdLine	<- plot_grid(
				umapClusterTreePlot, 
				dotPlotRes + theme( plot.margin = unit( c( 0, 0.5, 0, 0.5), "inches")), 
				nrow = 1,
				labels = c("E", "F"),
				label_size = 25,
				rel_widths = c(1, 3))
 
	gridPlot 	<- plot_grid( 
				firstLine  + theme( plot.margin = unit( c( 0, 0, 1, 1), "inches")),  
				secondLine + theme( plot.margin = unit( c( 0, 0, 1, 1), "inches")),  
				thirdLine  + theme( plot.margin = unit( c( 0, 0, 1, 0), "inches")),  
					labels = '',
					rel_heights = c(1, 1),
					ncol = 1) +
			   theme( plot.margin = unit( c( 1, 1, 5, 1), "inches"))


png( file.path( geneSetPlotDir, "geneSubsetClusteringPlots.png"), width = 1536, height = 2048)
	plot( gridPlot)
dev.off()



#levels(seuratObj@ident) <- c(levels(seuratObj@ident), "G")
#seuratObj@ident[ grep("general", names(seuratObj@ident))] <- "G"
#seuratObj@ident <- droplevels(seuratObj@ident)
#source("R/plotInitCellTypePCAs.r")
#png( file.path( PCAPlotDirName, "geneSpacePlotsDir.png"), width = 480, height = 640)
#	plotInitCellTypePCAs( seuratObj, 5)
#dev.off() 

#plotInitCellTypePCAs(seuratObj, 6)	#Plot PCA diagrams with cell colors, uses its own directorial structure

#remove values, that are too close to zero

return( seuratObj)

#}
