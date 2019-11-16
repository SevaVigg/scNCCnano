# Requiers seuratObj created by seuratNorm.r
#
#this snippet is use to make TSNE, UMAP and PCA plots for initial cell types and clusters in the gene space
#it requires to run first seuratNorm.r, which among others sets the appropriate experimentType ( allCells, WT, WT_sox10)
#in seuratObj@project.name
#
# Directory structure :     Res <- Plots <- geneSpace
# 
# before clustering for vizualization tSNE plot must be prepared by makeGeneSpacePlotsSubset
# which also prepares hpfIdent and generalCellTypes columns in seuratObj@meta.data


makeGeneSpacePlotsClustering <- function( seuratObj){

source("R/getClusterTypes.r")
source("R/setCellTypeColors.r")
source("R/setClusterColors.r")

library(ape)

if(!require("gridGraphics")){
install.packages("gridGraphics")
library("gridGraphics")}

require(cowplot)
resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypePlotDir <- file.path(plotDir, seuratObj@project.name)
dir.create( experimentTypePlotDir, showWarnings = FALSE)

geneSpacePlotDir <- file.path( experimentTypePlotDir, "geneSpacePlots")
dir.create( geneSpacePlotDir, showWarnings = FALSE)

seuratObj	<- SetAllIdent( seuratObj, id = "originalCellTypes")

levels(seuratObj@ident) <- c(levels(seuratObj@ident), "G")
seuratObj@ident[ grep("general", names(seuratObj@ident))] <- "G"
seuratObj@ident <- droplevels(seuratObj@ident)

seuratObj 	<- StashIdent( seuratObj, save.name = 'generalCellTypes')

#gene Space based clusters

seuratObj	<- FindClusters( seuratObj, genes.use = rownames(seuratObj@data), k.param = 12, print.output = FALSE, force.recalc = TRUE, resolution = 1.4)
seuratObj 	<- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = TRUE) 
#This functions sorts clusters accourding to their size, so we need to assign cluster types

clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj		<- StashIdent( seuratObj, save.name = 'geneClusters')

#pca based clusters

seuratObj 	<- FindClusters( seuratObj, reduction.type = 'pca', dims.use = 1:8, k.param = 10, print.output = FALSE, force.recalc = TRUE, resolution = 0.6)
seuratObj 	<- BuildClusterTree( seuratObj, pcs.use = 1:8, do.plot = FALSE, do.reorder = TRUE) 
#This functions sorts clusters accourding to their size, so we need to assign cluster types

clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj		<- StashIdent( seuratObj, save.name = 'pcaClusters')

#now we are going to cluster in HiD UMAP spaceall plots are related to clustering

umapDim			<- 5 

seuratObj		<- RunUMAP( seuratObj, genes.use = rownames( seuratObj@data), max.dim =  umapDim, reduction.name = 'umap', n_neighbors = 20L, 
				min_dist = 0.3, metric = "cosine", seed.use = 2, spread = 1)
seuratObj 		<- FindClusters( seuratObj, reduction.type = 'umap', dims.use = 1:umapDim, k.param = 15, print.output = FALSE, force.recalc = TRUE, resolution = .7) 
seuratObj 		<- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = TRUE) 
clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj		<- StashIdent( seuratObj, save.name = 'umapHiDClusters')

#may be some strange tSNE was calculated by other scripts, thus we recalculate tSNE as well

seuratObj		<- RunTSNE( seuratObj, genes.use = rownames( seuratObj@data), perplexity = 20, theta = 0.0, eta = 100, reduction.name = 'tsne', 
				max_iter = 500, seed.use = 142 )


#we have spoiled 2D umap needed for visualization, so we recalculate umap in 2D. Before that we keep the hi-dimensional umap to use it to make lineages

visSeed			<- 142 

seuratObjHiDUmap	<- seuratObj
seuratObj		<- RunUMAP( seuratObj, genes.use = rownames( seuratObj@data), max.dim = 2, reduction.name = 'umap', n_neighbors = 15L, 
				min_dist = 0.4, metric = "cosine", seed.use = visSeed, spread = 1 )

#now prepare plots

plotList <- list()

# 	1. Initial cell types with tSNE

seuratObj		<- SetAllIdent( seuratObj, id = 'generalCellTypes')

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
		) + xlab( label = "tSNE1")+ylab( label = "tSNE2") + ggtitle("Initial cell types")

#	2. Initial cell types with UMAP

seuratObj	<- SetAllIdent( seuratObj, id = 'generalCellTypes')

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
		)  + ggtitle("Initial cell types")

#	3. Visualize Gene Set Clusters with tSNE

seuratObj		<- SetAllIdent( seuratObj, id = 'geneClusters')

	tsnePlotGeneClusters <- DimPlot( object = seuratObj, reduction.use = 'tsne', cols.use = setClusterColors( seuratObj), pt.size = 2, do.return = TRUE)
	tsnePlotGeneClusters <- tsnePlotGeneClusters +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			axis.title.x = element_text( margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 0,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		) + xlab( label = "tSNE1")+ylab( label = "tSNE2") + ggtitle("Clustering by gene expressions")


#	4. Visualize Gene Set Clusters with UMAP

seuratObj	<- SetAllIdent( seuratObj, id = 'geneClusters')

	umapPlotGeneClusters <- DimPlot(object = seuratObj, reduction.use = 'umap',  cols.use = setClusterColors( seuratObj), pt.size = 2, do.return = TRUE)
	umapPlotGeneClusters <- umapPlotGeneClusters +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			axis.title.x = element_text( margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 0,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		)  + ggtitle("Clustering by gene expressions")

		  
#	5. Visualize PCA Clusters with tSNE

seuratObj		<- SetAllIdent( seuratObj, id = 'pcaClusters')

	tsnePlotPCAClusters <- DimPlot( object = seuratObj, reduction.use = 'tsne', cols.use = setClusterColors( seuratObj), pt.size = 2, do.return = TRUE)
	tsnePlotPCAClusters <- tsnePlotPCAClusters +
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
	xlab( label = "tSNE1")+ylab( label = "tSNE2") + ggtitle("Clustering by PCA")


#	6. Visualize PCA Clusters with UMAP

seuratObj	<- SetAllIdent( seuratObj, id = 'pcaClusters')

	umapPlotPCAClusters <- DimPlot(object = seuratObj, reduction.use = 'umap',  cols.use = setClusterColors( seuratObj), pt.size = 2, do.return = TRUE)
	umapPlotPCAClusters <- umapPlotPCAClusters +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			axis.title.x = element_text( margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 0,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		) + ggtitle("Clustering by PCA")


#	7. Visualize UMAP Clusters with tSNE

seuratObj		<- SetAllIdent( seuratObj, id = 'umapHiDClusters')

	tsnePlotUMAPClusters <- DimPlot( object = seuratObj, reduction.use = 'tsne', cols.use = setClusterColors( seuratObj), pt.size = 2, do.return = TRUE)
	tsnePlotUMAPClusters <- tsnePlotUMAPClusters +
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
	xlab( label = "tSNE1")+ylab( label = "tSNE2") + ggtitle( paste0("Clustering by UMAP_", umapDim))


#	8. Visualize UMAP Clusters with UMAP

seuratObj	<- SetAllIdent( seuratObj, id = 'umapHiDClusters')

	umapPlotUMAPClusters <- DimPlot(object = seuratObj, reduction.use = 'umap',  cols.use = setClusterColors( seuratObj), pt.size = 2, do.return = TRUE)
	umapPlotUMAPClusters <- umapPlotUMAPClusters +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			axis.title.x = element_text( margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 0,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		) + ggtitle( paste0("Clusteryg by UMAP_", umapDim))


#	9. Genes Cluster tree

	geneClusterTreePlot <- function() {	seuratObj <- SetAllIdent( seuratObj, id = 'geneClusters') 
					   	seuratObj <- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = FALSE)
					   	PlotClusterTree( seuratObj, type = "phylogram", cex = 2, main = "Raw gene expressions")
						par( mar = c(5, 5, 5, 5)); nodelabels( text = "  ")}
#	9. PCA Cluster tree

	pcaClusterTreePlot <- function() {	seuratObj <- SetAllIdent( seuratObj, id = 'pcaClusters') 
					   	seuratObj <- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = FALSE)
						PlotClusterTree( seuratObj, type = "phylogram", cex = 2)
						par( mar = c(5,5,5,5)); nodelabels( text = "  "); title( "PCA")}

#	10. UMAP Cluster tree

	umapClusterTreePlot <- function() {	seuratObj <- SetAllIdent( seuratObj, id = 'umapHiDClusters') 
					   	seuratObj <- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = FALSE)
					   	PlotClusterTree( seuratObj, type = "phylogram", cex = 2)
						par( mar = c(5,5,5,5)); nodelabels( text = "  "); title( paste0("UMAP_", umapDim))}
#	11. DotPlot for clusters

	noiseTol		<- log2(19)
	seuratObjDenoise	<- seuratObj
	seuratObjDenoise@data	<- apply( seuratObjDenoise@data, c(1,2), function(x) if(x>noiseTol) x else 0)


	seuratObjDenoise <- SetAllIdent( seuratObj, id = 'umapHiDClusters')
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
				tsnePlotCells + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				umapPlotCells + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				nrow = 1,
				align = "h",
				labels = c("A", "B"),
				label_size = 25,
				rel_widths = c( 1, 1)
				)
	secondLine	<- plot_grid( 
				tsnePlotGeneClusters + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				umapPlotGeneClusters + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				nrow = 1, 
				align = "h",
				labels = c("C", "D"),
				label_size = 25,
				rel_widths = c( 1, 1)
				)
	thirdLine	<- plot_grid( 
				tsnePlotPCAClusters + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				umapPlotPCAClusters + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				nrow = 1,
				align = "h", 
				labels = c("E", "F"),
				label_size = 25,
				rel_widths = c( 1, 1)
				)
	forthLine	<- plot_grid( 
				tsnePlotUMAPClusters + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				umapPlotUMAPClusters + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				nrow = 1, 
				align = "h",
				labels = c("G", "H"),
				label_size = 25,
				rel_widths = c( 1, 1)
				)
	fifthLine	<- plot_grid(
				geneClusterTreePlot,
				pcaClusterTreePlot,
				umapClusterTreePlot, 
				nrow = 1,
				align = "h",
				labels = c("E", "F", "G"),
				label_size = 25,
				rel_widths = c(1, 1, 1)
				)
	sixthLine	<- plot_grid(
				dotPlotRes + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				nrow = 1,
				align = "h",
				labels = c("G"),
				label_size = 25
#				rel_widths = (1, 1, 1))
 				)
	gridPlot 	<- plot_grid( 
				firstLine  + theme( plot.margin = unit( c( 0.3, 0, 0.3, 1), "inches")),  
				secondLine + theme( plot.margin = unit( c( 0.3, 0, 0.3, 1), "inches")),  
				thirdLine + theme( plot.margin = unit( c( 0.3, 0, 0.3, 1), "inches")),  
				forthLine + theme( plot.margin = unit( c( 0.3, 0, 0., 1), "inches")),  
				fifthLine + theme( plot.margin = unit( c( 0, 0, 0.3, 1), "inches")),  
				sixthLine  + theme( plot.margin = unit( c( 0.3, 0, 0.3, 1), "inches")),  
					align = "v",
					labels = '',
					rel_heights = c(1, 1),
					ncol = 1) +
			   theme( plot.margin = unit( c( 1, 1, 1, 1), "inches"))
			

png( file.path( geneSpacePlotDir, paste0("geneSpaceClusteringPlots_umap", umapDim, ".png")), width = 1536, height = 2048)
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

seuratObj@misc <- append( seuratObj@misc, visSeed)

return( seuratObj)

}
