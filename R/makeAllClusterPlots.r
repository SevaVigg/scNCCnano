# Requiers seuratAll created by seuratNorm.r
#
#this snippet is use to make TSNE, UMAP and PCA plots for initial cell types and clusters in the gene space
#it requires to run first seuratNorm.r, which among others sets the appropriate experimentType ( allCells, WT, WT_sox10)
#in seuratAll@project.name
#
# Directory structure :     Res <- Plots <- Clustering 
# 
# before clustering for vizualization tSNE plot must be prepared by makeGeneSpacePlotsSubset
# which also prepares hpfIdent and genCellTypeIdent columns in seuratAll@meta.data


source("R/seuratNorm.r")
source("R/getClusterTypes.r")
source("R/setCellTypeColors.r")
source("R/setClusterColors.r")
source("R/calcUmapGeneSpace.r")
source("R/calcTSNEGeneSpace.r")
source("R/makeUMAPclusters.r")
source("R/makePCAGeneClusters.r")

library(ape)

if(!require("gridGraphics")){
install.packages("gridGraphics")
library("gridGraphics")}

require(cowplot)
resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

clusteringPlotDir 	<- file.path( plotDir, "allClusteringPlots")
dir.create( clusteringPlotDir, showWarnings = FALSE)

clusterResName		<- file.path( resDir, "frozenClusters")
dir.create( clusterResName, showWarning = FALSE) 

seuratWT		<- seuratNorm("WT")
seuratAll		<- seuratNorm("allCells")

#make trees after dimension reduction with UMAP

umapDim			<- 5 
pcaDim			<- 6 

umapRes			<- calcUMAPGeneSpace( seuratAll, seurWT = seuratWT,  Dim = umapDim, myNeighbors = 15L, 
				minDist = 0.3,  UMAPRandSeed = 42, experimentType <- "allCondWT")

myResolutionPCAwt 	<- 0.9
myResolutionPCAall	<- 0.8
myResolutionUMAP	<- 1

seuratWT		<- makeUmapClusters( umapRes$WT,  umapDim, myResolutionUMAP)
seuratAll		<- makeUmapClusters( umapRes$All, umapDim, myResolutionUMAP)
seuratWT		<- makePCAGeneClusters( seuratWT, pcaDim, myResolutionPCAwt)
seuratAll		<- makePCAGeneClusters( seuratAll, pcaDim, myResolutionPCAall)
#seuratWT		<- makeBestPCAClusters( seuratWT)
#seuratAll		<- makeBestPCAClusters( seuratAll)



#may be some strange tSNE was calculated by other scripts, thus we recalculate tSNE as well

seuratWT		<- calcTSNEGeneSpace( seuratWT, TSNErandSeed = 42) 
seuratAll		<- calcTSNEGeneSpace( seuratAll, TSNErandSeed = 42, initSeurObj = seuratWT)

#we have spoiled 2D umap needed for visualization, so we recalculate umap in 2D. Before that we save the hi-dimensional umap to use it to make lineages

visSeed			<- 42 

seuratWTHiDUmap		<- seuratWT
umapRes			<- calcUMAPGeneSpace( seuratAll, seurWT = seuratWT, Dim = 2, myNeighbors = 20L, 
				minDist = 0.6, UMAPRandSeed = visSeed, experimentType <- "allCondWT")
seuratWT		<- umapRes$WT
seuratAll		<- umapRes$All

#now we need the results for all cells conditioned by WT


plotList <- list()

# 	1. Initial cell types with tSNE, WT only

seuratWT		<- SetAllIdent( seuratWT, id = 'genCellTypeIdent')

	tsnePlotCellsWT 	<- DimPlot( object = seuratWT, reduction.use = 'tsne', cols.use = setCellTypeColors( seuratWT), pt.size = 2, do.return = TRUE)
	tsnePlotCellsWT 	<- tsnePlotCellsWT +
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

#	2. Inital cell types with tSNE, All cells	

seuratAll		<- SetAllIdent( seuratAll, id = 'genCellTypeIdent')

	tsnePlotCellsAll 	<- DimPlot( object = seuratAll, reduction.use = 'tsne', cols.use = setCellTypeColors( seuratAll), pt.size = 2, do.return = TRUE)
	tsnePlotCellsAll 	<- tsnePlotCellsAll +
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

#	3. Initial cell types with UMAP, WT only 

seuratWT	<- SetAllIdent( seuratWT, id = 'genCellTypeIdent')

	umapPlotCellsWT	<- DimPlot(object = seuratWT, reduction.use = 'umap',  cols.use = setCellTypeColors( seuratWT), pt.size = 2, do.return = TRUE)
	umapPlotCellsWT 	<- umapPlotCellsWT +
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

#	4. Initial cell types with UMAP, All cells

seuratAll	<- SetAllIdent( seuratAll, id = 'genCellTypeIdent')

	umapPlotCellsAll	<- DimPlot(object = seuratAll, reduction.use = 'umap',  cols.use = setCellTypeColors( seuratAll), pt.size = 2, do.return = TRUE)
	umapPlotCellsAll 	<- umapPlotCellsAll +
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

#	5. UMAP Clusters with tSNE, WT

seuratWT		<- SetAllIdent( seuratWT, id = paste0( umapDim, "D_UMAP_res_", myResolutionUMAP))

	tsnePlotUMAPClustersWT <- DimPlot( object = seuratWT, reduction.use = 'tsne', cols.use = setClusterColors( seuratWT), pt.size = 2, do.return = TRUE)
	tsnePlotUMAPClustersWT <- tsnePlotUMAPClustersWT +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			axis.title.x = element_text( margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 0,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		) + xlab( label = "tSNE1")+ylab( label = "tSNE2") + ggtitle("Clustering by UMAP")

#	6. Clusters with tSNE, All

seuratAll		<- SetAllIdent( seuratAll, id = paste0( umapDim, "D_UMAP_res_", myResolutionUMAP))

	tsnePlotUMAPClustersAll <- DimPlot( object = seuratAll, reduction.use = 'tsne', cols.use = setClusterColors( seuratAll), pt.size = 2, do.return = TRUE)
	tsnePlotUMAPClustersAll <- tsnePlotUMAPClustersAll +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			axis.title.x = element_text( margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 0,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		) + xlab( label = "tSNE1")+ylab( label = "tSNE2") + ggtitle("Clustering by UMAP")


#	7. Clusters with UMAP, WT

seuratWT		<- SetAllIdent( seuratWT, id = paste0( umapDim, "D_UMAP_res_", myResolutionUMAP))

	umapPlotUMAPClustersWT <- DimPlot( object = seuratWT, reduction.use = 'umap', cols.use = setClusterColors( seuratWT), pt.size = 2, do.return = TRUE)
	umapPlotUMAPClustersWT <- umapPlotUMAPClustersWT +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			axis.title.x = element_text( margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 0,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		) + xlab( label = "UMAP1")+ylab( label = "UMAP2") + ggtitle("Clustering by UMAP")

#	8. Clusters with UMAP, All

seuratAll		<- SetAllIdent( seuratAll, id = paste0( umapDim, "D_UMAP_res_", myResolutionUMAP))

	umapPlotUMAPClustersAll <- DimPlot( object = seuratAll, reduction.use = 'umap', cols.use = setClusterColors( seuratAll), pt.size = 2, do.return = TRUE)
	umapPlotUMAPClustersAll <- umapPlotUMAPClustersAll +
#		xlim( -0.22, 0.12) +
#		ylim( -0.22, 0.12) + 
		theme(
			legend.position="none", 
			axis.text = element_text( size = 30),
			axis.title  = element_text( size = 25, face = "bold"),
			axis.title.x = element_text( margin = margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")),
			axis.title.y = element_text( margin = margin(t = 0,  r = 20, b = 0, l = 0, unit = "pt")),
			panel.background = element_rect(fill = "gray60")
		) + xlab( label = "UMAP1")+ylab( label = "UMAP2") + ggtitle("Clustering by UMAP")


#	9. Visualize PCA Clusters with tSNE, WT

seuratWT		<- SetAllIdent( seuratWT, id = paste0( pcaDim, "D_PCA_res_", myResolutionPCAwt))

	tsnePlotPCAClustersWT <- DimPlot( object = seuratWT, reduction.use = 'tsne', cols.use = setClusterColors( seuratWT), pt.size = 2, do.return = TRUE)
	tsnePlotPCAClustersWT <- tsnePlotPCAClustersWT +
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

		  
#	10. Visualize PCA Clusters with tSNE, All

seuratAll		<- SetAllIdent( seuratAll, id = paste0( pcaDim, "D_PCA_res_", myResolutionPCAall))

	tsnePlotPCAClustersAll <- DimPlot( object = seuratAll, reduction.use = 'tsne', cols.use = setClusterColors( seuratAll), pt.size = 2, do.return = TRUE)
	tsnePlotPCAClustersAll <- tsnePlotPCAClustersAll +
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


#	11. Visualize PCA Clusters with UMAP, WT

seuratWT	<- SetAllIdent( seuratWT, id = paste0( pcaDim, "D_PCA_res_", myResolutionPCAwt))

	umapPlotPCAClustersWT <- DimPlot(object = seuratWT, reduction.use = 'umap',  cols.use = setClusterColors( seuratWT), pt.size = 2, do.return = TRUE)
	umapPlotPCAClustersWT <- umapPlotPCAClustersWT +
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
	xlab( label = "UMAP1")+ylab( label = "UMAP2") + ggtitle("Clustering by PCA")



#	12. Visualize PCA Clusters with UMAP, All

seuratAll	<- SetAllIdent( seuratAll, id = paste0( pcaDim, "D_PCA_res_", myResolutionPCAall))

	umapPlotPCAClustersAll <- DimPlot(object = seuratAll, reduction.use = 'umap',  cols.use = setClusterColors( seuratAll), pt.size = 2, do.return = TRUE)
	umapPlotPCAClustersAll <- umapPlotPCAClustersAll +
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
	xlab( label = "UMAP1")+ylab( label = "UMAP2") + ggtitle("Clustering by PCA")

#	9. Genes Cluster tree

	geneClusterTreePlot <- function() {	seuratAll <- SetAllIdent( seuratAll, id = 'geneClusters') 
					   	seuratAll <- BuildClusterTree( seuratAll, genes.use = rownames(seuratAll@data), do.plot = FALSE, do.reorder = FALSE)
					   	PlotClusterTree( seuratAll, type = "phylogram", cex = 2, main = "Raw gene expressions")
						par( mar = c(5, 5, 5, 5)); nodelabels( text = "  ")}
#	9. PCA Cluster tree

	pcaClusterTreePlot <- function() {	seuratAll <- SetAllIdent( seuratAll, id = paste0( pcaDim, "D_PCA_res_", myResolutionPCAall)) 
					   	seuratAll <- BuildClusterTree( seuratAll, genes.use = rownames(seuratAll@data), do.plot = FALSE, do.reorder = FALSE)
						PlotClusterTree( seuratAll, type = "phylogram", cex = 2)
						par( mar = c(5,5,5,5)); nodelabels( text = "  "); title( "PCA")}

#	10. UMAP Cluster tree

	umapClusterTreePlot <- function() {	seuratAll <- SetAllIdent( seuratAll, id = 'umapHiDClusters') 
					   	seuratAll <- BuildClusterTree( seuratAll, genes.use = rownames(seuratAll@data), do.plot = FALSE, do.reorder = FALSE)
					   	PlotClusterTree( seuratAll, type = "phylogram", cex = 2)
						par( mar = c(5,5,5,5)); nodelabels( text = "  "); title( paste0("UMAP_", umapDim))}
#	11. DotPlot for clusters

	noiseTol		<- log2(19)
	seuratWTDenoise		<- seuratWT
	seuratWTDenoise@data	<- apply( seuratWTDenoise@data, c(1,2), function(x) if(x>noiseTol) x else 0)


	seuratWTDenoise  <- SetAllIdent( seuratWTDenoise, id = paste0( umapDim, "D_UMAP_res_", myResolutionUMAP))
	dotPlotRes 	 <- DotPlot(seuratWTDenoise, genes.plot = rownames(seuratWTDenoise@data), x.lab.rot = TRUE, dot.scale = 10, 
					plot.legend = TRUE, dot.min = 0, scale.by = "radius", do.return = TRUE)
	dotPlotRes	<- dotPlotRes +
		theme(
			legend.position="none", 
			axis.text.y = element_text( size = 20),
			axis.text.x = element_text( size = 20, angle = 90),
			axis.title  = element_text( size = 25, face = "bold"),
#panel.background = element_rect(fill = "gray90")
			)	 

#Prepare the final panel of plots

	tsneCellsLine	<- plot_grid( 
				tsnePlotCellsWT  + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				tsnePlotCellsAll + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				nrow = 1,
				align = "h",
				labels = c("A", "B"),
				label_size = 25,
				rel_widths = c( 1, 1)
				)
	tsneUMAPLine	<- plot_grid( 
				tsnePlotUMAPClustersWT + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				tsnePlotUMAPClustersAll + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				nrow = 1,
				align = "h", 
				labels = c("E", "F"),
				label_size = 25,
				rel_widths = c( 1, 1)
				)
	tsnePCALine	<- plot_grid( 
				tsnePlotPCAClustersWT + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				tsnePlotPCAClustersAll + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				nrow = 1,
				align = "h", 
				labels = c("E", "F"),
				label_size = 25,
				rel_widths = c( 1, 1)
				)
	umapCellsLine	<- plot_grid( 
				umapPlotCellsWT  + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				umapPlotCellsAll + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				nrow = 1, 
				align = "h",
				labels = c("C", "D"),
				label_size = 25,
				rel_widths = c( 1, 1)
				)
	umapUMAPLine	<- plot_grid( 
				umapPlotUMAPClustersWT  + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				umapPlotUMAPClustersAll + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				nrow = 1, 
				align = "h",
				labels = c("G", "H"),
				label_size = 25,
				rel_widths = c( 1, 1)
				)
	umapPCALine	<- plot_grid( 
				umapPlotPCAClustersWT  + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				umapPlotPCAClustersAll + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				nrow = 1, 
				align = "h",
				labels = c("G", "H"),
				label_size = 25,
				rel_widths = c( 1, 1)
				)
	clTreesLine	<- plot_grid(
				geneClusterTreePlot,
				pcaClusterTreePlot,
				umapClusterTreePlot, 
				nrow = 1,
				align = "h",
				labels = c("E", "F", "G"),
				label_size = 25,
				rel_widths = c(1, 1, 1)
				)
	dotPlotLine	<- plot_grid(
				dotPlotRes + theme( plot.margin = unit( c( 0, 0.3, 0, 0.3), "inches")), 
				nrow = 1,
				align = "h",
				labels = c("H"),
				label_size = 25
#				rel_widths = (1, 1, 1))
 				)
	gridPlot 	<- plot_grid( 
				tsneCellsLine  + theme( plot.margin = unit( c( 0.3, 0, 0.3, 1), "inches")),  
				tsneUMAPLine + theme( plot.margin = unit( c( 0.3, 0, 0.3, 1), "inches")),  
				tsnePCALine + theme( plot.margin = unit( c( 0.3, 0, 0.3, 1), "inches")),  
				umapCellsLine + theme( plot.margin = unit( c( 0.3, 0, 0., 1), "inches")),  
				umapUMAPLine + theme( plot.margin = unit( c( 0, 0, 0.3, 1), "inches")),  
				umapPCALine  + theme( plot.margin = unit( c( 0.3, 0, 0.3, 1), "inches")),  
					align = "v",
					labels = '',
					rel_heights = c(1, 1),
					ncol = 1) +
			   theme( plot.margin = unit( c( 1, 1, 1, 1), "inches"))
			

png( file.path( clusteringPlotDir, paste0("geneSpaceClusteringPlots_umap", umapDim, ".png")), width = 1536, height = 2048)
	plot( gridPlot)
dev.off()



#levels(seuratAll@ident) <- c(levels(seuratAll@ident), "G")
#seuratAll@ident[ grep("general", names(seuratAll@ident))] <- "G"
#seuratAll@ident <- droplevels(seuratAll@ident)
#source("R/plotInitCellTypePCAs.r")
#png( file.path( PCAPlotDirName, "geneSpacePlotsDir.png"), width = 480, height = 640)
#	plotInitCellTypePCAs( seuratAll, 5)
#dev.off() 

#plotInitCellTypePCAs(seuratAll, 6)	#Plot PCA diagrams with cell colors, uses its own directorial structure

#remove values, that are too close to zero

seuratAll@misc <- append( seuratAll@misc, visSeed)

if( askYesNo("Do you want to save seurat object to clusterFreeze folder ?")) save( seuratWT, seuratAll, file = file.path( clusterResName, "seurAllWTfreeze"))
#return( seuratAll)


