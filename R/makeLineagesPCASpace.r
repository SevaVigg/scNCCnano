#This snippet requires precomputed tSNE for ipmc object created with seuratNorm.r
#It requires precomputed tSNE in PCA space with makePCADimRedPlots.r


if(!require(slingshot)){
  install.packages("slingshot")
}

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypeDir <- file.path(plotDir, ipmc@project.name)
dir.create( experimentTypeDir, showWarnings = FALSE)

pcaPlotDir 	<- file.path( experimentTypeDir, "PCAdimReduction")
dir.create( pcaPlotDir, showWarnings = FALSE)

compsDir 	<- file.path( pcaPlotDir, paste0("comps", 6))
dir.create(compsDir, showWarnings = FALSE)

#lineagePCAPlotDir 	<- file.path( compsDir, "lineagesPCASpace")
#dir.create( lineagePCAPlotDir, showWarnings = FALSE)

tSNEValsMD	<- as.matrix(ipmc@dr$pca@cell.embeddings)
#these are initial coordinates

slingObjMD 	<- slingshot(tSNEValsMD, ipmc@ident, start.clus = clTypes["Tl"],end.clus=c(clTypes["I"], clTypes["M"]))
MC_linId	<- which( slingObjMD@slingParams$end.clus == clTypes["M"])
MC_linName	<- paste0("Lineage", MC_linId)
IP_linId	<- which( slingObjMD@slingParams$end.clus == clTypes["I"])
IP_linName	<- paste0("Lineage", IP_linId)

plotVals	<- ipmc@dr$tsne@cell.embeddings

source("R/getLineageCoords.r")
LineageTree	<- getLineageCoords(ipmc, slingObjMD)  

linPlotDir <- ( file.path( lineagePCAPlotDir, "lineagePlots"))
dir.create( linPlotDir,  showWarnings = FALSE)

source("R/plot2DallLineages.r")

png( file.path( compsDir, "lineagePlot.png"))
	plot2DallLineages( LineageTree, plotVals, clTypes)
dev.off()

#the following part plot lineage plots for all lineages

#source("R/lineageVlnPlot.r")
#source("R/plot2DidLineage.r")
#source("R/plot2Dcells.r")

#allGenes <- "sox10" # NB !!!

#	for (lineageId in seq(1, length( slingObjMD@lineages))){
#		linIdPlotDir <- file.path( linPlotDir, names( slingObjMD@lineages)[lineageId] )
#		dir.create( linIdPlotDir, showWarnings = FALSE)
#		png( file.path( linIdPlotDir, paste0( "Lineage", lineageId, "_plot.png")))
 #       		plot2Dcells( plotVals, ipmc@ident, clTypes, linIdPlotDir)
#			plot2DidLineage( LineageTree, lineageId)
#		dev.off()	
#Now plot violin plots for allGenes for the lineages
#		for (gene in allGenes) lineageVlnPlot( ipmcMD, slingObjMD, gene, lineageId, MC_linId, IP_linId, linIdPlotDir)
#	}

#} #comps
#} #resolDec 


cat("Now getting principle curves\n")

psTime		<- slingPseudotime(slingObjMD)
source("R/getPseudoOrder.r")
curves		<- getPseudoOrder(slingObjMD)



