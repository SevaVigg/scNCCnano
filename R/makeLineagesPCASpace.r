#This snippet requires precomputed tSNE for seuratObj object created with seuratNorm.r
#It requires precomputed tSNE in PCA space with makePCADimRedPlots.r

makeLineagesPCASpace <- function( seuratObj, comps, clResolution){

if(!require(slingshot)){
  install.packages("slingshot")
}

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypeDir <- file.path(plotDir, seuratObj@project.name)
dir.create( experimentTypeDir, showWarnings = FALSE)

pcaPlotDir 	<- file.path( experimentTypeDir, "PCADimReduction")
dir.create( pcaPlotDir, showWarnings = FALSE)

compsDir 	<- file.path( pcaPlotDir, paste0("comps", comps))
dir.create(compsDir, showWarnings = FALSE)

#resolDir	<- file.path( compsDir, paste0( "resol_", seuratObj@calc.params[tail(grep("FindClusters", names(seuratObj@calc.params), value = TRUE),1)][[1]]$resolution))
resolDir	<- file.path( compsDir, paste0( "resol_", clResolution))
dir.create( resolDir, showWarnings = FALSE)

#coordsMD	<- as.matrix(seuratObj@dr$pca@cell.embeddings)
#these are initial coordinates

#coordsMD	<- coordsMD + jitter(coordsMD)

source("R/createSlingShotObject.r")
seuratObjSling 	<- createSlingShotObject( seuratObj, "pca")

source("R/plot2DallLineages.r")

png( file.path( resolDir, paste0("TSNE_ClusterLineages", comps, "_res", clResolution, ".png")))
plot2DallLineages( seuratObjSling, seuratObj, "tsne")
dev.off() 

png( file.path( resolDir, paste0("UMAP_ClusterLineages", comps, "_res", clResolution, ".png")))
plot2DallLineages( seuratObjSling, seuratObj, "umap")
dev.off() 

heatmapPlotDir 	<- file.path(resolDir, "HeatMaps")
dir.create(heatmapPlotDir, showWarnings = FALSE)

source("R/getCurveHeatMap.r")

png( file.path(heatmapPlotDir, "IPnew.png"), width = 800, height = 600)
IP_id <- which(lapply( seuratObjSling@lineages, tail, 1) == "I")
getCurveHeatMap( seuratObj, seuratObjSling, IP_id, 2)
dev.off()

png( file.path(heatmapPlotDir, "MCnew.png"), width = 800, height = 600)
MC_id <- which(lapply( seuratObjSling@lineages, tail, 1) == "M")
getCurveHeatMap( seuratObj, seuratObjSling, MC_id, 2)
dev.off()


}
