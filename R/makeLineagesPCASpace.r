#This snippet requires precomputed tSNE for ipmc object created with seuratNorm.r
#It requires precomputed tSNE in PCA space with makePCADimRedPlots.r
makeLineagesPCASpace <- function( ipmc, comps, clResolution){

if(!require(slingshot)){
  install.packages("slingshot")
}

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypeDir <- file.path(plotDir, ipmc@project.name)
dir.create( experimentTypeDir, showWarnings = FALSE)

pcaPlotDir 	<- file.path( experimentTypeDir, "PCADimReduction")
dir.create( pcaPlotDir, showWarnings = FALSE)

compsDir 	<- file.path( pcaPlotDir, paste0("comps", comps))
dir.create(compsDir, showWarnings = FALSE)

#resolDir	<- file.path( compsDir, paste0( "resol_", ipmc@calc.params[tail(grep("FindClusters", names(ipmc@calc.params), value = TRUE),1)][[1]]$resolution))
resolDir	<- file.path( compsDir, paste0( "resol_", clResolution))
dir.create( resolDir, showWarnings = FALSE)

#coordsMD	<- as.matrix(ipmc@dr$pca@cell.embeddings)
#these are initial coordinates

#coordsMD	<- coordsMD + jitter(coordsMD)

source("R/createSlingShotObject.r")
ipmcSling 	<- createSlingShotObject( ipmc, "pca")

source("R/plot2DallLineages.r")

png( file.path( resolDir, paste0("TSNE_ClusterLineages", comps, "_res", clResolution, ".png")))
plot2DallLineages( ipmcSling, ipmc, "tsne")
dev.off() 

png( file.path( resolDir, paste0("UMAP_ClusterLineages", comps, "_res", clResolution, ".png")))
plot2DallLineages( ipmcSling, ipmc, "umap")
dev.off() 


png( file.path( resolDir, "clustersLineagesUMAP.png"))
plot2DallLineages( ipmcSling, ipmc, "umap")
dev.off() 

}
