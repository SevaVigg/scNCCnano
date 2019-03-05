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

pcaPlotDir 	<- file.path( experimentTypeDir, "PCADimReduction")
dir.create( pcaPlotDir, showWarnings = FALSE)

compsDir 	<- file.path( pcaPlotDir, paste0("comps", 6))
dir.create(compsDir, showWarnings = FALSE)

resolDir	<- file.path( compsDir, paste0( "resol_", ipmc@calc.params[tail(grep("FindClusters", names(ipmc@calc.params), value = TRUE),1)][[1]]$resolution))
dir.create( resolDir, showWarnings = FALSE)

coordsMD	<- as.matrix(ipmc@dr$pca@cell.embeddings)
#these are initial coordinates

source("R/createSlingShotObject.r")
ipmcSling 	<- createSlingShotObject( coordsMD, ipmc)

source("R/plot2DallLineages.r")

png( file.path( resolDir, "clustersLineages.png"))
plot2DallLineages( ipmcSling, ipmc)
dev.off() 

