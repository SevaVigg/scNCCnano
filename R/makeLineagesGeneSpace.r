#This snippet requires precomputed tSNE for ipmc object created with seuratNorm.r
#It requires precomputed dimension reduction like tsne or umap, which is selected by DimRed variable
# 

makeLineagesGeneSpace <- function(seuratObj, DimRed){

if(!require(slingshot)){
  install.packages("slingshot")
}

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypeDir <- file.path(plotDir, ipmc@project.name)
dir.create( experimentTypeDir, showWarnings = FALSE)

geneSpacePlotsDir <- file.path( experimentTypeDir, "geneSpacePlots")
dir.create( geneSpacePlotsDir, showWarnings = FALSE)

#resolDir	<- file.path( geneSpacePlotDir, paste0( "resol_", ipmc@calc.params[tail(grep("FindClusters", names(ipmc@calc.params), value = TRUE),1)][[1]]$resolution))
#resolDir	<- geneSpacePlotDir
#dir.create( resolDir, showWarnings = FALSE)

geneSubsetPlotDir <- file.path( geneSpacePlotsDir, seuratObj@misc)
dir.create( geneSubsetPlotDir, showWarnings = FALSE)

source("R/createSlingShotObject.r")
slingShotObject	<- createSlingShotObject( seuratObj, DimRed)

#source("R/plot2DallLineages.r")

#png( file.path( resolDir, "clustersLineagesTSNE.png"))
#plot2DallLineages( ipmcSling, ipmc, "tsne")
#dev.off() 

#png( file.path( resolDir, "clustersLineagesUMAP.png"))
#plot2DallLineages( ipmcSling, ipmc, "umap")
#dev.off() 

return( slingShotObject)
}

