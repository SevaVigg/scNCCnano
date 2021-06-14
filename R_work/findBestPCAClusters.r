findBestUmapClusters	<- function( seuratObj){

require("Seurat")
source("R/makeUmapClusters.r")
source("R/getClusterTypes.r")
source("R/calcTargetClusterQuals.r")
source("R/calcUmapGeneSpace.r")

#initialize parameters
minPCADim 	<- 3
maxPCADim	<- 10
minMyResolution	<- 3
maxMyResolution	<- 5

PCADims	<- seq( minPCADim, maxPCADim)

#initialize variables
currSeurObj 	<- StashIdent( seuratObj, save.name = "currentClust")
bestSeurObj	<- StashIdent( seuratObj, save.name = "bestClust")

currSeurObj	<- makeUmapClusters( currSeurObj, PCADim = minPCADim, myResolution = minMyResolution)

currParams	<- list( Resolution = minMyResolution, PCAdim = minPCADim) 
bestParams	<- currParams
bestClusterQual <- calcTargetClusterQuals( currSeurObj, targetCellType = "M")["M"] * calcTargetClusterQuals( currSeurObj, targetCellType = "I")["I"]

   for (PCADim in PCADims) { 
	currSeurObj <- calcUmapGeneSpace( currSeurObj, Dim = PCADim, myNeighbors = 20L, 
				minDist = minDist,  UMAPRandSeed = umapRandSeed, experimentType <- "allCells", mySpread = Spread)$All

	for( myResolution in resolutions){
		cat( "Finding best UMAP cluster, PCADim = ", PCADim, " myResolution = ", myResolution, minDist, "minDist", "\n")
		currParams$Resolution = myResolution; currParams$PCADim = PCADim; currParams$minDist = minDist
		currSeurObj 	<- makeUmapClusters( currSeurObj, PCADim, myResolution)
		currClusterQual <- calcTargetClusterQuals( currSeurObj, targetCellType = "M")["M"] * calcTargetClusterQuals( currSeurObj, targetCellType = "I")["I"]
		cat( "Cluster Quality = ", currClusterQual, "\n")
	if( currClusterQual > bestClusterQual) { 
		bestClusterQual <- currClusterQual
		bestParams 	<- currParams
		bestSeurObj	<- StashIdent( currSeurObj, save.name = paste0( "Best_D", bestParams$PCADim, "_R", bestParams$Resolution, "_minDist", bestParams$minDist))
		cat("Current best D = ", bestParams$PCADim, "Resolution", bestParams$Resolution, "minDist ", bestParams$minDist, "\n")}else{ cat("Weak combination\n")} 
		}	#for myResolution
	}	#PCADim

cat( "Best cluster quality ", bestClusterQual, "for D = ", bestParams$PCADim, " R = ", bestParams$Resolution, " minDist = ", minDist, "\n")

return( bestSeurObj) 
} #makeBestPcaClusters
