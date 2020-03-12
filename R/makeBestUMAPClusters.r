makeBestUMAPClusters	<- function( seuratObj){

require("Seurat")
source("R/makeUmapClusters.r")
source("R/getClusterTypes.r")
source("R/calcTargetClusterQuals.r")

minUmapDim 	<- 2
maxUmapDim	<- 12
minMyResolution	<- 0.2
maxMyResolution	<- 1.2

umapDims	<- seq( minUmapDim, maxUmapDim)
resolutions	<- seq( minMyResolution, maxMyResolution, by = 0.05)

#initialize variables
currSeurObj 	<- StashIdent( seuratObj, save.name = "currentClust")
bestSeurObj	<- StashIdent( seuratObj, save.name = "bestClust")

currSeurObj 	<- calcUMAPGeneSpace( currSeurObj, Dim = minUmapDim, myNeighbors = 15L, 
				minDist = 0.3,  UMAPRandSeed = 42, experimentType <- "allCells")$All

currSeurObj	<- makeUmapClusters( currSeurObj, umapDim = minUmapDim, myResolution = minMyResolution)

currParams	<- list( Resolution = minMyResolution, UMAPdim = minUmapDim) 
bestParams	<- currParams
bestClusterQual <- max( calcTargetClusterQuals( currSeurObj, targetCellType = "M"))*max( calcTargetClusterQuals( currSeurObj, targetCellType = "I"))

for (umapDim in umapDims) { 
	currSeurObj <- calcUMAPGeneSpace( currSeurObj, Dim = umapDim, myNeighbors = 15L, 
				minDist = 0.3,  UMAPRandSeed = 42, experimentType <- "allCells")$All

	for( myResolution in resolutions){
		cat( "Finding best UMAP cluster, umapDim = ", umapDim, " myResolution = ", myResolution, "\n")
		currParams$Resolution = myResolution; currParams$UMAPDim = umapDim
		currSeurObj 	<- makeUmapClusters( currSeurObj, umapDim, myResolution)
		currClusterQual <- max( calcTargetClusterQuals( currSeurObj, targetCellType = "M"))*max( calcTargetClusterQuals( currSeurObj, targetCellType = "I"))
		cat( "Cluster Quality = ", currClusterQual, "\n")
	if( currClusterQual > bestClusterQual) { 
		bestClusterQual <- currClusterQual
		bestParams 	<- currParams
		bestSeurObj	<- StashIdent( currSeurObj, save.name = paste0( "Best_D", bestParams$UMAPDim, "_R", bestParams$Resolution))
		cat("Current best D = ", bestParams$PCADim, "Resolution", bestParams$Resolution, "\n")}else{ cat("Weak combination\n")} 
		}	#for myResolution
	}	#umapDim
return( bestSeurObj) 
} #makeBestPcaClusters
