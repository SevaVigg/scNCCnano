findBestUmapClusters	<- function( seuratObj){

require("Seurat")
source("R/makeUmapClusters.r")
source("R/getClusterTypes.r")
source("R/calcTargetClusterQuals.r")
source("R/calcUmapGeneSpace.r")

#initialize parameters
Spread		<- 8
minUmapDim 	<- 3
maxUmapDim	<- 7
minMyResolution	<- 3
maxMyResolution	<- 5
minMinDist	<- 3
maxMinDist	<- Spread

umapDims	<- seq( minUmapDim, maxUmapDim)
minDists	<- seq( minMinDist, maxMinDist, by = 0.1) 
resolutions	<- seq( minMyResolution, maxMyResolution, by = 0.2)
umapRandSeed	<- 2

#initialize variables
currSeurObj 	<- StashIdent( seuratObj, save.name = "currentClust")
bestSeurObj	<- StashIdent( seuratObj, save.name = "bestClust")

currSeurObj 	<- calcUmapGeneSpace( currSeurObj, Dim = minUmapDim, myNeighbors = 20L, mySpread = Spread,  
				minDist = minMinDist,  UMAPRandSeed = umapRandSeed, experimentType <- "allCells")$All

currSeurObj	<- makeUmapClusters( currSeurObj, umapDim = minUmapDim, myResolution = minMyResolution)

currParams	<- list( Resolution = minMyResolution, UMAPdim = minUmapDim, minDist = minMinDist) 
bestParams	<- currParams
bestClusterQual <- calcTargetClusterQuals( currSeurObj, targetCellType = "M")["M"] * calcTargetClusterQuals( currSeurObj, targetCellType = "I")["I"]

for (minDist in minDists) {	
   for (umapDim in umapDims) { 
	currSeurObj <- calcUmapGeneSpace( currSeurObj, Dim = umapDim, myNeighbors = 20L, 
				minDist = minDist,  UMAPRandSeed = umapRandSeed, experimentType <- "allCells", mySpread = Spread)$All

	for( myResolution in resolutions){
		cat( "Finding best UMAP cluster, umapDim = ", umapDim, " myResolution = ", myResolution, minDist, "minDist", "\n")
		currParams$Resolution = myResolution; currParams$UMAPDim = umapDim; currParams$minDist = minDist
		currSeurObj 	<- makeUmapClusters( currSeurObj, umapDim, myResolution)
		currClusterQual <- calcTargetClusterQuals( currSeurObj, targetCellType = "M")["M"] * calcTargetClusterQuals( currSeurObj, targetCellType = "I")["I"]
		cat( "Cluster Quality = ", currClusterQual, "\n")
	if( currClusterQual > bestClusterQual) { 
		bestClusterQual <- currClusterQual
		bestParams 	<- currParams
		bestSeurObj	<- StashIdent( currSeurObj, save.name = paste0( "Best_D", bestParams$UMAPDim, "_R", bestParams$Resolution, "_minDist", bestParams$minDist))
		cat("Current best D = ", bestParams$PCADim, "Resolution", bestParams$Resolution, "minDist ", bestParams$minDist, "\n")}else{ cat("Weak combination\n")} 
		}	#for myResolution
	}	#umapDim
     }	#minDist
cat( "Best cluster quality ", bestClusterQual, "for D = ", bestParams$UMAPDim, " R = ", bestParams$Resolution, " minDist = ", minDist, "\n")

return( bestSeurObj) 
} #makeBestPcaClusters
