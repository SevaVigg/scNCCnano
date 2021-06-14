findBestUmapClusters	<- function( seuratObj, Spread){

# this snippets enumerates umap parameters to find the best clustering calculated by makeUmapClusters.r
# clustering quality is calculated with the reference to the control cell types M and I in calcTargetClusterQuals.r
# the entire process is rather time consuming - umually more than 20h at my MacBook Pro
# Written by Vsevolod J. Makeev, 2017 - 2021


require("Seurat")
source("R/makeUmapClusters.r")
source("R/getClusterTypes.r")
source("R/calcTargetClusterQuals.r")
source("R/calcUmapGeneSpace.r")

#initialize parameters
minUmapDim 	<- 2
maxUmapDim	<- 8
minMyResolution	<- 0.5
maxMyResolution	<- 6
minMinDist	<- 1
maxMinDist	<- 4
nReplicates	<- 3


umapDims	<- seq( maxUmapDim, minUmapDim)
minDists	<- seq( minMinDist, maxMinDist, by = 0.1) 
resolutions	<- seq( minMyResolution, maxMyResolution, by = 0.2)
#umapRandSeed	<- 42
umapRandSeed	<- as.numeric(as.POSIXct(Sys.time()))

#initialize variables
currSeurObj 	<- StashIdent( seuratObj, save.name = "currentClust")
bestSeurObj	<- StashIdent( seuratObj, save.name = "bestClust")

currSeurObj 	<- calcUmapGeneSpace( currSeurObj, Dim = minUmapDim, myNeighbors = 25L, mySpread = Spread,  
				minDist = minMinDist,  UMAPRandSeed = umapRandSeed, experimentType <- "allCells")$All

currSeurObj	<- makeUmapClusters( currSeurObj, umapDim = minUmapDim, myResolution = minMyResolution)

currParams	<- list( Resolution = minMyResolution, umapDim = minUmapDim, minDist = minMinDist) 
bestParams	<- currParams
bestClusterQual <- calcTargetClusterQuals( currSeurObj, targetCellType = "M")["M"] * calcTargetClusterQuals( currSeurObj, targetCellType = "I")["I"]


for (minDist in minDists) {	
   for (umapDim in umapDims) { 
      for (replicate in 1:nReplicates){ cat( "Replicate ", replicate, "\n")
	currSeurObj <- calcUmapGeneSpace( currSeurObj, Dim = umapDim, myNeighbors = 25L, 
				minDist = minDist,  UMAPRandSeed = as.numeric(as.POSIXct(Sys.time()))
, experimentType <- "allCells", mySpread = Spread)$All

	for( myResolution in resolutions){
		cat( "Finding best UMAP cluster, umapDim = ", umapDim, " myResolution = ", myResolution, minDist, "minDist", "\n")
		currParams$Resolution = myResolution; currParams$umapDim = umapDim; currParams$minDist = minDist
		currSeurObj 	<- makeUmapClusters( currSeurObj, umapDim, myResolution)
		currClusterQual <- calcTargetClusterQuals( currSeurObj, targetCellType = "M")["M"] * calcTargetClusterQuals( currSeurObj, targetCellType = "I")["I"]
		cat( "Cluster Quality = ", currClusterQual, " Best Cluster Quality = ", bestClusterQual,  "\n")
	if( currClusterQual > bestClusterQual) { 
		cat("do switching  ")
		bestClusterQual <- currClusterQual
		bestParams 	<- currParams
		bestSeurObj	<- StashIdent( currSeurObj, save.name = paste0( "Best_D$", bestParams$umapDim, "_R$", bestParams$Resolution, "_minDist$", bestParams$minDist))
		cat("New current best D = ", bestParams$umapDim, "for Resolution = ", bestParams$Resolution, "minDist =", bestParams$minDist, "\n")
		}else{
	 	cat("Weak combination, the current best D = ", bestParams$umapDim, "for Resolution = ", bestParams$Resolution, "minDist =", bestParams$minDist, "\n")} 
		}	#for myResolution
	}	#nReplicates
     }	#umapDim
}	#minDist
cat( "Best cluster quality ", bestClusterQual, "for D = ", bestParams$umapDim, " R = ", bestParams$Resolution, " minDist = ", bestParams$minDist, "\n")

return( bestSeurObj) 
} #makeBestPcaClusters
