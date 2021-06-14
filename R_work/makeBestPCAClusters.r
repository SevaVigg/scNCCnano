makeBestPCAClusters	<- function( seuratObj){

require("Seurat")
source("R/makePCAGeneClusters.r")
source("R/getClusterTypes.r")
source("R/calcTargetClusterQuals.r")

minPcaDim 	<- 2
maxPcaDim	<- 7
minMyResolution	<- 0.6
maxMyResolution	<- 1

pcaDims		<- seq( minPcaDim, maxPcaDim)
resolutions	<- seq( minMyResolution, maxMyResolution, by = 0.1)

#initialize variables
currSeurObj 	<- StashIdent( seuratObj, save.name = "currentClust")
bestSeurObj	<- StashIdent( seuratObj, save.name = "bestClust")


currSeurObj	<- makePCAGeneClusters( seuratObj, pcaDim = minPcaDim, myResolution = minMyResolution)

currParams	<- list( Resolution = minMyResolution, PCAdim = minPcaDim) 
bestParams	<- currParams
bestClusterQual <- max( calcTargetClusterQuals( currSeurObj, targetCellType = "M"))*max( calcTargetClusterQuals( currSeurObj, targetCellType = "I"))

for (pcaDim in pcaDims) { 
	for( myResolution in resolutions){
	cat( "Finding best PCA cluster, pcaDim = ", pcaDim, " myResolution = ", myResolution, "\n")
	currParams$Resolution = myResolution; currParams$PCADim = pcaDim
	currSeurObj 	<- makePCAGeneClusters( currSeurObj, pcaDim, myResolution)
	currClusterQual <- max( calcTargetClusterQuals( currSeurObj, targetCellType = "M"))*max( calcTargetClusterQuals( currSeurObj, targetCellType = "I"))
	cat( "Cluster Quality = ", currClusterQual, "\n")
	if( currClusterQual > bestClusterQual) { 
		bestClusterQual	<- currClusterQual
		bestParams 	<- currParams
		bestSeurObj	<- StashIdent( currSeurObj, save.name = paste0( "Best_D", bestParams$PCADim, "_R", bestParams$Resolution))
		cat("Current best D = ", bestParams$PCADim, "Resolution", bestParams$Resolution, "\n")}else{ cat("Weak combination\n")} 
		} 		# for myResolution
	}	# pdaDims

return( bestSeurObj) 
} #makeBestPcaClusters
