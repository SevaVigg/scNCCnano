calcTargetClusterQuals <- function( seuratObj, targetCellType = "M"){

source("R/isTargetCell.r")

cellClusters 	<- sapply( levels( seuratObj@ident), function(x) WhichCells( seuratObj, ident = x)) 
clusterQuals	<- sapply( cellClusters, function(x) { nCells <- length(which( sapply( x, isTargetCell, seuratObj = seuratObj, targetCellType = targetCellType))); if(nCells ==0) 0 else ( 1 + nCells)/sqrt(1+ length(x))})

return( clusterQuals)
}
