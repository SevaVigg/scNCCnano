calcTargetClusterQuals <- function( seuratObj, targetCellType = "M"){

cellClusters 	<- sapply( levels( seuratObj@ident), function(x) WhichCells( seuratObj, ident = x)) 
#clusterQuals	<- sapply( cellClusters, function(x) if(length(x[grep( targetCellType, x)])==0) 0 else ( 1 + length( grep( targetCellType, x)))/sqrt(1+length(x)))
clusterQuals	<- sapply( cellClusters, function(x) if(length(x[grep( targetCellType, x)])==0) 0 else ( 1 + length( grep( targetCellType, x)))/(1+length(x)))
return( clusterQuals)
}
