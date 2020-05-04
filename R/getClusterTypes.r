getClusterTypes <- function( seuratObj){

#This funciton takes a factor whose names are target cell types and whose values are clusters

source("R/calcTargetClusterQuals.r")

clFactor     <- seuratObj@ident

cellClusters <- sapply( levels( clFactor), function(x) WhichCells( seuratObj, ident = x))
nClust	     <- length( levels( clFactor))

cl_IP	<- which.max( calcTargetClusterQuals( seuratObj, "I")) 
cl_MC	<- which.max( calcTargetClusterQuals( seuratObj, "M")) 
cl_tail	<- which.max( calcTargetClusterQuals( seuratObj, "Tl")) 
cl_X	<- which.max( calcTargetClusterQuals( seuratObj, "X"))
cl_L	<- which.max( calcTargetClusterQuals( seuratObj, "L"))

clusterTypes 			<- levels( clFactor) 
if (length(unique(c(cl_IP, cl_MC, cl_tail, cl_X, cl_L))) < 5) {cat( "Some reference clusters conincide \n")
							 attr(clusterTypes, which = "success") <- FALSE
							 return( clusterTypes)}
names( clusterTypes)		<- levels( clFactor)
names( clusterTypes)[ which( names(cl_IP) == levels(clFactor))]		<- "I"
names( clusterTypes)[ which( names(cl_MC) == levels(clFactor))]		<- "M"
names( clusterTypes)[ which( names(cl_tail) == levels(clFactor))]	<- "E"
names( clusterTypes)[ which( names(cl_X) == levels(clFactor))]		<- "X"
names( clusterTypes)[ which( names(cl_L) == levels(clFactor))]		<- "L"

attr(clusterTypes, which = "success") <- TRUE

return(clusterTypes)
}
