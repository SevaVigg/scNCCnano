getClusterTypes <- function( seuratObj){

#This funciton takes a factor whose names are target cell types and whose values are clusters

source("R/calcTargetClusterQuals.r")

clFactor     <- seuratObj@ident

cellClusters <- sapply( levels( clFactor), function(x) WhichCells( seuratObj, ident = x))
nClust	     <- length( levels( clFactor))

cl_IP	<- which.max( calcTargetClusterQuals( seuratObj, "I")) 
cl_MC	<- which.max( calcTargetClusterQuals( seuratObj, "M")) 

clusterTypes 			<- levels( clFactor) 
if (length(unique(c(cl_IP, cl_MC))) < 2) {cat( "Some reference clusters conincide \n")
							 attr(clusterTypes, which = "success") <- FALSE
							 return( clusterTypes)}
names( clusterTypes)		<- levels( clFactor)
names( clusterTypes)[ which( names(cl_IP) == levels(clFactor))]		<- "I"
names( clusterTypes)[ which( names(cl_MC) == levels(clFactor))]		<- "M"
attr(clusterTypes, which = "success") <- TRUE

return(clusterTypes)
}
