getClusterTypes <- function( seuratObj){

# This funciton takes a seurat objects, takes its @ident field and identifies clusters with M and I control cell types as cell types containing 
# the maximal number of control cells
# written by Vsevolod J. Makeev 2017 - 2021

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
