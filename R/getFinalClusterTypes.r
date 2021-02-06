getFinalClusterTypes <- function( seuratObj){

#This funciton takes a factor whose names are target cell types and whose values are clusters

source("R/calcTargetClusterQuals.r")

clFactor     	<- seuratObj@ident
clusterTypes	<- levels( clFactor)
nClust	     	<- length( clusterTypes)
	

cellClusters <- sapply( clusterTypes, function(x) WhichCells( seuratObj, ident = x))

cl_IP	<- which.max( calcTargetClusterQuals( seuratObj, "I")) 
cl_MC	<- which.max( calcTargetClusterQuals( seuratObj, "M")) 
cl_tail	<- which.max( calcTargetClusterQuals( seuratObj, "Tl")) 
cl_X	<- which.max( calcTargetClusterQuals( seuratObj, "X"))
cl_HMP	<- which.max( calcTargetClusterQuals( seuratObj, "HMP"))

clusterTypes 			<- levels( clFactor) 
if (length(unique(c(cl_IP, cl_MC))) < 2) {cat( "Some reference clusters conincide \n")
							 attr(clusterTypes, which = "success") <- FALSE
							 return( clusterTypes)}
#if (length(unique(c(cl_IP, cl_MC, cl_tail, cl_X, cl_HMP))) < 5) {cat( "Some reference clusters conincide \n")
#							 attr(clusterTypes, which = "success") <- FALSE
#							 return( clusterTypes)}

#We must clear old assigned names 

names( clusterTypes)							<- clusterTypes
names( clusterTypes)[ which( names(clusterTypes) == "I")]		<- "old_I"
names( clusterTypes)[ which( names(clusterTypes) == "M")]		<- "old_M"
names( clusterTypes)[ which( names(clusterTypes) == "eNCC")]		<- "old_eNCC"
names( clusterTypes)[ which( names(clusterTypes) == "X")]		<- "old_X"
names( clusterTypes)[ which( names(clusterTypes) == "HMP")]		<- "old_HMP"



names( clusterTypes)[ which( names(cl_IP) == clusterTypes)]		<- "I"
names( clusterTypes)[ which( names(cl_MC) == clusterTypes)]		<- "M"
names( clusterTypes)[ which( names(cl_tail) == clusterTypes)]		<- "eNCC"
names( clusterTypes)[ which( names(cl_X) == clusterTypes)]		<- "X"
names( clusterTypes)[ which( names(cl_HMP) == clusterTypes)]		<- "HMP"

attr(clusterTypes, which = "success") <- TRUE

return(clusterTypes)
}
