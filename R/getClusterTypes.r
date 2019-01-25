getClusterTypes <- function( seuratObj){

#This funciton takes a factor whose names are cells and whose values are clusters and assigns colors to clusters with predominant cell types
if(!require(data.table)){
  install.packages("data.table")}

clFactor     <- seuratObj@ident

cellClusters <- sapply( unique(levels(clFactor)), function(x) unlist( transpose( strsplit( names( clFactor[clFactor == x]), "\\." ))[[1]]))
nClust	     <- length(levels( clFactor))


cl_IP	<- which.max( sapply(cellClusters, function(x) if(length(x[grep("I", x)])==0) 0 else ( 1 + length( grep("I", x)))/sqrt(1+length(x))))
cl_MC	<- which.max( sapply(cellClusters, function(x) if(length(x[grep("M", x)])==0) 0 else ( 1 + length( grep("M", x)))/sqrt(1+length(x))))
cl_tail	<- which.max( sapply(cellClusters, function(x) if(length(x[grep("Tl", x)])==0) 0 else ( 1 + length( grep("Tl", x)))/sqrt(1+length(x))))

clusterTypes 			<- levels( clFactor) 
if (length(unique(c(cl_IP, cl_MC, cl_tail))) < 3) {cat( "Some reference clusters conincide \n")
							 attr(clusterTypes, which = "success") <- TRUE
							 return( clusterTypes)}
names( clusterTypes)		<- levels( clFactor)
names( clusterTypes)[ which( names(cl_IP) == levels(clFactor))]		<- "I"
names( clusterTypes)[ which( names(cl_MC) == levels(clFactor))]		<- "M"
names( clusterTypes)[ which( names(cl_tail) == levels(clFactor))]	<- "Tl"

attr(clusterTypes, which = "success") <- TRUE

return(clusterTypes)
}
