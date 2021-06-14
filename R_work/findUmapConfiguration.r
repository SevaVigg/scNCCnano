findUmapConfiguration <- function( seuratObj, umapDim = 6, minDist = 4.5, myResolution = 1.4 ){

source("R/calcUmapGeneSpace.r")
source("R/makeUmapClusters.r")
source("R/getFinalClusterTypes.r")

exitCond	<- TRUE

while(exitCond){
Spread		<- 9.5
umapRandSeed 	<- as.numeric(Sys.time())
currSeurObj <- calcUmapGeneSpace( seuratObj, Dim = umapDim, myNeighbors = 25L, minDist = minDist,  UMAPRandSeed = umapRandSeed, experimentType <- "allCells", mySpread = Spread)$All
currSeurObj 	<- makeUmapClusters( currSeurObj, umapDim, myResolution)
levels( currSeurObj@ident) <- names(getFinalClusterTypes( currSeurObj))

slingWT <- createSlingShotObject( currSeurObj, "umap")

logicalResI <- any(unlist(lapply( slingWT@lineages, function(x) x[which(x == "HMP")+1] == "I")))
logicalResM <- any(unlist(lapply( slingWT@lineages, function(x) x[which(x == "HMP")+1] == "M")))
cat( "ResI =", logicalResI, "ResM = ", logicalResM, "\n")

logicalResI <- if( is.na( logicalResI)) FALSE else logicalResI
logicalResM <- if( is.na( logicalResI)) FALSE else logicalResM


exitCond <- !(logicalResI&logicalResM) #break while if both are true

cat("exitCond = ", exitCond)

}
umapConfigRes <- list( logI = logicalResI, logM = logicalResM, lineages = slingWT@lineages, clWT = currSeurObj)

return(umapConfigRes)
}

