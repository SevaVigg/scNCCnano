plotInitCellTypePCAs <- function( seuratObj, Ncomps){

source("R/setClusterColors.r")

clIdent <- seuratObj@ident

if (Ncomps = 6) 

plotList <- list()
for (c1 in 1:(Ncomps-1)){
	for (c2 in (c1+1):Ncomps){
	plotS	<- PCAPlot( seuratObj, dim.1 = c1, dim.2 = c2, cols.use = setClusterColors( seuratObj))	
	plotList <- c( plotList, plotS)
	}
}
   
}


