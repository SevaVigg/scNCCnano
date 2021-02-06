isTargetCell <- function( seuratObj, cellName, targetCellType){
	XgeneList <- c("pax7a", "pax7b")
	HMPgeneList <- c( "ltk", "foxo1a")
	if( targetCellType %in% c("M", "I", "Tl")) {return( grepl( targetCellType, cellName))}else{ 	
	if( targetCellType == "X") { sum( seuratObj@data[ XgeneList, cellName]) > 4*length(XgeneList)} else{ 
	if( targetCellType == "HMP") { sum( seuratObj@data[ HMPgeneList, cellName]) > 3*length(HMPgeneList)} else return( FALSE)}}
}
