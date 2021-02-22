isTargetCell <- function( seuratObj, cellName, targetCellType){
	XgeneList <- c("pax7a", "pax7b")
	ltHMPgeneList <- c( "ltk", "foxo1a")
	if( targetCellType %in% c("M", "I", "Tl")) {return( grepl( targetCellType, cellName))}else{ 	
	if( targetCellType == "X") { sum( seuratObj@data[ XgeneList, cellName]) > 4*length( XgeneList)} else{ 
	if( targetCellType == "ltHMP") { sum( seuratObj@data[ ltHMPgeneList, cellName]) > 3*length( ltHMPgeneList)} else return( FALSE)}}
}
