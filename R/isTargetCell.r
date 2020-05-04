isTargetCell <- function( seuratObj, cellName, targetCellType){
	XgeneList <- c("pax7a", "pax7b")
	LgeneList <- c( "ltk", "pax7b", "tfec","phox2b", "tyr")
	if( targetCellType %in% c("M", "I", "Tl")) {return( grepl( targetCellType, cellName))}else{ 	
	if( targetCellType == "X") { sum( seuratObj@data[ XgeneList, cellName]) > 8} else{ 
	if( targetCellType == "L") { sum( seuratObj@data[ LgeneList, cellName]) > 4*length(LgeneList)} else return( FALSE)}}
}
