isTargetCell <- function( seuratObj, cellName, targetCellType){
	XgeneList <- c("pax7a", "pax7b")
	LgeneList <- c( "ltk", "foxp4", "mitfa" )
	if( targetCellType %in% c("M", "I", "Tl")) {return( grepl( targetCellType, cellName))}else{ 	
	if( targetCellType == "X") { sum( seuratObj@data[ XgeneList, cellName]) > 4*length(XgeneList)} else{ 
	if( targetCellType == "L") { sum( seuratObj@data[ LgeneList, cellName]) > 3*length(LgeneList)} else return( FALSE)}}
}
