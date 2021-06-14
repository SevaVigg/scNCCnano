isTargetCell <- function( seuratObj, cellName, targetCellType){

#This snippet identifies the cell type by expression of target genes or presence of control cell types

geneMax	<- apply( seuratObj@data, 1, max)

if( seuratObj@project.name == "taqman"){
	XgeneList 	<- c("xdh", "pnp4a")
	MgeneList	<- c("tyrp1b", "pnp4a", "mbpa", "sox10")
	IgeneList	<- c("ltk", "pnp4a")
	eHMPgeneList	<- c("sox9b", "snai1b")
	ltHMPgeneList	<- c("ltk", "mitfa", "snai1b")
	
	if( targetCellType == "I") { sum( seuratObj@data[ IgeneList, cellName] / geneMax[IgeneList]) > 0.3 * length( IgeneList)} else{ 
	if( targetCellType == "M") { sum( seuratObj@data[ MgeneList, cellName] / geneMax[MgeneList]) > 0.3 * length( MgeneList)} else{ 
	if( targetCellType == "X") { sum( seuratObj@data[ XgeneList, cellName] / geneMax[XgeneList]) > 0.3 * length( XgeneList)} else{ 
	if( targetCellType == "ltHMP") { sum( seuratObj@data[ ltHMPgeneList, cellName] / geneMax[ltHMPgeneList]) > 0.3 * length( ltHMPgeneList)} else{ 
	if( targetCellType == "eHMP")  { sum( seuratObj@data[  eHMPgeneList, cellName] / geneMax[eHMPgeneList])  > 0.3 * length( eHMPgeneList)} else return( FALSE)}}}}
}else if (seuratObj@project.name == "sox10_mutants"){
	ltHMPgeneList 	<- c( "ltk", "foxo1a")
	eHMPgeneList	<- c( "sox9b", "tfap2e")
	if( targetCellType %in% c("M", "I", "Tl")) {return( grepl( targetCellType, cellName))}else{ 
	if( targetCellType == "ltHMP") { sum( seuratObj@data[ ltHMPgeneList, cellName]) > 3*length( ltHMPgeneList)} else{
	if( targetCellType == "eHMP")  { sum( seuratObj@data[  eHMPgeneList, cellName]) > 3*length(  eHMPgeneList)} else return( FALSE)}}
}else{
	XgeneList 	<- c("pax7a", "pax7b")
	ltHMPgeneList 	<- c( "ltk", "foxo1a")
	if( targetCellType %in% c("M", "I", "Tl")) {return( grepl( targetCellType, cellName))}else{ 	
	if( targetCellType %in% c("eHMP")) {return( grepl( "Tl", cellName))}else{ 	
	if( targetCellType == "X") { sum( seuratObj@data[ XgeneList, cellName]) > 4*length( XgeneList)} else{ 
	if( targetCellType == "ltHMP") { sum( seuratObj@data[ ltHMPgeneList, cellName]) > 3*length( ltHMPgeneList)} else return( FALSE)}}}}
}
