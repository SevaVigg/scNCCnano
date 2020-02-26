initWTtsneAll	<- function( seuratWT, seuratAll){

initWT 		<- seuratWT@dr$tsne@cell.embeddings
tsneMin 	<- apply(initWT, 2, min)
tsneMax		<- apply(initWT, 2, max)

newCells 	<- setdiff( colnames(seuratAll@data), colnames(seuratWT@data))

initAll		<- data.frame( 
			initCellsX = runif( ncol( seuratAll@data), min = tsneMin[1], max = tsneMax[1]),
			initCellsY = runif( ncol( seuratAll@data), min = tsneMin[2], max = tsneMax[2]), 
			row.names = colnames( seuratAll@data)	
			)

initAll[ rownames(initWT) , ] <- initWT

return( initAll)
}
