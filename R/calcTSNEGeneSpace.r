#calc TSNE for Seurat object
# This function takes initial data from dr@tsne in initSeurObj and initiate by those the coordinates in newSeurObj;
# the cells that are not found in initSeurObj are put at random positions


calcTSNEGeneSpace <- function( newSeurObj, TSNErandSeed = 25, Norm = FALSE, initSeurObj = NULL){

require( proxy)

if(!require("Rtsne")){
install.packages("Rtsne")}
library("Rtsne")

if(!require("beepr")){
install.packages("beepr")}
library("beepr")




if( is.null(initSeurObj)){ initNew <- NULL}else{


initCells 	<-  GetDimReduction(object = initSeurObj, reduction.type = "tsne", 
            		slot = "cell.embeddings")
tsneMin 	<- apply(initCells, 2, min)
tsneMax		<- apply(initCells, 2, max)

newCells 	<- setdiff( colnames(newSeurObj@data), colnames(initSeurObj@data))

#first we put all cells in newSeurObj at random positions
initNew		<- data.frame( 
			initCellsX = runif( ncol( newSeurObj@data), min = tsneMin[1], max = tsneMax[1]),
			initCellsY = runif( ncol( newSeurObj@data), min = tsneMin[2], max = tsneMax[2]), 
			row.names = colnames( newSeurObj@data)	
			)
#and then update those present in initSeurObj with their positions
initNew[ rownames(initCells) , ] <- initCells
}



#now call tsne
distanceMatrix 	<- dist( t( newSeurObj@data), method = "cosine", pairwise = TRUE)
rtsneRes	<- Rtsne( distanceMatrix, is_distance = TRUE, 
				perplexity = 30, 
				theta = 0, 
				eta = 500, 
				pca = FALSE, 
				max_iter = 50000,
				pca_center = FALSE,
				normalize = Norm,
				Y_init = initNew
			)
#and prepare Seurat Obj
tsneData		<- rtsneRes$Y
rownames( tsneData) 	<- colnames( newSeurObj@data)
colnames( tsneData)	<- c("tSNE_1", "tSNE_2")
newSeurObj		<- SetDimReduction( newSeurObj, reduction.type = "tsne", slot = "cell.embeddings", new.data = tsneData)
newSeurObj		 <- SetDimReduction(object = newSeurObj, reduction.type = "tsne", slot = "key", new.data = "tSNE")



#newSeurObj <- RunTSNE( newSeurObj, genes.use = rownames(newSeurObj@data), seed.use = TSNErandSeed, 
#	theta = 0, eta = 10, max_iter = 3000, perplexity = 15, verbose = FALSE)

beep(3)

return( newSeurObj)
}




