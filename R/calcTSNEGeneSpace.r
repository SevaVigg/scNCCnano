#calc TSNE for Seurat object

calcTSNEGeneSpace <- function( seuratObj, TSNErandSeed, Norm = FALSE, initTSNE = NULL){

require( proxy)
require( Rtsne)
require( beepr)

distanceMatrix 	<- dist( t( seuratObj@data), method = "cosine", pairwise = TRUE)
rtsneRes	<- Rtsne( distanceMatrix, is_distance = TRUE, 
				perplexity = 15, 
				theta = 0, 
				eta = 5000, 
				pca = FALSE, 
				max_iter = 50000,
				pca_center = FALSE,
				normalize = Norm,
				Y_init = initTSNE
			)

tsneData		<- rtsneRes$Y
rownames( tsneData) 	<- colnames( seuratObj@data)
colnames( tsneData)	<- c("tSNE_1", "tSNE_2")
seuratObj		<- SetDimReduction( seuratObj, reduction.type = "tsne", slot = "cell.embeddings", new.data = tsneData)
seuratObj		 <- SetDimReduction(object = seuratObj, reduction.type = "tsne", slot = "key", new.data = "tSNE")



#seuratObj <- RunTSNE( seuratObj, genes.use = rownames(seuratObj@data), seed.use = TSNErandSeed, 
#	theta = 0, eta = 10, max_iter = 3000, perplexity = 15, verbose = FALSE)

beep(3)

return( seuratObj)
}




