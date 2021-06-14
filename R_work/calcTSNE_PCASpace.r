#Draw initial TSNEPlot

calcTSNE_PCASpace <- function( seuratObj, comps, TSNErandSeed){

seuratObj <- RunTSNE( seuratObj, reduction.use = "pca", dims.use = 1:comps, seed.use = TSNErandSeed, 
	theta = 0, eta = 10, max_iter = 3000, perplexity = 20, verbose = FALSE)

return( seuratObj)
}




