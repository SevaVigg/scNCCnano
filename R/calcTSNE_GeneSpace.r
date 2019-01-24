#Draw initial TSNEPlot

calcTSNE_GeneSpace <- function( seuratObj, TSNErandSeed){

seuratObj <- RunTSNE( seuratObj, genes.use = rownames(seuratObj@data), seed.use = TSNErandSeed, 
	theta = 0, eta = 10, max_iter = 3000, perplexity = 20, verbose = FALSE)

return( seuratObj)
}




