#Draw initial TSNEPlot

calcUMAP_PCASpace <- function( seuratObj, comps, randSeed, dimRed){

seuratObj <- RunUMAP( seuratObj, reduction.use = "pca", dims.use = 1:comps, seed.use = randSeed, 
	n_neighbors = 40, min_dist = 0.01, max.dim = dimRed)

return( seuratObj)
}
