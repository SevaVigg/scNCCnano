#Draw initial TSNEPlot

calcUMAP_PCASpace <- function( seuratObj, comps, randSeed){

seuratObj <- RunUMAP( seuratObj, reduction.use = "pca", dims.use = 1:comps, seed.use = randSeed, 
	n_neighbors = 30, min_dist = 0.6)

return( seuratObj)
}
