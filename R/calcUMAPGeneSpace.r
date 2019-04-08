calcUMAPGeneSpace <- function( seuratObj, UMAPRandSeed, Dim){

seuratObj <- RunUMAP( seuratObj, genes.use = rownames(seuratObj@data), max.dim = Dim, seed.use = UMAPRandSeed, 
	n_neighbors = 40, min_dist = 0.3)

return( seuratObj)
}

