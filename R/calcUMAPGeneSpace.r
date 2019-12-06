calcUMAPGeneSpace <- function( seuratObj, UMAPRandSeed, Dim){

seuratObj <- RunUMAP( seuratObj, genes.use = rownames(seuratObj@data), max.dim = Dim, seed.use = UMAPRandSeed, 
	n_neighbors = 20L, min_dist = 0.005, metric = "cosine")

return( seuratObj)
}

