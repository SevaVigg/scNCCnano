calcUMAPGeneSpace <- function( seuratObj, UMAPRandSeed){

seuratObj <- RunUMAP( seuratObj, genes.use = rownames(seuratObj@data), seed.use = UMAPRandSeed, 
	n_neighbors = 20, min_dist = 0.2)

return( seuratObj)
}

