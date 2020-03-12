makeUmapClusters	<- function( seuratObj, umapDim, myResolution){

require("Seurat")
source("R/getClusterTypes.r")


seuratObj 		<- FindClusters( seuratObj, reduction.type = 'umap', 
					dims.use = 1:umapDim, 
					k.param = 15, print.output = FALSE, force.recalc = TRUE, 
					resolution = myResolution) 
seuratObj 		<- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = TRUE) 
clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj		<- StashIdent( seuratObj, save.name = paste0( umapDim, "D_UMAP_res_", myResolution))
return(seuratObj)

}
