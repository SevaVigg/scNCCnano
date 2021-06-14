makeUmapClusters	<- function( seuratObj, umapDim, myResolution){

# this snipper runs Seurat tool FindClusters in the UMAP image space and identify control cluster types by
# control cells or marker genes 

require("Seurat")
source("R/getClusterTypes.r")

#pcaDim			<- 6
seuratObj 		<- FindClusters( seuratObj, reduction.type = 'umap', 
					dims.use = 1:umapDim, 
					k.param = 15, print.output = FALSE, force.recalc = TRUE, save.SNN = TRUE,
					resolution = myResolution) 

clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj		<- StashIdent( seuratObj, save.name = paste0( umapDim, "D_UMAP_res_", myResolution))
return(seuratObj)

}
