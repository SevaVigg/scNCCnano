makeUmapClusters	<- function( seuratObj, umapDim, myResolution){

require("Seurat")
source("R/getClusterTypes.r")

pcaDim			<- seuratObj@calc.params$RunPCA$pcs.compute
seuratObj 		<- FindClusters( seuratObj, reduction.type = 'umap', 
					dims.use = 1:umapDim, 
					k.param = 15, print.output = FALSE, force.recalc = TRUE, save.SNN = TRUE,
					resolution = myResolution) 
#seuratObj 		<- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = TRUE) 
clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
#seuratObj		<- ValidateClusters( seuratObj, pc.use = 1:pcaDim, top.genes = 3, acc.cutoff = 0.8, min.connectivity = 0) 
#seuratObj 		<- BuildClusterTree( seuratObj, pcs.use = 1:pcaDim, do.plot = FALSE, do.reorder = TRUE)
#clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj		<- StashIdent( seuratObj, save.name = paste0( umapDim, "D_UMAP_res_", myResolution))
return(seuratObj)

}
