makePCAClusters	<- function( seuratObj, PCADim, myResolution){

require("Seurat")
source("R/getClusterTypes.r")

seuratObj 		<- FindClusters( seuratObj, reduction.type = 'pca', 
					dims.use = 1:PCADim, 
					k.param = 15, print.output = FALSE, force.recalc = TRUE, save.SNN = TRUE,
					resolution = myResolution) 
#seuratObj 		<- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = TRUE) 
clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
#seuratObj		<- ValidateClusters( seuratObj, pc.use = 1:pcaDim, top.genes = 4, acc.cutoff = 0.82, min.connectivity = 0) 
#seuratObj 		<- BuildClusterTree( seuratObj, pcs.use = 1:pcaDim, do.plot = FALSE, do.reorder = TRUE)
clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj		<- StashIdent( seuratObj, save.name = paste0( PCADim, "D_PCA_res_", myResolution))
return(seuratObj)

}
