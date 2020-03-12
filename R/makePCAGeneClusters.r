makePCAGeneClusters <- function( seuratObj, pcaDim, myResolution){
# This snippet performs clustering in seuratObj in gene and PCA space

require("Seurat")

source("R/getClusterTypes.r")
	
seuratObj	<- SetAllIdent( seuratObj, id = "originalCellTypes")

#gene Space based clusters
seuratObj	<- FindClusters( seuratObj, genes.use = rownames(seuratObj@data), k.param = 6, print.output = FALSE, force.recalc = TRUE, resolution = 1.4)
seuratObj 	<- BuildClusterTree( seuratObj, genes.use = rownames(seuratObj@data), do.plot = FALSE, do.reorder = TRUE) 

#BuildClusterTree names clusters accourding to their size, so we need to assign cluster types

clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj		<- StashIdent( seuratObj, save.name = 'geneClusters')

#pca based clusters

seuratObj 	<- FindClusters( seuratObj, reduction.type = 'pca', dims.use = 1:pcaDim, k.param = 5, 
			print.output = FALSE, force.recalc = TRUE, resolution = myResolution, save.SNN = TRUE)
seuratObj 	<- BuildClusterTree( seuratObj, pcs.use = 1:pcaDim, do.plot = FALSE, do.reorder = TRUE)
seuratObj	<- ValidateClusters( seuratObj, pc.use = 1:pcaDim, top.genes = 3, acc.cutoff = 0.8, min.connectivity = 0) 
#BuildClusterTree names clusters accourding to their size, so we need to reassign cluster types

clTypes 		<- getClusterTypes(seuratObj)
levels(seuratObj@ident) <- names(clTypes)
seuratObj		<- StashIdent( seuratObj, save.name = paste0( pcaDim, "D_PCA_res_", myResolution))

return( seuratObj)
}
