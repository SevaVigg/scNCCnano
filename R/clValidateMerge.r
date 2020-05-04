clValidateMerge <- function( seuratObj, accuracy = 0.8, genes = 3){
#this snippet is used to merge clusters with the same marker genes using ValidateClusters function

source("R/getClusterTypes.r")

seuratObj 	<- ValidateClusters( seuratObj, pc.use = 1:6, top.genes = genes, min.connectivity = 0, acc.cutoff = accuracy, verbose = TRUE)
levels(seuratObj@ident) <- 1:length( levels( seuratObj@ident))
levels(seuratObj@ident) <- names( getClusterTypes( seuratObj))

return( seuratObj)

}
