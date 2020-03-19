clValidateMerge <- function( seuratObj){
#this snippet is used to merge clusters with the same marker genes using ValidateClusters function

source("R/getClusterTypes.r")

seuratObj 	<- ValidateClusters( seuratVal, pc.use = 1:6, top.genes = 3, min.connectivity = 0, acc.cutoff = 0.85, verbose = TRUE)
levels(seuratObj@ident) <- 1:length( levels( seuratObj@ident))
levels(seuratObj@ident) <- names( getClusterTypes( seuratObj))

return( seuratObj)

}
