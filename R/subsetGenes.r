subsetGenes 		<- function( seuratObj, genesSubset){

subsetMatrix		<- seuratObj@data[ genesSubset, ]	
seuratGeneSubset   	<- CreateSeuratObject( raw.data = subsetMatrix)

origIdent 		<- seuratObj@ident
names(origIdent) 	<- colnames(ipmc@data)

seuratGeneSubset	<- AddMetaData( object = seuratGeneSubset, metadata = origIdent, col.name = "orig.ident") # Add the idents to the meta.data slot
seuratGeneSubset	<- SetAllIdent( object = seuratGeneSubset, id = "orig.ident") # Assign identities for the new Seurat object

seuratGeneSubset@misc	<- "subsetGenes"

seuratGeneSubset	<- ScaleData( seuratGeneSubset, do.scale = TRUE, do.center = TRUE)
seuratGeneSubset 	<- RunPCA( seuratGeneSubset, pc.genes = rownames(seuratGeneSubset@data), weight.by.var = FALSE, do.print = FALSE)
seuratGeneSubset	<- StashIdent(object = seuratGeneSubset, save.name = "originalCellTypes")

seuratGeneSubset@project.name <- seuratObj@project.name

return( seuratGeneSubset)

}
