renameSeuratGenes <- function( seuratObj){

source("R/correctGeneNames.r")

seuratNew 	<- seuratObj
allGenes	<- rownames( seuratObj@data)

rownames( seuratNew@raw.data) 		<- correctGeneNames( allGenes)
rownames( seuratNew@data) 		<- correctGeneNames( allGenes)
rownames( seuratNew@scale.data) 	<- correctGeneNames( allGenes)
rownames( seuratNew@dr$pca@gene.loadings)<- correctGeneNames( allGenes) 

return( seuratNew)

}
