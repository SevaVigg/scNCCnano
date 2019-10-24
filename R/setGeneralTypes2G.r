# Requiers seuratObj created by seuratNorm.r
#
#this snippet is use to make TSNE, UMAP and PCA plots for initial cell types and clusters in the gene space
#it requires to run first seuratNorm.r, which among others sets the appropriate experimentType ( allCells, WT, WT_sox10)
#in seuratObj@project.name
#
# Directory structure :     Res <- Plots <- geneSpace

setGeneralTypes2G <- function( seuratObj){

source("R/getClusterTypes.r")
source("R/calcTSNEGeneSpace.r")
source("R/setCellTypeColors.r")
source("R/setClusterColors.r")
source("R/calcUMAPGeneSpace.r")

if( length(setdiff( c("18", "21", "24", "30", "36", "48", "60", "72"), levels( seuratObj@ident)))!=0) { cat("The seuratObj does not contain general types"); return(seuratObj)}

seuratObj 	<- StashIdent( seuratObj, save.name = "generalIdents")
levels(seuratObj@ident) <- c(levels(seuratObj@ident), "G")
seuratObj@ident[ grep("general", names(seuratObj@ident))] <- "G"
seuratObj@ident <- droplevels(seuratObj@ident)

return(seuratObj)
}
