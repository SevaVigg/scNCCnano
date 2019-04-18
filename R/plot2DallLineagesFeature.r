plot2DallLineagesFeature <- function( slingShotObj, seuratObj, gene, dimRed){

#seuratObject must contain results of clustering and tSNE pre-calculated

source("R/setClusterColors.r")
source("R/getLineageCoords.r")
source("R/dimPlot3DFeature.r")
source("R/getLineageCoords.r")

LineageTree	<- getLineageCoords( seuratObj, slingShotObj, dimRed)  

plotVals 	<- seuratObj@dr[[ dimRed ]]@cell.embeddings

featureValues <- seuratObj@data[ gene, ]
cellColors <- cut( featureValues/max(featureValues),
breaks= seq(0, 1, 0.1))

colors 			<- sequential_hcl(11, palette = "YlGnBu")
levels(cellColors)	<- colors
#cellColors <- setClusterColors( seuratObj)[ seuratObj@ident]


par( bg = "darkblue", fg = "white", col.axis = "white", col.lab = "white")
plot( plotVals[,2]~ plotVals[,1], cex = 0.7, pch = 16, col = cellColors, xlab = paste0(dimRed, "_1"), ylab = paste0(dimRed, "_2"))

sapply(LineageTree, function(x) {lines(t(x), col = "white"); points(t(x), pch = 16, col = "white")})

}


