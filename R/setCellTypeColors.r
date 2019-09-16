#This snippet assignes cell colors accordingly to their cell types
setCellTypeColors <- function(seuratObj, testCellType = "G"){


if(!require("colorspace")){
install.packages("colorspace")
library(colorspace)}


#basic colors are shades of green, according to the clustering ident


cellColors		<- integer( length(levels( seuratObj@ident)))
names(cellColors)	<- levels( seuratObj@ident) 
standColors		<- grep( "[0-9][0-9]", names(cellColors), value = TRUE)
if (length(standColors)) cellColors[standColors] <- sequential_hcl( length( standColors) , "BluGrn", alpha = 0.5, rev = TRUE)

cellColors[ "I" ] 	= "cyan"
cellColors[ "M" ]   	= "black"
cellColors[ "sox10-" ]  = "gold"
cellColors[ "mitfa-" ]  = "dodgerblue"
cellColors[ "Tl" ]  	= "red"
cellColors[ testCellType  ]	= "green"

#make tails red according to their "snail2" expression

return(cellColors)

}
