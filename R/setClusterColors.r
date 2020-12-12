setClusterColors <- function( seuratObj){

if(!require("colorspace")){
install.packages("colorspace")
library(colorspace)}

#This snippet sets a color code for clusters used in TSNEPlot

nClust		<- length( levels(seuratObj@ident))

clColors	<- integer(nClust)
names(clColors) <- levels( seuratObj@ident)

numericColors	<- grep("^[0-9]+", names(clColors))

clColors[ grep( "Tl", names(clColors))]		<- "red"
clColors[ grep( "eNCC", names(clColors))]	<- "red"
clColors[ grep( "X", names(clColors))]		<- "gold"
clColors[ grep( "M", names(clColors))]		<- "black"
clColors[ grep( "I", names(clColors))]		<- "cyan"
clColors[ grep( "sox", names(clColors))]	<- "mediumorchid3"
clColors[ grep( "G", names(clColors))]		<- "yellowgreen"
clColors[ grep( "R", names(clColors))]		<- "yellowgreen"
clColors[ grep( "HMP", names(clColors))]	<- "magenta"

if ( length(numericColors)) clColors[numericColors] <- sequential_hcl(n = length(numericColors), h1 = 250, h2 = 90, c1 = 40, c2 = 55, l1 = 50, l2 = 90, p1 = .5, p2 = 1.3)

names( clColors) <- NULL
return(clColors)
}
