setClusterColors <- function( seuratObj){

if(!require("colorspace")){
install.packages("colorspace")
library(colorspace)}

#This snippet sets a color code for clusters used in TSNEPlot

nClust		<- length( levels(seuratObj@ident))

clColors	<- integer(nClust)
names(clColors) <- levels( seuratObj@ident)

numericColors	<- grep("^[0-9]+" , names(clColors))

clColors[ which( "Tl" == names(clColors))]		<- "red"
clColors[ which( "eHMP" == names(clColors))]		<- "red"
clColors[ which( "X" == names(clColors))]		<- "gold"
clColors[ which( "M" == names(clColors))]		<- "black"
clColors[ which( "I" == names(clColors))]		<- "cyan"
clColors[ which( "sox10" == names(clColors))]		<- "mediumorchid3"
clColors[ which( "G" == names(clColors))]		<- "yellowgreen"
clColors[ which( "R" == names(clColors))]		<- "yellowgreen"
clColors[ which( "ltHMP" == names(clColors))]		<- "magenta"

if ( length(numericColors)) clColors[numericColors] <- sequential_hcl(n = length(numericColors), h1 = 250, h2 = 90, c1 = 40, c2 = 55, l1 = 50, l2 = 90, p1 = .5, p2 = 1.3)

names( clColors) <- NULL
return(clColors)
}
