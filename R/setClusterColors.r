setClusterColors <- function( seuratObj){

if(!require("colorspace")){
install.packages("colorspace")
library(colorspace)}

#This snippet sets a color code for clusters used in TSNEPlot

nClust		<- length( levels(seuratObj@ident))

clColors	<- integer(nClust)
names(clColors) <- levels( seuratObj@ident)
unassignedColors	<- grep("[0-9][0-9]*", names(clColors))

clColors[unassignedColors] <- sequential_hcl(n = length(unassignedColors), h1 = 250, h2 = 90, c1 = 40, c2 = 55, l1 = 50, l2 = 90, p1 = .5, p2 = 1.3)


clColors[ "I" ] 	<- "cyan"
clColors[ "M" ]   	<- "black"
clColors[ "E" ]  	<- "red"

names( clColors) <- NULL
return(clColors)
}
