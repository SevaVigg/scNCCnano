setClusterColors <- function( clustType){

if(!require("colorspace")){
install.packages("colorspace")
library(colorspace)}

#This snippet sets a color code for clusters used in TSNEPlot

nClust		<- length( clustType)

clColors	<- integer(nClust)
names(clColors) <- names( clustType)
clColors[grep("[0-9][0-9]*", names(clustType))] <- qualitative_hcl( nClust-3, "Dynamic", alpha = 1)


clColors[ "I" ] 	<- "cyan"
clColors[ "M" ]   	<- "black"
clColors[ "Tl" ]  	<- "red"

names( clColors) <- NULL
return(clColors)
}
