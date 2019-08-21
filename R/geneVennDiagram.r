# This script constructs gene Venn diagram from gene expression data
# expMatrix is the matrix of gene log expressions with expression values in columns
# geneList is the list of two, three, or four genes


geneVennDiagram <- function( geneList, expMatrix){

if(!require("VennDiagram")){
  install.packages("VennDiagram")
}

require(VennDiagram)
require(gtools)

nGenes	<- length( geneList)

dens		<- density( t(expMatrix))
expThreshold	<- optimize(approxfun(dens$x,dens$y),interval=c(5,15))$minimum 

geneIndexList 	<- sapply( geneList,  function(gene) which(expMatrix[ gene, ] > expThreshold))
geneIndexNums 	<- lapply(geneIndexList, length)
combs		<- if(nGenes == 3) data.frame( "C1" = geneList[c(1,2)], "C2" = geneList[c(2,3)], "C3" = geneList[c(1,3)], 
		stringsAsFactors = FALSE) else data.frame(geneList[c(1,2)])				#required by draw.triple.venn - which requires 12, 23, 13 
intersectIndex	<- lapply( combs, function( genePair) do.call( intersect, unname(geneIndexList[ genePair])))
intersectIndexNums <- unname(lapply( intersectIndex, length)) 



if (nGenes == 2){ 
	do.call( draw.pairwise.venn, c( unname( geneIndexNums), intersectIndexNums, list(scaled = TRUE), list(category = names( geneIndexList))))
	}else{
		if (nGenes == 3){
		do.call( draw.triple.venn,  c( unname( geneIndexNums), 
						intersectIndexNums, 
						length(Reduce(intersect, geneIndexList)), 
						list( scaled = FALSE,
						category = names( geneIndexList),
						fill = c("yellow", "blue", "green"),
						cat.cex = 1.8,
						euler.d = FALSE)
						)
			)
		}else{ 
			if (nGenes == 4){
				}else{ 
				return("incorrect number of genes")}
			}
		}	
}





