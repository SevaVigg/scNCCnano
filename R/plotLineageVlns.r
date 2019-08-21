plotLineageVlns <- function(seuratObj, slingObj, gene, lineageId){

source("R/setClusterColors.r")

lineage 	<- slingObj@lineages[[lineageId]]

#lind 		<- data.frame(l = lineage, sl = sort( as.numeric (lineage)), stringsAsFactors = FALSE) # create a transposition
#colorInd 	<- sapply( lind$l, function(x) which( lind$sl ==x)  )

vln <- VlnPlot(seuratObj, gene, ident.include = slingObj@lineages[[lineageId]], do.return = TRUE)
clColors	<- setClusterColors( SubsetData(seuratObj, ident.use = lineage))


plot(vln 
		+ scale_x_discrete(limits = slingObj@lineages[[lineageId]])
		+ scale_fill_manual(values = clColors)
		)

}
