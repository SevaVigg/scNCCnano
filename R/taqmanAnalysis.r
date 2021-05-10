require( tidyverse )
require( viridis )
require( Seurat )
require( ComplexHeatmap)

workDir <- getwd()

resDir		 <- file.path( workDir, "Res")
dir.create( resDir, showWarnings = FALSE)

scTablesDir	<- file.path( resDir, "scTables")
dir.create( scTablesDir, showWarnings = FALSE)

taqmanSourceDir	 	<- file.path( workDir, "SourceData", "Taqman")

taqmanResDir	<- file.path( resDir, "taqman")
dir.create( taqmanResDir, showWarnings = FALSE)

plotDir 	<- file.path( resDir, "Plots")
dir.create( plotDir, showWarnings = FALSE)  

taqmanData		<- read_delim( file = file.path( taqmanSourceDir, "30hpf_TAQMAN sum_paper.csv"), na = c("NA"), delim = ";")
taqmanGenes		<- setdiff( colnames( taqmanData), c("Sample", "rpl13", "Kan "))

taqmanData		<- as_tibble(sapply( taqmanData, function(x) ifelse( x %in% c("Undetermined", "undetermined", "Underemined", "underemined"), 40, x)))
taqmanData		<- type_convert( taqmanData)
taqmanData		<- rename( taqmanData, Kan = 'Kan ')

taqmanData 		<- mutate_at(taqmanData, setdiff( colnames( taqmanData), "Sample"), .funs = list( ~(40 - . )*log10(2)))
taqmanData		<- filter( taqmanData, rpl13+Kan > quantile( apply(taqmanData[ , c("rpl13", "Kan")], 1, sum), 0.05))# re			
taqmanData_Norm		<- mutate_at(taqmanData, .vars = vars(taqmanGenes), .funs = list( ~ .-Kan)) #we work in the log domain
#taqmanData_Norm		<- mutate_at(taqmanData, .vars = vars(taqmanGenes), .funs = list( ~ .-0.5*(Kan + rpl13))) #we work in the log domain
#taqmanData_Norm	<- taqmanData

taqmanMatrix	<- as.matrix(taqmanData_Norm[ , taqmanGenes])
rownames(taqmanMatrix) <- taqmanData$Sample

taqmanMatrix	<- apply( taqmanMatrix, 1:2, function(x) ifelse( is.na(x), 0, x))

keep <- apply(taqmanMatrix, 1, sum) > quantile( apply(taqmanMatrix, 1, sum), 0.05) 
taqmanMatrixF <- taqmanMatrix[keep, ]

taqmanMatrixF <- taqmanMatrixF - min( taqmanMatrixF)

dens			<- density(as.matrix( taqmanMatrixF)) 
expLogMinimal		<- optimize(approxfun(dens$x,dens$y),interval=c(2,8))$minimum
taqmanMatrixF 		<- apply( taqmanMatrixF, 1:2, function(x) ifelse( x > expLogMinimal-1, x, 0)) 


#now normalize for the minimum - which means substract the minimal value in the log scale


taqmanScTableDir 	<- file.path( scTablesDir, "taqman")
dir.create( taqmanScTableDir , showWarnings = FALSE)

write.table( taqmanMatrixF, file = file.path( taqmanScTableDir, "taqmanLogExps.csv"), sep = "\t")

source("R/seuratNorm.r")
source("R/getFinalClusterTypes.r")

clusterTreeDir		<- file.path( plotDir, "clusterTrees")
dir.create( clusterTreeDir, showWarnings = FALSE)

seuratTaqman 		<- seuratNorm( "taqman" )
genePos			<- sapply( taqmanGenes, function(x) WhichCells(seuratTaqman, subset.name = x, accept.low = 0))

Pigment			<- union( union( genePos[["ltk"]], genePos[["mitfa"]]), genePos[["pax7b"]])
phox2bPigment		<- intersect( genePos[["phox2b"]], Pigment)
phox2bPigmentRat	<- length( phox2bPigment)/length( genePos[["phox2b"]])

neurog1Pigment		<- intersect( genePos[["neurog1"]], Pigment)
neurog1PigmentRat	<- length( neurog1Pigment)/length( genePos[["neurog1"]])



seuratTaqman 		<- FindClusters( seuratTaqman, dims.use = 1:10, k.param = 4, prune.SNN = 0, save.SNN = TRUE, resolution = 0.8)
seuratTaqman 		<- BuildClusterTree( seuratTaqman, pcs.use = 1:10, do.reorder = TRUE, reorder.numeric = TRUE)
clusterOrderedCells	<- colnames( seuratTaqman@data[ , order( t(seuratTaqman@data)[ ,"sox10"], t(seuratTaqman@data)[ ,"mitfa"], t(seuratTaqman@data)[, "neurog1"], t(seuratTaqman@data)[, "phox2b"], t(seuratTaqman@data)[, "ltk"], t(seuratTaqman@data)[, "pax7b"], decreasing = TRUE )] )
clusterOrderedCells	<- intersect( clusterOrderedCells, ltkPos)
	#numbers are easier to sort; we will need this for the heatmap
levels( seuratTaqman@ident)	<- names( getFinalClusterTypes( seuratTaqman))
seuratTaqman		<- BuildClusterTree( seuratTaqman, pcs.use = 1:10, do.reorder = FALSE)

png( file.path( clusterTreeDir, "taqmanClusterTree.png"))
	PlotClusterTree( seuratTaqman) 
dev.off()

source("R/makeVlnPlots.r")
#makeVlnPlots( seuratTaqman, name = "taqmanVlnPlots")

source("R/calcUmapGeneSpace.r")
#seuratTaqman <- calcUmapGeneSpace( seuratTaqman, minDist = 0.1, myNeighbors = 15L)$All
   
source("R/makeFeaturePlots.r")
#makeFeaturePlots( seuratTaqman, 0.5, "umap", name = "taqmanFeaturePlots_")

heatMapDir <- file.path( plotDir, "heatMaps")
dir.create( heatMapDir, showWarnings = FALSE)

source("R/drawTaqmanHeatMap.r")

png( file.path( heatMapDir, "taqmanClustHeatMap.png"), width = 800, height = 600)
	draw( drawTaqmanHeatMap( seuratTaqman, clusterOrderedCells), annotation_legend_side = "bottom")  
dev.off()

stop()
png( filename = file.path( taqmanPlotDir, "taqmanHeatMap.png"))
	heatmap( taqmanMatrixF, col = viridis(1024))
dev.off()

png( filename = file.path( taqmanPlotDir, "taqmanGeneCors.png"))
	heatmap( cor(taqmanMatrixF), col = viridis(1024))
dev.off()

ltkPos		<- which(taqmanMatrixF[ , "ltk"] > 0.5)
ltkPosGenes	<- head( sort ( apply(taqmanMatrixF[ ltkPos, ], 2, mean), decreasing = TRUE), 5)

mB		<- which(taqmanMatrixF[ , "tyrp1b"] > 0.5)
mBGenes		<- head( sort ( apply(taqmanMatrixF[ mB, ], 2, mean), decreasing = TRUE), 5)

eHMP		<- which(taqmanMatrixF[ , "sox9b"] > 0.5)
eHMPGenes	<- head( sort ( apply(taqmanMatrixF[ eHMP, ], 2, mean), decreasing = TRUE), 5)

xP		<- which(taqmanMatrixF[ , "xdh"] > 0.5)
xPGenes		<- head( sort ( apply(taqmanMatrixF[ xP, ], 2, mean), decreasing = TRUE), 5)

MIX		<- which(taqmanMatrixF[ , "pnp4a"] > 0.5)
MIXGenes	<- head( sort ( apply(taqmanMatrixF[ MIX, ], 2, mean), decreasing = TRUE), 5)

neuro		<- which(taqmanMatrixF[ , "neurog1"] > 0.5)
neuroGenes	<- head( sort ( apply(taqmanMatrixF[ neuro, ], 2, mean), decreasing = TRUE), 5)

phox2bPos	<- which(taqmanMatrixF[ , "phox2b"] > 0.5)
neuroGenes	<- head( sort ( apply(taqmanMatrixF[ phox2bPos, ], 2, mean), decreasing = TRUE), 5)


seuratTaqman 	<- CreateSeuratObject( raw.data <- t(taqmanMatrixF), is.expr = log10(5))
seuratTaqman	<- ScaleData( seuratTaqman, do.scale = TRUE, do.center = TRUE)
seuratTaqman 	<- RunPCA( seuratTaqman, pcs.compute = 10, pc.genes = rownames(seuratTaqman@data), weight.by.var = FALSE, do.print = FALSE)

seuratTaqman 	<- FindClusters( seuratTaqman, dims.use = 1:10, k.param = 4, prune.SNN = 0, resolution = 0.8)











#normCoeff	<- apply( taqmanData %>% select( c("Kan ", "rpl13")), 1, function(x) exp(mean(log(x))))
#normData	<- apply( taqmanData %>% select( taqmanGenes), 2, function(x) x/normCoeff)



		

