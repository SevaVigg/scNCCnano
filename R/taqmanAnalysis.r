require( tidyverse )
require( viridis )
require( Seurat )
require( ComplexHeatmap)

plotDPI <- 600

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

taqmanPlotDir	<- file.path( plotDir, "taqman")
dir.create( taqmanPlotDir, showWarnings = FALSE)

 
taqmanPlotDir100dpi	<- file.path( taqmanPlotDir, "100dpi")
dir.create( taqmanPlotDir100dpi, showWarnings = FALSE)

taqmanPlotDir600dpi	<- file.path( taqmanPlotDir, "600dpi")
dir.create( taqmanPlotDir600dpi, showWarnings = FALSE)


plotHeight 	<- 10
plotWidth 	<- 10
Margin		<- 2

source("R/correctGeneNames.r")

taqmanData		<- read_delim( file = file.path( taqmanSourceDir, "30hpf_TAQMAN sum_paper.csv"), na = c("NA"), delim = ";")
colnames(taqmanData)	<- correctGeneNames( colnames(taqmanData))

taqmanGenes		<- setdiff( colnames( taqmanData), c("Sample", "rpl13", "Kan "))

taqmanData		<- as_tibble(sapply( taqmanData, function(x) ifelse( x %in% c("Undetermined", "undetermined", "Underemined", "underemined"), 40, x)))
taqmanData		<- type_convert( taqmanData)
taqmanData		<- rename( taqmanData, Kan = 'Kan ')

taqmanData 		<- mutate_at(taqmanData, setdiff( colnames( taqmanData), "Sample"), .funs = list( ~(40 - . )*log10(2)))
taqmanData		<- filter( taqmanData, Kan > quantile( apply(taqmanData[ , c("Kan")], 1, sum), 0.05))# re			
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


densPlot	<- ggplot( data = data.frame( dens$x, dens$y), mapping = aes( dens$x, dens$y)) + 
			geom_line() +
			geom_vline( xintercept = expLogMinimal - 1 , color = "red") +
			scale_x_continuous( name = "log10 probe expression" ) +
			scale_y_continuous( name = "Probe expression density") +
			theme(	plot.title 	= element_text( size = 50, hjust = 0),
				axis.title	= element_text( size = 30), 
				axis.text	= element_text( size = 30),
				axis.text.y	= element_text( margin = margin( l = 0, t = 0, r = 1, b = 0))
				) 
			#ggtitle( "taqman probes")


if (plotDPI == 600) {

  ggsave( paste0( "taqmanDensityPlot", ".png"), path = taqmanPlotDir600dpi, device = "png" , plot = densPlot, width = plotWidth, height = plotHeight, units = "cm", dpi = 600, scale = 4)

}

if (plotDPI == 100) {

  ggsave( paste0( "taqmanDensityPlot", ".png"), path = taqmanPlotDir600dpi, device = "png" , plot = densPlot, width = plotWidth, height = plotHeight, units = "cm", dpi = 100, scale = 4)

}


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
phox2bbPigment		<- intersect( genePos[["phox2bb"]], Pigment)
phox2bbPigmentRat	<- length( phox2bbPigment)/length( genePos[["phox2bb"]])

neurog1Pigment		<- intersect( genePos[["neurog1"]], Pigment)
neurog1PigmentRat	<- length( neurog1Pigment)/length( genePos[["neurog1"]])



clusterOrderedCells	<- colnames( seuratTaqman@data[ , order( t(seuratTaqman@data)[ ,"sox10"], t(seuratTaqman@data)[ ,"mitfa"], t(seuratTaqman@data)[, "neurog1"], t(seuratTaqman@data)[, "phox2bb"], t(seuratTaqman@data)[, "ltk"], t(seuratTaqman@data)[, "pax7b"], decreasing = TRUE )] )
clusterOrderedCells	<- intersect( clusterOrderedCells, genePos["ltk"])
	#numbers are easier to sort; we will need this for the heatmap
levels( seuratTaqman@ident)	<- names( getFinalClusterTypes( seuratTaqman))

heatMapHeight 	<- 7
heatMapWidth 	<- 10
Margin		<- 2

source("R/drawTaqmanHeatMap.r")

if (plotDPI == 600) {
png( file = file.path( taqmanPlotDir600dpi, "taqmanBiclustHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawTaqmanHeatMap( seuratTaqman,  clusterOrderedCells = colnames( seuratTaqman@data),  heatMapHeight, heatMapWidth, showCellNames = FALSE))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( taqmanPlotDir100dpi, "taqmanBiclustHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawTaqmanHeatMap( seuratTaqman,  clusterOrderedCells = colnames( seuratTaqman@data), heatMapHeight, heatMapWidth, showCellNames = FALSE))
dev.off()}

seuratTaqmanLtk <- SubsetData( seuratTaqman, cells.use = genePos$ltk)

if (plotDPI == 600) {
png( file = file.path( taqmanPlotDir600dpi, "taqmanLtkPosHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	draw( drawTaqmanHeatMap( seuratTaqmanLtk,  clusterOrderedCells = colnames( seuratTaqmanLtk@data),  heatMapHeight, heatMapWidth, showCellNames = FALSE))
dev.off()}

if (plotDPI == 100) {
png( file = file.path( taqmanPlotDir100dpi, "taqmanLtkPosHeatMap.png"),
	height = heatMapHeight + Margin + 1,
	width =  heatMapWidth + 2*Margin + 1,
	units = "in",
	res = plotDPI, 
	pointsize = 2 
)
	
	draw( drawTaqmanHeatMap( seuratTaqmanLtk,  clusterOrderedCells = colnames( seuratTaqmanLtk@data),  heatMapHeight, heatMapWidth, showCellNames = FALSE))

dev.off()}

png( filename = file.path( taqmanPlotDir, "taqmanGeneCors.png"))
	heatmap( cor(taqmanMatrixF), col = viridis(1024))
dev.off()

