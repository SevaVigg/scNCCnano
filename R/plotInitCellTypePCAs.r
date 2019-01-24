plotInitCellTypePCAs <- function( seuratObj, Ncomps, plotDir){
source("R/setCellTypeColors.r")

clIdent <- seuratObj@ident
cellColors <- setCellTypeColors(ipmc)

PCAPlotDirName <- file.path( plotDir, "PCAPlots")
dir.create( PCAPlotDirName, showWarnings = FALSE)

for (c1 in 1:(Ncomps-1)){
	for (c2 in (c1+1):Ncomps){
		png( file.path( PCAPlotDirName, paste0( "PCA_", c1, "_", c2, ".png")))
			PCAPlot( seuratObj, dim.1 = c1, dim.2 = c2, cols.use = cellColors)
		dev.off()
	}

}
   
}


