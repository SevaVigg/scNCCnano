makeFeaturePlots	<- function( seuratObj, dimRed){

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

experimentTypePlotDir <- file.path(plotDir, seuratObj@project.name)
dir.create( experimentTypePlotDir, showWarnings = FALSE)

geneSpacePlotDir <- file.path( experimentTypePlotDir, "geneSpacePlots")
dir.create( geneSpacePlotDir, showWarnings = FALSE)


png( file.path( geneSpacePlotDir, paste0("geneDistributioni_", dimRed, ".png")), height = 1536, width = 2048)
	FeaturePlot( seuratObj, c( "sox10", "sox9b", "tfap2e", "snail2", "phox2b",  
			 	"pax7a", "pax7b",  "tfap2a", "foxd3", "impdh1b",
				"tfec", "mitfa", "ltk", "mbpa", "pnp4a", 
				"tyrp1b", "mlphb", "oca2", "silva", "slc24a5"), nCol = 5, 
	  	do.return = FALSE, pt.size = 3, cols.use = c("blue", "red"), reduction.use = dimRed)   

dev.off()

}

