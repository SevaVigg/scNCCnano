makeFeaturePlots	<- function( seuratObj, minCutoff, dimRed, name = "featurePlots_"){

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

featurePlotDir <- file.path( plotDir, "featurePlots")
dir.create( plotDir, showWarnings = FALSE)


png( file.path( featurePlotDir, paste0( name, minCutoff, ".png")), height = 2048, width = 1536)
		FeaturePlot( seuratObj, c( 	"tfap2e", "sox9b", "tfap2a", "snail2","sox10", 
						"foxo1b", "foxg1b", "pax7b", "otx2",  "hmx4", 
						"mc1r",    "mbpa",  "alx4b", "foxd3", "impdh1b", 
						"id2a",   "her9", "ednrba","kita", "mitfa", 
						"tyr",    "tyrp1b", "slc24a5", "oca2", "mlphb", 
						"silva", "phox2b", "ltk",      "foxo1a", "tfec", 
						"hbp1", "ets1a", "mycl1a", "gapdh", "fgfr3_v2", 
						"pnp4a", "pax7a", "myo5aa","foxp4", "smad9"),
	nCol = 5, do.return = FALSE, pt.size = 3, min.cutoff = minCutoff, cols.use = c("blue", "red"), reduction.use = dimRed)   
dev.off()

}

