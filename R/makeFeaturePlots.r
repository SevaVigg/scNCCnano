makeFeaturePlots	<- function( seuratObj, minCutoff, dimRed){

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

featurePlotDir <- file.path( plotDir, "featurePlots")
dir.create( plotDir, showWarnings = FALSE)


png( file.path( featurePlotDir, paste0("featurePlots_", minCutoff, ".png")), height = 2048, width = 1536)
		FeaturePlot( seuratObj, c( 	"tfap2e", "her9", "mc1r",   "foxg1b", "snail2", 
						"sox9b",  "hmx1", "tfap2a", "hmx4",   "foxo1b",
						"mbpa",   "otx2",  "foxd3", "sox10",  "kita", 
						"mitfa",   "tyr",  "tfec",   "ltk", "foxo1a",
					        "tyrp1b",  "slc24a5", "oca2", "mlphb", "silva",
					        "pnp4a",   "ednrba", "alx4b", "impdh1b", "id2a",
					        "pax7a",   "pax7b", "myo5aa", "foxp4", "smad9"),    
	nCol = 5, do.return = FALSE, pt.size = 3, min.cutoff = minCutoff, cols.use = c("blue", "red"), reduction.use = dimRed)   
dev.off()

}

