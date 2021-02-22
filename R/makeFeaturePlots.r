makeFeaturePlots	<- function( seuratObj, minCutoff, dimRed, name = "featurePlots_"){

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

featurePlotDir <- file.path( plotDir, "featurePlots")
dir.create( featurePlotDir, showWarnings = FALSE)


png( file.path( featurePlotDir, paste0( name, minCutoff, ".png")), height = 2048, width = 1536)
		FeaturePlot( seuratObj, 
				c( 	"tfap2e"	, "tfap2a"	, "her9"	, "sox9b"	, "foxg1b" 	,
					"snail2"	, "alx4b"	, "hmx1"	, "otx2"	, "sox10"	, 
					"impdh1b"	, "foxo1b"	, "tyr"		, "pax7b"	, "mc1r"	,						
					"id2a"		, "hmx4"	, "foxd3"	, "ednrba"	, "kita"	,
					"mbpa"		, "phox2b"	, "pax7b"	, "tfec"	, "mitfa"	,
					"foxp4" 	, "foxo1a"	, "hbp1"	, "ltk"		, "mycl1a"	, 
					"pax7a"		, "tyrp1b"	, "slc24a5"	, "oca2"	, "mlphb"	, 
					"silva"		, "myo5aa"	, "pnp4a"	, "ets1a"	, "fgfr3_v2"	),
,
	nCol = 5, do.return = FALSE, pt.size = 3, min.cutoff = minCutoff, cols.use = c("blue", "red"), reduction.use = dimRed)   
dev.off()

}

