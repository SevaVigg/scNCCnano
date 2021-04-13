makeVlnPlots	<- function( seuratObj,  name = "vlnPlots"){

source("R/setClusterColors.r")

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

vlnPlotDir <- file.path( plotDir, "vlnPlots")
dir.create( vlnPlotDir, showWarnings = FALSE)


vlnPlot <- VlnPlot( seuratObj, c( 	"tfap2e"	, "tfap2a"	, "her9"	, "sox9b"	, "foxg1b" 	,
					"snail2"	, "alx4b"	, "hmx1"	, "otx2"	, "sox10"	, 
					"impdh1b"	, "foxo1b"	, "tyr"		, "pax7b"	, "mc1r"	, 
					"id2a"		, "hmx4"	, "foxd3"	, "ednrba"	, "kita"	, 
					"mbpa"		, "phox2b"	, "tfec"	, "mitfa"	, "foxp4" 	, 
					"foxo1a"	, "hbp1"	, "ltk"		, "mycl1a"	, "pax7a"	, 
					"tyrp1b"	, "slc24a5"	, "oca2"	, "mlphb"	, "silva"	, 
					"myo5aa"	, "pnp4a"	, "ets1a"	, "fgfr3_v2"	, "pax3_v2"	, 
					"dpf3"		, "smad9"),

	nCol = 5,  cols.use = setClusterColors(seuratObj), do.return = TRUE)

ggsave( paste0( name, ".png"), path = vlnPlotDir, device = "png" , plot = vlnPlot, width = 13, height = 18, units = "cm", dpi = 600, scale = 5) 

}

