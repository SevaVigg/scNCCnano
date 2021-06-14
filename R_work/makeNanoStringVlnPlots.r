makeVlnPlots	<- function( seuratObj, name = "vlnPlotsNanostring"){

source("R/setClusterColors.r")

orientation <- "landscape"
#orientation <- "portrait"

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

vlnPlotDir <- file.path( plotDir, "vlnPlots")
dir.create( vlnPlotDir, showWarnings = FALSE)


vlnPlot <- VlnPlot( seuratObj, c( 	"tfap2e"	, "tfap2a"	, "her9"	, "sox9b"	, "foxg1b" 	,
					"snai1b"	, "alx4b"	, "hmx1"	, "otx2b"	, "sox10"	, 
					"impdh1b"	, "foxo1b"	, "tyr"		, "pax7b"	, "mc1r"	, 
					"id2a"		, "hmx4"	, "foxd3"	, "ednrba"	, "kita"	, 
					"mbpa"		, "phox2bb"	, "tfec"	, "mitfa"	, "foxp4" 	, 
					"foxo1a"	, "hbp1"	, "ltk"		, "mycla"	, "pax7a"	, 
					"tyrp1b"	, "slc24a5"	, "oca2"	, "mlphb"	, "pmela"	, 
					"myo5aa"	, "pnp4a"	, "ets1"	, "fgfr3"	, "pax3a"	, 
					"dpf3"		, "smad9"),

	nCol = 5,  cols.use = setClusterColors(seuratObj), do.return = TRUE)

if (orientation == "landscape") { pageWidth = 18; pageHeight = 13 }
if (orientation == "portrait") { pageWidth = 18; pageHeight = 13 }


ggsave( paste0( name, ".png"), path = vlnPlotDir, device = "png" , plot = vlnPlot, width = pageWidth, height = pageHeight, units = "cm", dpi = 600, scale = 5) 

}

