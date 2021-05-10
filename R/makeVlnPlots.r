makeVlnPlots	<- function( seuratObj, name = "vlnPlots", orientation = "portrait"){

source("R/setClusterColors.r")

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

vlnPlotDir <- file.path( plotDir, "vlnPlots")
dir.create( vlnPlotDir, showWarnings = FALSE)


if (orientation == "landscape") { pageWidth = 18; pageHeight = 13 }
if (orientation == "portrait") { pageWidth = 18; pageHeight = 13 }

if( seuratObj@project.name == "taqman"){ 
	
vlnPlot <- VlnPlot( seuratObj, c( 	"sox9b"		, "snail2"	, "sox10"	, "pax7b"	, "mbpa"	,
					"phox2b"	, "mitfa"	, "ltk"		, "neurog1"	, "pnp4a"	,
					 "tyrp1b"	, "xdh"		, "elavl3"),
			nCol = 5, cols.use = setClusterColors( seuratObj), do.return = TRUE) 

}else{

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
}

if( plotDPI == 600){
	
  vlnPlotDir600dpi <- file.path( vlnPlotDir, "600dpi")
  dir.create( vlnPlotDir600dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = vlnPlotDir600dpi, device = "png" , plot = vlnPlot, width = pageWidth, height = pageHeight, units = "cm", dpi = 600, scale = 4)

}else if( plotDPI == 100){

  vlnPlotDir100dpi <- file.path( vlnPlotDir, "100dpi")
  dir.create( vlnPlotDir100dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = vlnPlotDir100dpi, device = "png" , plot = vlnPlot, width = pageWidth, height = pageHeight, units = "cm", dpi = 100, scale = 4)

}else{ cat("Select plotDPI 600 or plotDPI 100")}


ggsave( paste0( name, ".png"), path = vlnPlotDir, device = "png" , plot = vlnPlot, width = pageWidth, height = pageHeight, units = "cm", dpi = 600, scale = 5) 

}

