makeVlnPlots	<- function( seuratObj,  name = "vlnPlots"){

source("R/setClusterColors.r")

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

vlnPlotDir <- file.path( plotDir, "vlnPlots")
dir.create( vlnPlotDir, showWarnings = FALSE)


		vlnPlot <- VlnPlot( seuratObj, c( "tfap2e", "tfap2a", "sox9b", "snail2","sox10", 
						"ednrba", "foxd3", "id2a", "impdh1b", "otx2",						
						"alx4b", "foxg1b", "her9", "hmx1", "hmx4",
						"kita", "mbpa", "mc1r", "pax7b", "mitfa",						
						"tyr",   "tyrp1b", "slc24a5", "oca2", "mlphb", 
						"silva", "phox2b", "ltk",  "foxo1a", "tfec", 
						"hbp1", "ets1a", "mycl1a", "foxo1b", "fgfr3_v2", 
						"pnp4a", "pax7a", "myo5aa","foxp4", "smad9"),
	nCol = 5,  cols.use = setClusterColors(seuratObj), do.return = TRUE)

ggsave( paste0( name, ".png"), path = vlnPlotDir, device = "png" , plot = vlnPlot, width = 13, height = 18, units = "cm", dpi = 300, scale = 5) 

}

