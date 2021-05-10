makeFeaturePlots	<- function( seuratObj, minCutoff, dimRed, plotDPI, orientation = "portrait", name = "featurePlots_"){

require( cowplot)
require( Seurat)
source("R/setClusterColors.r")

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

featurePlotDir 	<- file.path( plotDir, "featurePlots")
dir.create( featurePlotDir, showWarnings = FALSE)
 
if (orientation == "landscape") { pageWidth = 11.7; pageHeight = 8.3 }
if (orientation == "portrait") { pageWidth = 8.3; pageHeight = 11.7 }

clusterPlot	    <- DimPlot( seuratObj, reduction.use = "umap", pt.size = 1, cols.use = setClusterColors( seuratObj), no.legend = TRUE, label.size = 0.5) +
			 	theme( plot.title = element_text( size = 12, face = "bold"), 
					axis.title.y = element_text( margin = margin( t = 0, r = 0, b = 0, l = 10)),  
					axis.title.x = element_text( margin = margin( t = 0, r = 0, b = 10, l = 0))) + 
				ggtitle( "clusters")  
source("R/setClusterColors.r")

if( seuratObj@project.name == "taqman"){ 
	
	featurePlotList <- FeaturePlot( seuratObj, c( 	"sox9b"		, "snail2"	, "sox10"	, "pax7b"	, "mbpa"	,
							"phox2b"	, "mitfa"	, "ltk"		, "neurog1"	, "pnp4a"	,
					 		"tyrp1b"	, "xdh"		, "elavl3"),
			nCol = 5, pt.size = 2, cols.use = c("blue", "red"), reduction.use = dimRed, do.return = TRUE) 

	featurePlotList$clusters <- clusterPlot
	featurePlotList <- lapply( featurePlotList, function(cPlot) cPlot + theme( 	plot.margin = unit( c( 0, 0, 0, 0), "cm"),
										plot.title  = element_text( size = 12), 
										axis.title  = element_text( size = 8),
										axis.text.x = element_text(
											angle = 90,
											hjust = 0.1,
											vjust = 0.1,
											size =	6), 
										axis.text.y  = element_text( size = 6),
										axis.title.y = element_text( margin = margin( t = 0,r = 0, b = 0, l = 10)), 
										axis.title.x = element_text( margin = margin( t = 0,r = 0, b = 10, l = 0)), 
										legend.text = element_text( size = 6),
										legend.key.size = unit( 1, "points")
										 )
			
				)
	featurePlot <- plot_grid( plotlist = featurePlotList, ncol = 5, nrow = 3, label_size = 10) + theme_nothing() + theme( plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "inches"))
}else{
	featurePlotList <- FeaturePlot( seuratObj, 
				c( 	"tfap2e"	, "tfap2a"	, "her9"	, "sox9b"	, "foxg1b" 	,
					"snail2"	, "alx4b"	, "hmx1"	, "otx2"	, "sox10"	, 
					"impdh1b"	, "foxo1b"	, "tyr"		, "pax7b"	, "mc1r"	, 
					"id2a"		, "hmx4"	, "foxd3"	, "ednrba"	, "kita"	, 
					"mbpa"		, "phox2b"	, "tfec"	, "mitfa"	, "foxp4" 	, 
					"foxo1a"	, "hbp1"	, "ltk"		, "mycl1a"	, "pax7a"	, 
					"tyrp1b"	, "slc24a5"	, "oca2"	, "mlphb"	, "silva"	, 
					"myo5aa"	, "pnp4a"	, "ets1a"	, "fgfr3_v2"	, "pax3_v2"	, 
					"smad9"
),
	do.return = TRUE, nCol = 5, pt.size = 1,  min.cutoff = minCutoff, cols.use = c("blue", "red"), reduction.use = dimRed)   

	featurePlotList$clusters <- clusterPlot


	featurePlotList <- lapply( featurePlotList, function(cPlot) cPlot + theme( 	plot.margin = unit( c( 0, 0, 0, 0), "cm"),
										plot.title  = element_text( size = 12), 
										axis.title  = element_text( size = 8),
										axis.text.x = element_text(
											angle = 90,
											hjust = 0.1,
											vjust = 0.1,
											size =	6), 
										axis.text.y  = element_text( size = 6),
										axis.title.y = element_text( margin = margin( t = 0,r = 0, b = 0, l = 10)), 
										axis.title.x = element_text( margin = margin( t = 0,r = 0, b = 10, l = 0)), 
										legend.text = element_text( size = 6),
										legend.key.size = unit( 1, "points")
										 )
			
				)
	featurePlot <- plot_grid( plotlist = featurePlotList, ncol = 6, nrow = 7, label_size = 10) + theme_nothing() + theme( plot.margin = unit( c(0.5, 0.5, 0.5, 0.5), "inches"))
}

if( plotDPI == 600){
	
  featurePlotDir600dpi <- file.path( featurePlotDir, "600dpi")
  dir.create( featurePlotDir600dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = featurePlotDir600dpi, device = "png" , plot = featurePlot, width = pageWidth, height = pageHeight, units = "cm", dpi = 600, scale = 4)

}else if( plotDPI == 100){

  featurePlotDir100dpi <- file.path( featurePlotDir, "100dpi")
  dir.create( featurePlotDir100dpi, showWarnings = FALSE)

  ggsave( paste0( name, ".png"), path = featurePlotDir100dpi, device = "png" , plot = featurePlot, width = pageWidth, height = pageHeight, units = "cm", dpi = 100, scale = 4)

}else{ cat("Select plotDPI 600 or plotDPI 100")}

}

