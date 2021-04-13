makeFeaturePlots	<- function( seuratObj, minCutoff, dimRed, name = "featurePlots_"){

require( cowplot)
require( Seurat)
source("R/setClusterColors.r")

resDir		<- file.path(getwd(), "Res")

plotDir		<- file.path(resDir, "Plots")
dir.create(plotDir, showWarnings = FALSE)

featurePlotDir 	<- file.path( plotDir, "featurePlots")
dir.create( featurePlotDir, showWarnings = FALSE)

featurePlotFileName <- file.path( featurePlotDir, paste0(name, minCutoff, ".png"))

clusterPlot	    <- DimPlot( seuratObj, reduction.use = "umap", pt.size = 0.1, cols.use = setClusterColors( seuratObj), no.legend = TRUE, label.size = 0.5) +
			 	theme( plot.title = element_text( size = 12, face = "bold"), 
					axis.title.y = element_text( margin = margin( t = 0, r = 0, b = 0, l = 10)),  
					axis.title.x = element_text( margin = margin( t = 0, r = 0, b = 10, l = 0))) + 
				ggtitle( "clusters")  

#png( featurePlotFileName, height = 7016, width = 4960)
	featurePlotList <- FeaturePlot( seuratObj, 
				c( 	"tfap2e"	, "tfap2a"	, "her9"	, "sox9b"	, "foxg1b" 	,
					"snail2"	, "alx4b"	, "hmx1"	, "otx2"	, "sox10"	, 
					"impdh1b"	, "foxo1b"	, "tyr"		, "pax7b"	, "mc1r"	, 
					"id2a"		, "hmx4"	, "foxd3"	, "ednrba"	, "kita"	, 
					"mbpa"		, "phox2b"	, "tfec"	, "mitfa"	, "foxp4" 	, 
					"foxo1a"	, "hbp1"	, "ltk"		, "mycl1a"	, "pax7a"	, 
					"tyrp1b"	, "slc24a5"	, "oca2"	, "mlphb"	, "silva"	, 
					"myo5aa"	, "pnp4a"	, "ets1a"	, "fgfr3_v2"	, "pax3_v2"	, 
				#	"dpf3"		, 
					"smad9"
),
	do.return = TRUE, pt.size = 0.2,  min.cutoff = minCutoff, cols.use = c("blue", "red"), reduction.use = dimRed)   
#dev.off()

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
save_plot( filename = featurePlotFileName, plot = featurePlot, device = "png", dpi = 600, base_height = 11.7, base_width = 8.3) 

}

