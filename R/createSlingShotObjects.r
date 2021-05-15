createSlingShotObjects <- function( seurHD, seur2D, dimRedHD, dimRed2D, genesUseHD = NULL, useDims = 1:10, startClust = "eHMP", endClust = c("I", "M", "X"), distFun = cosineClusterDist){

#this snippet creates sling shot objects in HD and 2D trajectories based on seurat clustering obtain in HD needed for presentation. 
#This is done by substituting 2D lineages for HD lineages. HD can have dimensions of 2, the main point here is that one is trying
#use view plain saved in seur2D. Yet, not always this trick works. Sometimes the trajectories in the 2D become not very nice, 
#going into neighboring clusters
#it returns two slingshot objects as components of a list 

if(!require(slingshot)){
  install.packages("slingshot")
}

library("slingshot")	
library("matrixStats")


#check if the HD contains clusters for lineage starts and terminals
if(  length( setdiff( startClust, levels( seurHD@ident)) != 0) || ( length( setdiff( endClust, levels( seurHD@ident))) != 0 ))
	{ cat( "createSlingShotObject: Wrong cluster types HD"); return()} 


if ( !is.null( genesUseHD)) {cellsHD 	<- t( as.matrix(seurHD@data[ genesUseHD, ]))
}else{ cellsHD  <- GetCellEmbeddings( seurHD, reduction.type = dimRedHD)}

cells2D 	<- GetCellEmbeddings( seur2D, reduction.type = dimRed2D)

#Lineages connecting centers of the clusters are estimated in HD
slingShotHD 	<- slingshot( cellsHD, seurHD@ident, start.clus = startClust ,end.clus = endClust, reassign = TRUE, dist.fun = distFun )

#We don't have correct ident in 2D, so we ask for lineages in HD, this command is only to create the structure
slingShot2D 	<- slingshot( cells2D, seurHD@ident, start.clus = startClust, end.clus = endClust, dist.fun = distFun)

#now replace the lineage thee. The correct tree must be estimated in HD but we need curves in 2D for presentations
slingShot2D@lineages 	<- slingShotHD@lineages

#now estimate the curves
slingShot2D 		<- getCurves( slingShot2D, extend = "n", reassign = TRUE, stretch = 0, thresh = 0.005, shrink = 0.4)
slingShotHD		<- getCurves( slingShotHD, extend = "n", reassign = TRUE, stretch = 0, thresh = 0.005, shrink = 0.4)


#lineageIds <- which(unlist(lapply( slingLins2D@lineages, function(x)
#tail(x, 1) %in% c("M", "I"))))
#slingshot to here

slingShotObj		<- list()
slingShotObj$DimH	<- slingShotHD
slingShotObj$Dim2	<- slingShot2D


return( slingShotObj) 

}


