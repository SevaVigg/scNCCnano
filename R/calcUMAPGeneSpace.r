calcUmapGeneSpace <- function( seurAll, seurWT = NULL, experimentType = "allCells",
			UMAPRandSeed = 42L, Dim = 2, minDist = 0.65, myNeighbors = 20L, myMetric = "cosine",  
			assay.use = "RNA", reduction.key = "UMAP", reduction.name = "umap", mySpread = 1){

#experimentType = c("allCells", "allCondWT", "WToutAll")
# 'allCells' performs standard UMAP of seurAll
# 'WToutAll' takes out sox- cells and performs UMAP only WT
# 'allCondWT' uses WT to fit UMAP and performs transform for all cells 
# we cannot subset to WT because WT depends on imputation on the WT cells only
#UMAP parameters

umapRes	<- list( All = list(), WT = list())

if(!require(reticulate)){install.packages("reticulate")}
library(reticulate)

#np <- import(module = "numpy")
#np$np.random.seed(42)

SetIfNull <- function(x, default) {
  if (is.null(x)) {
    return(default)
  } else {
    return(x)
  }
}

SetCalcParams <- function(object, calculation, time = TRUE, ...) {
  object@calc.params[calculation] <- list(...)
  object@calc.params[[calculation]]$object <- NULL
  object@calc.params[[calculation]]$object2 <- NULL
  if(time) {
    object@calc.params[[calculation]]$time <- Sys.time()
  }
  return(object)
}

cellsWT <- setdiff( colnames( seurAll@data), grep( "sox", colnames(seurAll@data), value = TRUE))
allData	<- GetAssayData( object = seurAll, assay.type = assay.use, 
            		slot = "scale.data")
transformCells	<- colnames( x = seurAll@data)
allData		<- t(allData)

#now calculate UMAP for init data
parameters.to.store <- as.list(x = environment(), all = TRUE)[names(formals("RunUMAP"))]

#fetch UMAP
    umap_import <- import(module = "umap", delay_load = TRUE)
    umap <- umap_import$UMAP(n_neighbors = as.integer(x = myNeighbors), 
        n_components = as.integer(x = Dim), metric = myMetric, spread = as.numeric( x = mySpread),  
        min_dist = minDist, transform_seed = as.integer( x = UMAPRandSeed))

if( experimentType == "allCells") 
{ 
		transformAll 	<- umap$fit_transform( as.matrix( x = allData))
    		allResObj 	<- seurAll
    		colnames(x = transformAll) <- paste0(reduction.key, 1:ncol(x = transformAll))
    		rownames(x = transformAll) <- transformCells
    		allResObj 	<- SetCalcParams(object = allResObj, calculation = "RunUMAP", 
        				... = parameters.to.store)
    		allResObj 	<- SetDimReduction(object = allResObj, reduction.type = reduction.name, 
        				slot = "cell.embeddings", new.data = as.matrix(x = transformAll))
    		allResObj 	<- SetDimReduction(object = allResObj, reduction.type = reduction.name, 
        				slot = "key", new.data = reduction.key)
		umapRes$All <- allResObj
		return( umapRes)
	}else{
	if( experimentType == "WToutAll"){
		seurWT <- SubsetData( seurAll, cells.use = cellsWT)
		wtData <- allData[ cellsWT, ]
    		transformWT	<- umap$fit_transform( as.matrix( x = wtData))
    		wtResObj <- seurWT
    		colnames(x = transformWT) <- paste0(reduction.key, 1:ncol(x = transformWT))
    		rownames(x = transformWT) <- cellsWT
    		wtResObj <- SetCalcParams(object = wtResObj, calculation = "RunUMAP", 
        		... = parameters.to.store)
    		wtResObj <- SetDimReduction(object = wtResObj, reduction.type = reduction.name, 
        		slot = "cell.embeddings", new.data = as.matrix(x = transformWT))
    		wtResObj <- SetDimReduction(object = wtResObj, reduction.type = reduction.name, 
        		slot = "key", new.data = reduction.key)
		umapRes$WT <- wtResObj
	return( umapRes)
}else{

#now we have "allCondWT"

if( experimentType == "allCondWT"){
#if there is a separate WT object use it to take train data

if( !is.null( seurWT)){ wtData <- GetAssayData(seurWT, slot = "scale.data"); wtData <- t(wtData)}
#otherwise take data by subsetting the allData
		       else{ wtData <- allData[ cellsWT, ]}

#we fit on the wtData and transform all data
    umap 		<- umap$fit( as.matrix( x = wtData))

#transform wt data
    transformWT		<- umap$transform( as.matrix( x = wtData))
    
    wtResObj <- seurWT
    colnames(x = transformWT) <- paste0(reduction.key, 1:ncol(x = transformWT))
    rownames(x = transformWT) <- cellsWT
    wtResObj <- SetCalcParams(object = wtResObj, calculation = "RunUMAP", 
        ... = parameters.to.store)
    wtResObj <- SetDimReduction(object = wtResObj, reduction.type = reduction.name, 
        slot = "cell.embeddings", new.data = as.matrix(x = transformWT))
    wtResObj <- SetDimReduction(object = wtResObj, reduction.type = reduction.name, 
        slot = "key", new.data = reduction.key)
 
#transform All data
    transformAll 	<- umap$transform( as.matrix( x = allData))
    allResObj <- seurAll
    colnames(x = transformAll) <- paste0(reduction.key, 1:ncol(x = transformAll))
    rownames(x = transformAll) <- transformCells
    allResObj <- SetCalcParams(object = allResObj, calculation = "RunUMAP", 
        ... = parameters.to.store)
    allResObj <- SetDimReduction(object = allResObj, reduction.type = reduction.name, 
        slot = "cell.embeddings", new.data = as.matrix(x = transformAll))
    allResObj <- SetDimReduction(object = allResObj, reduction.type = reduction.name, 
        slot = "key", new.data = reduction.key)

#prepare return object

    umapRes$All <- allResObj; umapRes$WT <- wtResObj
    return( umapRes)
}else{ cat("Unknown experiment type")}
}	#allCondWT
}	
}
