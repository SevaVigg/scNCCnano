#Developed by Leonid Uroshlev
#This function allows getting transfer data for UMAP

if(!require(reticulate)){install.packages("reticulate")}
library(reticulate)

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

RunUMAP2<-function (object, cells.use = NULL, dims.use = 1:5, reduction.use = "pca", 
    genes.use = NULL, assay.use = "RNA", max.dim = 2L, reduction.name = "umap", 
    reduction.key = "UMAP", n_neighbors = 30L, min_dist = 0.3, 
    metric = "correlation", seed.use = 42, ...) 
{
    if (!py_module_available(module = "umap")) {
        stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
    }
    if (!is.null(x = seed.use)) {
        set.seed(seed = seed.use)
        py_set_seed(seed = seed.use)
    }
    cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
    if (is.null(x = genes.use)) {
        dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
            slot = "key")
        dim.codes <- paste0(dim.code, dims.use)
        data.use <- GetDimReduction(object = object, reduction.type = reduction.use, 
            slot = "cell.embeddings")
        data.use <- data.use[cells.use, dim.codes, drop = FALSE]
    }
    else {
        data.use <- GetAssayData(object = object, assay.type = assay.use, 
            slot = "scale.data")
        genes.use <- intersect(x = genes.use, y = rownames(x = data.use))
        if (!length(x = genes.use)) {
            stop("No genes found in the scale.data slot of assay ", 
                assay.use)
        }
        data.use <- data.use[genes.use, cells.use, drop = FALSE]
        data.use <- t(x = data.use)
    }
    parameters.to.store <- as.list(x = environment(), all = TRUE)[names(formals("RunUMAP"))]
    object <- SetCalcParams(object = object, calculation = "RunUMAP", 
        ... = parameters.to.store)
    umap_import <- import(module = "umap", delay_load = TRUE)
    umap <- umap_import$UMAP(n_neighbors = as.integer(x = n_neighbors), 
        n_components = as.integer(x = max.dim), metric = metric, 
        min_dist = min_dist)
    umap_output <- umap$fit_transform(as.matrix(x = data.use))
    colnames(x = umap_output) <- paste0(reduction.key, 1:ncol(x = umap_output))
    rownames(x = umap_output) <- cells.use
    object <- SetDimReduction(object = object, reduction.type = reduction.name, 
        slot = "cell.embeddings", new.data = as.matrix(x = umap_output))
    object <- SetDimReduction(object = object, reduction.type = reduction.name, 
        slot = "key", new.data = reduction.key)
    return(c(object=object, umap=umap))
}

# Вызывать надо вот так umap_output <- umap$fit_transform(as.matrix(x = А_вот_сюда_новые_данные))
