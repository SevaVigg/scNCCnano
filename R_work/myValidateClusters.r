myValidateClusters <- function (object, pc.use = NULL, top.genes = 30, min.connectivity = 0.01, 
    acc.cutoff = 0.9, verbose = TRUE, iterations = 5) 
{
    PackageCheck("caret")
    if (length(x = object@snn) > 1) {
        SNN.use <- object@snn
    }
    else {
        stop("SNN matrix required. Please run BuildSNN() to save the SNN matrix in\n         the object slot")
    }
    if (is.null(pc.use)) {
        stop("pc.use not set. Please choose PCs.")
    }
    num.clusters.orig <- length(x = unique(x = object@ident))
    still_merging <- TRUE
    if (verbose) {
        connectivity <- CalcConnectivity(object = object)
        end <- length(x = connectivity[connectivity > min.connectivity])
        progress <- end
        status <- 0
    }
    while (still_merging) {
        connectivity <- CalcConnectivity(object = object)
        merge.done <- FALSE
        while (!merge.done) {
            m <- max(connectivity, na.rm = TRUE)
            mi <- which(x = connectivity == m, arr.ind = TRUE)
            c1 <- rownames(x = connectivity)[mi[, 1]]
            c2 <- rownames(x = connectivity)[mi[, 2]]
            if (m > min.connectivity) {
                acc <- RunClassifier(object = object, group1 = c1, 
                  group2 = c2, pcs = pc.use, num.genes = top.genes)
                if (acc < acc.cutoff) {
                  object <- SetIdent(object = object, cells.use = WhichCells(object = object, 
                    ident = c1), ident.use = c2)
                  if (verbose) {
                    progress <- length(x = connectivity[connectivity > 
                      min.connectivity])
                    message(paste0(sprintf("%3.0f", (1 - progress/end) * 
                      100), "% complete --- merge clusters ", 
                      c1, " and ", c2, ", classification accuracy of ", 
                      sprintf("%1.4f", acc)))
                  }
                  merge.done <- TRUE
                }
                else {
                  if (verbose & status == iterations) {
                    message(paste0(sprintf("%3.0f", (1 - progress/end) * 
                      100), "% complete --- Last 5 cluster comparisons failed to merge, ", 
                      "still checking possible merges ..."))
                    status <- 0
                  }
                  status <- status + 1
                  connectivity[c1, c2] <- 0
                  connectivity[c2, c1] <- 0
                }
            }
            else {
                still_merging <- FALSE
                break
            }
        }
    }
    if (verbose) {
        message(paste0("100% complete --- started with ", num.clusters.orig, 
            " clusters, ", length(x = unique(x = object@ident)), 
            " clusters remaining"))
    }
    return(object)
}

