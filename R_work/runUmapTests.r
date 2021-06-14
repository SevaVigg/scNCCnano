source("R/makeUMAPClusteringPlots.r")
for (umapDim in c(10, 15, 17)){
  for (clResolution in c(0.1, 0.2, 1, 1.5, 2)){

	makeUMAPClusteringPlots(ipmc, umapDim, clResolution)

}}
