qqPlotGeneExpressionBarChart <- funciton( geneExpMtx){

if(!require("tidyverse")){
  install.packages("tidyverse", dependencies = TRUE)
}

if(!require("ggplot2")){
  install.packages("ggplot2", dependencies = TRUE)
}



geneExpDf <- as.tibble( geneExpMtx)

}
