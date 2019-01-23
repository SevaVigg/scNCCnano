#Imputation for dropouts for logExpTable

imputeDropouts <- function( logExpsTable, randomSeed){

if(!require(DrImpute)){
  library(devtools)
  install_github('gongx030/DrImpute')		
  library("DrImpute")
}

set.seed(randomSeed)
logExpsTableImp <- DrImpute(logExpsTable)

return( logExpsTableImp)
}
