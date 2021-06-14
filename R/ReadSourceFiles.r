# This snippets reads source files for NanoString data
#
# called by ProcessData.r
#
# the files must be placed in the SourcePath which is the program parameter
#
# returns object CellTable, which contains three data frames:
#	Genes  - expression data
#	Cells  - cell meta data
#	Probes - NanoString probe data 
#
# NanoString NCC project 2017 - 2021, written by Vsevolod J. Makeev


ReadSourceFiles <- function(SourcePath){

# Genes containing gene descriptions. 

source("R/ReadNanoStringFile.r")
source("R/findDuplicated.r")
source("R/correctGeneNames.r")

rawPath		<- SourcePath
badPath		<- file.path( SourcePath, "BadFiles")

dir.create( badPath, showWarnings = FALSE)

filenames <- list.files(path = rawPath, pattern=".csv$", full.names = TRUE )

for( FileName in filenames){								#find the first correct file
	CellTable       <- ReadNanoStringFile( FileName)
	filenames	<- filenames[-1]	
	if(is.na( CellTable)){
		file.rename( from = FileName,  to = file.path( badPath, basename(FileName))); 
		next }else{ break}
}


for( FileName in filenames){
	NewTable	<- ReadNanoStringFile( FileName)
	filenames	<- filenames[-1]
	if(is.na(NewTable)){file.rename(from = FileName,  to = file.path( badPath, basename(FileName))); next }

	
	if( !all(unlist(lapply( rownames(NewTable$Cells), function(Var) {
			if( all(!is.element(NewTable$Cells[Var,], CellTable$Cells[Var,])) 
			    & identical(CellTable$Genes, NewTable$Genes)){
			return( TRUE)}else{return(FALSE)}
								})))){
			CellTable$Genes <- cbind(CellTable$Genes, NewTable$Genes)
			CellTable$Cells <- cbind(CellTable$Cells, NewTable$Cells)
	    }else{
		cat(FileName, " already present", "\n")
	}
}

DescNames		<- rownames(CellTable$Cells)

#some gene names are not those from ZFIN. correctGeneNames.r verifies genenames, and make corrections

CellTable$Probes[ ,"Gene Name"]	<- correctGeneNames( CellTable$Probes[ , "Gene Name"])

Genes 			<- as.data.frame( lapply(CellTable$Genes, unlist), stringsAsFactors = FALSE)
rownames(CellTable$Genes)		<- CellTable$Probes[, "Gene Name"] 
CellTable$Cells 	<- as.data.frame( lapply(CellTable$Cells, unlist), stringsAsFactors = FALSE)
CellTable$Cells		<- as.data.frame( lapply(CellTable$Cells, as.character), stringsAsFactors = FALSE)
CellTable$Probes	<- as.data.frame( lapply(CellTable$Probes, as.character), stringsAsFactor = FALSE)
rownames(CellTable$Cells)	<- DescNames

ans <- CellTable

} #main



