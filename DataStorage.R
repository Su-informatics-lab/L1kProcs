DataStorage <- function(Software){
# Processing all KD data
setwd(Software)
rm(list=ls())
lstAllFileNames <- unlist(strsplit(list.files(path="../data",pattern="*.lxb",recursive=TRUE),split="/"))
lstDataNames <- lstAllFileNames[grep(pattern="*.lxb",lstAllFileNames)]
lstNames <- sub(".lxb","",lstDataNames)
dir.create("Names")
tmp <- lapply(lstNames,function(x) file.create(paste("Names/",x,sep="")))

save(lstNames,file = "lstNames.RData")
lstPlates <- unique(lstAllFileNames[grep(pattern="HI..LO$",lstAllFileNames)])
save(lstPlates,file = "lstPlates.RData")
return(c(length(lstPlates),length(lstNames)))
}
