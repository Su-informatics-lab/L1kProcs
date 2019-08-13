args <- commandArgs(TRUE)
iStart <- as.integer(args[1])
iEnd <- as.integer(args[2])
Path <- args[3]
nThread <- as.integer(args[4])
Overwrite <- as.logical(as.integer(args[5]))
script <- paste(Path,"/PipeLine.R",sep="")
source(script)
PipeLine(iStart,iEnd,Path,nThread,Overwrite)

