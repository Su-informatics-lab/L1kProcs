args <- commandArgs(TRUE)
iStart <- as.integer(args[1])
iEnd <- as.integer(args[2])
Path <- args[3]
Overwrite <- as.logical(as.integer(args[4]))
script <- paste(Path,"/Preprocess.R",sep="")
source(script)
Preprocess(iStart,iEnd,Path,Overwrite)

