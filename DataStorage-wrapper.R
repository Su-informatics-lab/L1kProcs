args <- commandArgs(TRUE)
script <- paste(args,"/DataStorage.R",sep="")
source(script)
lstLen <- DataStorage(args)
cat(paste(as.character(lstLen[1]),"\t",
          as.character(lstLen[2])))
