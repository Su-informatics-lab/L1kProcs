rm(list=ls())
ParseDataName <-  function(lstDataName){
  ParseDataName <- matrix(NA,nrow = 7,ncol = length(lstDataName))
  rownames(ParseDataName) <- c("Plate","Cell","Time","Rep",
                               "BeadSet","DetectMode","Well")
  colnames(ParseDataName) <- lstDataName
  for (i in 1:length(lstDataName)){
    tmpPList <- as.matrix(strsplit(lstDataName[i],"_")[[1]])
    if (8 == length(tmpPList)) {
      ParseDataName[,i] <- as.matrix(c(tmpPList[1:3],paste(tmpPList[4],tmpPList[5],sep="_"),
                                       tmpPList[6:8]))
    } else
    ParseDataName[,i] <- tmpPList
  }
  return (ParseDataName)
}

perPlot <- function(DataName,iStart,iEnd,lfc,result,GeneList,lambda,Overwrite=FALSE){
  for (i in iStart:iEnd){
    Analyte <- lfc[i == lfc[,"RID"],"EXP"]
    if(result[i,"n_peak"]==1 || (sum(is.na(result[i,]))!=0)) {
      cat("Analyte",i, "does not have enough beads support!\n")
      next
    } else {
      png(paste("../Plots/",DataName,"-Analyte",i,".png",sep=""))
      x1 <- 0.9*min(Analyte)
      x2 <- 1.1*max(Analyte)
      suppressWarnings(histdata <- hist(Analyte,freq=FALSE,breaks=20,plot=F))
      maxY <- 1.2*max(histdata$density)
      HistAnalyte <- hist(Analyte,breaks=20,freq=F,
                          xlim=range(x1,x2),ylim=range(0,maxY),
                          col="grey",
                          main="GMM Peak Calling",
                          xlab="Expression",ylab="Density")
      expA <- formatC(result[i,"c1"],format="f",digits=2)
      expB <- formatC(result[i,"c2"],format="f",digits=2) 
      sub.txt <- paste("n=",result[i,"n_tot"],"; exp=[",expA,",",expB,"]; supportNo=[",result[i,"n1"],",",result[i,"n2"],"]; err=",formatC(result[i,"err"],format="f",digits=2),sep="")
      mtext(sub.txt,side=3,line=0)
      inX <- seq(x1,x2,0.05)
      y1 <- lambda*dnorm(inX,mean=result[i,"c1"],sd=result[i,"sd1"])
      y2 <- (1-lambda)*dnorm(inX,mean=result[i,"c2"],sd=result[i,"sd2"])
      points(inX,y1,col="red",type="l",lwd=2)
      points(inX,y2,col="blue",type="l",lwd=2)
      pAx <- result[i,"c1"]
      pBx <- result[i,"c2"]
      pAy <-lambda*dnorm(pAx,mean=result[i,"c1"],sd=result[i,"sd1"]) 
      pBy <-(1-lambda)*dnorm(pBx,mean=result[i,"c2"],sd=result[i,"sd2"]) 
      points(pAx,pAy,type="p",col="red",pch=19)
      points(pBx,pBy,type="p",col="blue",pch=19)
      text(pAx,(pAy+0.05),paste("(A) ",GeneList[((i-10)*2-1)],": ",expA,sep=""),font=1.8,col="red")
      text(pBx,(pBy+0.05),paste("(B) ",GeneList[((i-10)*2)],":",expB,sep=""),font=1.8,col="blue")
      dev.off()
    }
  }
}
WPlot <- function(DataName,iStart,iEnd,Overwrite=FALSE){
  PList <- ParseDataName(DataName)
  DirName <- paste(PList["Plate",],
                   PList["Cell",],
                   PList["Time",],
                   PList["Rep",],
                   PList["BeadSet",],
                   PList["DetectMode",],
                   sep="_")
  ifLXB <- paste("../data/",DirName,"/",DataName,".lxb",sep="")
  ofCSV <- paste("../output/",DataName,".csv",sep="")
  ofGCT <- paste("../output/",DataName,".gct",sep="")
  GeneList <- read.delim(ofGCT,sep="\t",skip=2,stringsAsFactors=F)[,1]

  suppressWarnings(dat<-readFCS(ifLXB))
  data <- dat@exprs
  rawlfc <- cbind(data[,"RID"], log2(data[,"RP1"]))
  colnames(rawlfc) <- c("RID","EXP")
  lfc <- rawlfc[!is.infinite(rawlfc[,"EXP"]),]

  result <- read.delim(ofCSV,sep=",",stringsAsFactors=F)
  lambda <- 1.25/2
  perPlot(DataName,iStart,iEnd,lfc,result,GeneList,lambda)
}

args <- commandArgs(TRUE)
DataName <- args[1]
Mode <- as.integer(args[2])
cat("Loading library ...\n")
suppressMessages(library(prada))
if(Mode>10 & Mode <501){
  iStart <- Mode
  iEnd <- Mode
  No <- Mode
  }else{
    iStart <- 11
    iEnd <- 500
    No <- "all"
  }
cat("GMM Peak calling: ",DataName," Analyte ",No,"\n",sep="")
WPlot(DataName,iStart,iEnd)
cat("Finished.\n")
