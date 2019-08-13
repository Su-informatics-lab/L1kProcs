########################################
#              Utilities               #
########################################

## Author: Jing Su
## Date: May 2012
## License: GPL-2

########################################
#            FindControls              #
########################################
FindControls <- function(lstDataNames,PlateMap,KeyWord,flag=0) {
  FindControls <- data.frame(Sample=character(0),Control=character(0),
                             stringsAsFactors = FALSE)

  for (DataName in lstDataNames){
    PList <- ParseDataName(DataName)
    lstControlWells <- c()
    if(flag==1){
      tmp <- paste(PList["Plate",],PList["Cell",],PList["Time",],
                   PList["Rep",],PList["BeadSet",],PList["DetectMode",],
                   sep="_")
      lstControlWells <- (PlateMap[PlateMap$PerturbagenPlate==tmp
                                 & PlateMap$PerturbagenID==KeyWord,])$Well
    }else if(flag ==0){
      lstControlWells <- (PlateMap[PlateMap$PerturbagenPlate==PList["Plate",]
                                   & PlateMap$PerturbagenDesc==KeyWord,])$Well}
    for (ControlWell in lstControlWells){
      FindControls <- rbind(FindControls,
                            data.frame("Sample"=DataName,
                                       "Control"=paste(PList["Plate",],
                                         PList["Cell",],
                                         PList["Time",],
                                         PList["Rep",],
                                         PList["BeadSet",],
                                         PList["DetectMode",],
                                         ControlWell,
                                         sep="_"),
                                       row.names=c(paste(DataName,
                                         ControlWell,
                                         sep="-")
                                         ),
                                       stringsAsFactors = FALSE
                                       )
                            )
    }
  }
   return (FindControls)
}

########################################
#          FindDrugControls            #
########################################
FindDrugControls <- function(lstDataNames,PlateMap) 
  return(FindControls(lstDataNames,PlateMap,"DMSO"))

########################################
#            FindKDControls            #
########################################
FindKDControls <- function(lstDataNames,PlateMap)
  return(FindControls(lstDataNames,PlateMap,"EMPTY_VECTOR"))

########################################
#            FindOEControls            #
########################################
FindOEControls <- function(lstDataNames,PlateMap) 
  return(FindControls(lstDataNames,PlateMap,"UnTrt"))

########################################
#          FindDrugControlsbac            #
########################################
FindDrugControlsbac <- function(lstDataNames,PlateMap){
  return(FindControls(lstDataNames,PlateMap,"DMSO",flag=1))
}

########################################
#            FindKDControlsbac            #
########################################
FindKDControlsbac <- function(lstDataNames,PlateMap)
  return(FindControls(lstDataNames,PlateMap,"EMPTY_VECTOR",flag=1))

########################################
#               LXB2Stats              #
########################################
GMM <- function(inX, inP) {
  lambda <- 1.25/2
  return(lambda*pnorm(inX,mean=inP[1],sd=inP[2])
         + (1-lambda)*pnorm(inX,mean=inP[3],sd=inP[4]))
}

FGMM <- function(inP,inX,inY) sum((inY-GMM(inX,inP))^2)

GMM2 <- function(inX,inP){
  return(inP[5]*dnorm(inX,mean=inP[1],sd=inP[2])
         + (1-inP[5])*dnorm(inX,mean=inP[3],sd=inP[4]))
}

GGMM2 <- function(inX,inP){
  return(inP[5]*pnorm(inX,mean=inP[1],sd=inP[2])
         + (1-inP[5])*pnorm(inX,mean=inP[3],sd=inP[4]))
}
FGMM2 <- function(inP,inX,inY) {
  ##  cat(length(inP),length(inX),length(inY),"\n")
  sum((inY-GMM2(inX,inP))^2)
}


LXB2Stats <- function(DataName,Overwrite=FALSE){
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
  cat ("\t\tFrom", ifLXB,"\n\t\tto",ofCSV,"\n")

  if (file.exists(ofCSV)) {
    if (Overwrite)
      file.remove(ofCSV)
    else {
      cat(ofCSV," exists. Skip. \n")
      return()
    }
  }

  dat<-readFCS(ifLXB)
  data <- dat@exprs
  rawlfc <- cbind(data[,"RID"], log2(data[,"RP1"]))
  colnames(rawlfc) <- c("RID","EXP")

  lfc <- rawlfc[!is.infinite(rawlfc[,"EXP"]),]

  results <- mat.or.vec(500,11)
  colnames(results) <- c("AnalyteN","n1","c1","sd1","n2","c2","sd2","n_peak","n_tot","err","conv")
  results[,"AnalyteN"] <- 1:500

  lambda <- 1.25/2
  for (i in 1:500){
    Analyte <- lfc[i == lfc[,"RID"],"EXP"]
    results[i,"n_tot"] <- length(Analyte)
    if ( i <= 10) {
      results[i,"n_peak"] <- 1
      results[i,"n1"] <- length(Analyte)
      results[i,"c1"] <- median(Analyte)
      results[i,"sd1"] <- sd(Analyte)
      results[i,"n2"] <- NA
      results[i,"c2"] <- NA
      results[i,"sd2"] <- NA
    } else if (1==length(lfc[i == lfc[,"RID"],"EXP"]))
      results[i,"n_peak"] <- 1
    else if (0==length(lfc[i == lfc[,"RID"],"EXP"])) {
      results[i,"n1"] <- NA
      results[i,"c1"] <- NA
      results[i,"n2"] <- NA
      results[i,"c2"] <- NA
      next
    } else if (2==length(lfc[i == lfc[,"RID"],"EXP"])) {
      results[i,"c1"] <- max(lfc[i == lfc[,"RID"],"EXP"])
      results[i,"n1"] <- 1
      results[i,"sd1"] <- 0
      results[i,"c2"] <- min(lfc[i == lfc[,"RID"],"EXP"])
      results[i,"n2"] <- 1
      results[i,"sd2"] <- 0
      next
    } else if (1==length(unique(Analyte))) {
      results[i,"n_peak"] <- 2
      results[i,"n1"] <- round(results[i,"n_tot"]*lambda)
      results[i,"n2"] <- results[i,"n_tot"]-results[i,"n1"]
      results[i,"c1"] <- median(Analyte)
      results[i,"sd1"] <- sd(Analyte)
      results[i,"c2"] <- median(Analyte)
      results[i,"sd2"] <- sd(Analyte)
      results[i,"err"] <- NA
      results[i,"conv"] <- NA
    } else {
      results[i,"n_peak"] <- 2
      HistAnalyte <- hist(Analyte,breaks=500,plot=FALSE)
      CDFAnalyte <- cumsum(HistAnalyte$counts/sum(HistAnalyte$counts))
      ## fit <- optim(c(2*median(Analyte)/3,sd(Analyte),4*median(Analyte)/3,sd(Analyte)),
      ##              FGMM,NULL,HistAnalyte$mids,CDFAnalyte,control=list(maxit=10000))
      km <- kmeans(Analyte,centers=2)

      ## fit <- optim(c(ifelse(km$size[1]>=km$size[2], km$centers[1], km$centers[2]),
      ##                sd(Analyte)/2,
      ##                ifelse(km$size[1]>=km$size[2], km$centers[2], km$centers[1]),
      ##                sd(Analyte)/2),
      ##              FGMM,NULL,HistAnalyte$mids,CDFAnalyte,control=list(maxit=10000))
      fit <- optim(c(ifelse(km$size[1]>=km$size[2], km$centers[1], km$centers[2]),
                     sd(Analyte)/2,
                     ifelse(km$size[1]>=km$size[2], km$centers[2], km$centers[1]),
                     sd(Analyte)/2,
                     1.25/2),
                   FGMM,NULL,
                   HistAnalyte$mids,CDFAnalyte,
                   control=list(maxit=10000))
      results[i,"n1"] <- round(results[i,"n_tot"]*lambda)
      results[i,"n2"] <- results[i,"n_tot"]-results[i,"n1"]
      results[i,"c1"] <- fit$par[1]
      results[i,"sd1"] <- fit$par[2]
      results[i,"c2"] <- fit$par[3]
      results[i,"sd2"] <- fit$par[4]
      results[i,"err"] <- fit$value
      results[i,"conv"] <- fit$conv
      fit <- optim(c(ifelse(km$size[1]<km$size[2], km$centers[1], km$centers[2]),
                     sd(Analyte)/2,
                     ifelse(km$size[1]<km$size[2], km$centers[2], km$centers[1]),
                     sd(Analyte)/2,
                     1.25/2),
                   FGMM,NULL,
                   HistAnalyte$mids,CDFAnalyte,
                   control=list(maxit=10000))
      if (results[i,"err"] > fit$value) {
        results[i,"c1"] <- fit$par[1]
        results[i,"sd1"] <- fit$par[2]
        results[i,"c2"] <- fit$par[3]
        results[i,"sd2"] <- fit$par[4]
        results[i,"err"] <- fit$value
        results[i,"conv"] <- fit$conv
      }
    }
  }
#  save(results,dat,rawlfc,lfc,file=outdatafile,compress=TRUE)
  write.csv(results,file=ofCSV)
}
LXB2StatsBAC <- function(DataName){
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
  cat ("From", ifLXB,"\n\tto",ofCSV,"\n")

  dat<-readFCS(ifLXB)
  data <- dat@exprs
  rawlfc <- cbind(data[,"RID"], log2(data[,"RP1"]))
  colnames(rawlfc) <- c("RID","EXP")

  lfc <- rawlfc[!is.infinite(rawlfc[,"EXP"]),]

  results <- mat.or.vec(500,7)
  colnames(results) <- c("AnalyteN","n1","c1","n2","c2", "n_peak","n_tot")
  results[,"AnalyteN"] <- 1:500
  
  for (i in 1:500){
#    cat("Processing", i, "of 500... ")
    if ( i <= 10) results[i,"n_peak"] <- 1 else results[i,"n_peak"] <- 2
    if (1==length(lfc[i == lfc[,"RID"],"EXP"])) results[i,"n_peak"] <- 1
    if (0==length(lfc[i == lfc[,"RID"],"EXP"])) {
      results[i,"n1"] <- 0
      results[i,"c1"] <- 0
      results[i,"n2"] <- NA
      results[i,"c2"] <- NA
      next
    }
    if (2==length(lfc[i == lfc[,"RID"],"EXP"])) {
      results[i,"n1"] <- max(lfc[i == lfc[,"RID"],"EXP"])
      results[i,"c1"] <- 1
      results[i,"n2"] <- min(lfc[i == lfc[,"RID"],"EXP"])
      results[i,"c2"] <- 1
      next
    }
    fit <- kmeans(lfc[i == lfc[,"RID"],"EXP"],results[i,"n_peak"], nstart=10)
    if (length(fit$centers) == 1){
      results[i,"n1"] <- fit$size[1]
      results[i,"c1"] <- fit$centers[1]
    } else {
      results[i,c("n1","c1","n2","c2")] <-
        if(fit$centers[1] > fit$centers[2]){
          c(fit$size[1],fit$centers[1],fit$size[2],fit$centers[2])
        } else c(fit$size[2],fit$centers[2],fit$size[1],fit$centers[1])
    }
#    cat("\tDone! \n")
  }
  results[,"n_tot"] <- results[,"n1"] + results[,"n2"];
  results[1:10,"n_tot"] <- results[1:10,"n1"];
#  save(results,dat,rawlfc,lfc,file=outdatafile,compress=TRUE)
  write.csv(results,file=ofCSV)
  cat("\tDone! \n")
}
########################################
#               Stats2gct              #
########################################
Stats2gct <- function(statsDataName, Overwrite=FALSE){
  ifBead2Gene <- "../Annotation/bead_gene_map.csv"
  ifCSV <- paste("../output/",statsDataName,".csv",sep="")
  ofGCT <- paste("../output/",statsDataName,".gct",sep="")

  if (file.exists(ofGCT)) {
    if (Overwrite)
      file.remove(ofGCT)
    else {
      cat(ofGCT," exists. Skip. \n")
      return()
    }
  }

  ifPlateMap <- "../Annotation/plate_maps.csv"
  PlateMap <- read.csv(ifPlateMap, header = TRUE, stringsAsFactors = FALSE)
  ifU133A <- "../Annotation/HG_U133A.chip"
  U133A <- read.delim(ifU133A, header = TRUE, stringsAsFactors = FALSE)
  rownames(U133A) <- U133A$pr_id
  ifBead2Gene <- "../Annotation/bead_gene_map.csv"
  Bead2Gene <- read.csv(ifBead2Gene, header = TRUE, stringsAsFactors = FALSE)
  DPeaks <- read.csv(ifCSV, header = TRUE, stringsAsFactors = FALSE)
  PList <- as.matrix(strsplit(statsDataName,"_")[[1]])
  tmpPList <- as.matrix(strsplit(statsDataName,"_")[[1]])
  if (8 == length(tmpPList)) {
    PList <- as.matrix(c(tmpPList[1:3],paste(tmpPList[4],tmpPList[5],sep="_"),tmpPList[6:8]))
  } else
  PList <- tmpPList
  rownames(PList) <- c("Plate","Cell","Time","Rep",
                       "BeadSet","DetectMode","Well")
  GeneAnno <- Bead2Gene[grep(PList["BeadSet",], Bead2Gene$BeadSetBatch,fixed=TRUE),]
  gctFile <- file(ofGCT,open = "a")
  cat(paste("#1.2","\n",sep=""),
      file = gctFile,append=TRUE)
  cat(paste("1000","\t","1","\n",sep=""),
      file = gctFile,append=TRUE)
  cat(paste("Gene","\t",
            "Name","\t",
            "Description","\t",
            statsDataName,
            "\n",sep=""),
      file = gctFile,append=TRUE)
  for (i in 1:nrow(DPeaks)) {
    # Handling the HI probes
    if ("DUO52HI53LO"==PList["DetectMode",]) {
      tmpGeneID <- GeneAnno[GeneAnno$LMNXAnalyteID ==
                            paste("Analyte ",as.character(i),sep="") &
                            GeneAnno$BeadSet == "DP52",]
      tmpProbe <- "DP52"
    } else if ("DUO53HI52LO"==PList["DetectMode",]) {
      tmpGeneID <- GeneAnno[GeneAnno$LMNXAnalyteID ==
                            paste("Analyte ",as.character(i),sep="") &
                            GeneAnno$BeadSet == "DP53",]
      tmpProbe <- "DP53"
    } else if ("DUO45HI44LO"==PList["DetectMode",]) {
      tmpGeneID <- GeneAnno[GeneAnno$LMNXAnalyteID ==
                            paste("Analyte ",as.character(i),sep="") &
                            GeneAnno$BeadSet == "DP45",]
      tmpProbe <- "DP45"
    } else if ("DUO44HI45LO"==PList["DetectMode",]) {
      tmpGeneID <- GeneAnno[GeneAnno$LMNXAnalyteID ==
                            paste("Analyte ",as.character(i),sep="") &
                            GeneAnno$BeadSet == "DP44",]
      tmpProbe <- "DP44"
    }
    tmpU133AID <- (U133A$pr_id[U133A$pr_gene_id ==
                               tmpGeneID[1,"EntrezGeneID"]])[1]


    if (!is.na(tmpU133AID)) {
      if (DPeaks[i,"c1"] <=0 | DPeaks[i,"c1"] > 20 | is.na(DPeaks[i, "c1"])|DPeaks[i,"n_tot"]<20) {
        tmpDPeaks <- "NA"
      } else if (1 == DPeaks[i,"c1"] & 1 == DPeaks[i,"c2"] & 2 == DPeaks[i,"n_tot"]){
        tmpDPeaks <- DPeaks[i,"n1"] # Due to a fixed bug in LXB2Stats()
      } else {
        tmpDPeaks <- DPeaks[i,"c1"]
      }
      
      cat(paste(U133A[tmpU133AID,"pr_gene_symbol"],"\t",
                tmpU133AID,"\t",
                paste("Analyte ",as.character(i),", ", tmpProbe, ", ",
                      U133A[tmpU133AID,"pr_gene_title"],sep=""),"\t",
                tmpDPeaks,
                "\n",sep=""),
          file = gctFile,append=TRUE)
    }
    # Handling the LO probes
    if ("DUO52HI53LO"==PList["DetectMode",]) {
      tmpGeneID <- GeneAnno[GeneAnno$LMNXAnalyteID ==
                            paste("Analyte ",as.character(i),sep="") &
                            GeneAnno$BeadSet == "DP53",]
      tmpProbe <- "DP53"
    } else if ("DUO53HI52LO"==PList["DetectMode",]) {
      tmpGeneID <- GeneAnno[GeneAnno$LMNXAnalyteID ==
                            paste("Analyte ",as.character(i),sep="") &
                            GeneAnno$BeadSet == "DP52",]
      tmpProbe <- "DP52"
    } else if ("DUO45HI44LO"==PList["DetectMode",]) {
      tmpGeneID <- GeneAnno[GeneAnno$LMNXAnalyteID ==
                            paste("Analyte ",as.character(i),sep="") &
                            GeneAnno$BeadSet == "DP44",]
      tmpProbe <- "DP44"
    } else if ("DUO44HI45LO"==PList["DetectMode",]) {
      tmpGeneID <- GeneAnno[GeneAnno$LMNXAnalyteID ==
                            paste("Analyte ",as.character(i),sep="") &
                            GeneAnno$BeadSet == "DP45",]
      tmpProbe <- "DP45"
    }
    tmpU133AID <- (U133A$pr_id[U133A$pr_gene_id ==
                               tmpGeneID[1,"EntrezGeneID"]])[1]
    if (!is.na(tmpU133AID)) {
      if (DPeaks[i,"c2"] <=0 | DPeaks[i,"c2"] > 20 |is.na(DPeaks[i, "c2"])|DPeaks[i,"n_tot"]<20) {
        tmpDPeaks <- "NA"
      } else if (1 == DPeaks[i,"c1"] & 1 == DPeaks[i,"c2"] & 2 == DPeaks[i,"n_tot"]){
        tmpDPeaks <- DPeaks[i,"n2"] # Due to a fixed bug in LXB2Stats()
      } else {
        tmpDPeaks <- DPeaks[i,"c2"]
      }
      cat(paste(U133A[tmpU133AID,"pr_gene_symbol"],"\t",
                tmpU133AID,"\t",
                paste("Analyte ",as.character(i),", ", tmpProbe, ", ",
                      U133A[tmpU133AID,"pr_gene_title"],sep=""),"\t",
                tmpDPeaks,
                "\n",sep=""),
          file = gctFile,append=TRUE)
    }
  }
  close(gctFile)
}

########################################
#             QTargetGenerate          #
########################################
QTargetGenerate <- function(lstTargetNames){
  ProbeNumSet <- mat.or.vec(nr=980,nc=0)
  for(TargetName in lstTargetNames){
    ifGCT <- paste("../output/",TargetName,".gct",sep="")
    tmpGCT <- read.delim(ifGCT,sep="\t",
                         header = TRUE,skip = 2,
                         stringsAsFactors = FALSE)
    if (nrow(tmpGCT)==978)
      tmpDataSet <- c(tmpGCT[[TargetName]],"NA","NA")
    else 
      tmpDataSet <- tmpGCT[[TargetName]]
    ProbeNumSet <- cbind(ProbeNumSet, tmpDataSet)
  }
  numProbeNumSet <- apply(ProbeNumSet,2,as.numeric)
  colnames(numProbeNumSet) <- lstTargetNames
  tmpNA <- is.na(numProbeNumSet)
  colBad <- colnames(tmpNA)[apply(tmpNA,2,sum)>50]
  colGood <- setdiff(colnames(numProbeNumSet),colBad)
  QTarget <- normalize.quantiles.determine.target(numProbeNumSet[,colGood])
  return(QTarget)
}
########################################
#            QNormgctfile              #
########################################
QNormgctfile <- function(lstDataNames,nthread,QTarget){
  lsfProbeSet <- list()
  lstProbe <- list()
  for(RawDataName in lstDataNames){
    ifGCT <- paste("../output/",RawDataName,".gct",sep="")
    tmpGCT <- read.delim(ifGCT,sep="\t",
                         header = TRUE,skip = 2,
                         stringsAsFactors = FALSE)
    ProbeName <- ParseDataName(RawDataName)["ShortModeID",]
    if(!is.element(ProbeName,names(lsfProbeSet))){
      lsfProbeSet[[ProbeName]] <- tmpGCT
      lstProbe[ProbeName] <- ProbeName
    }else{
      tmpDataSet<- data.frame(
                     tmpGCT[match(lsfProbeSet[[ProbeName]]$Gene,
                                  tmpGCT$Gene),RawDataName],
                     stringsAsFactors = FALSE)
      colnames(tmpDataSet) <- RawDataName
      lsfProbeSet[[ProbeName]] <- cbind(lsfProbeSet[[ProbeName]],tmpDataSet)
    }
  }
  for(ProbeName in lstProbe){
    IndexSet <- lsfProbeSet[[ProbeName]][,1:3]
    RawDataSet <- lsfProbeSet[[ProbeName]][,4:ncol(lsfProbeSet[[ProbeName]])]
    MatmpDataSet <- as.matrix(RawDataSet)
    NormDataSet <- normalize.quantiles.use.target(MatmpDataSet,QTarget,copy=TRUE)
    NormDataSet <- MatmpDataSet
    NormDataSet <- data.frame(NormDataSet)
    colnames(NormDataSet) <- colnames(lsfProbeSet[[ProbeName]])[4:ncol(lsfProbeSet[[ProbeName]])]
    cat(colnames(lsfProbeSet[[ProbeName]])[4:ncol(lsfProbeSet[[ProbeName]])])
    ProbeNormSet <- cbind(IndexSet,NormDataSet)

    ofGCT <- paste("../output/",ProbeName,"_QNorm_Thread_",nthread,".gct",sep="")
    WriteGCT(ProbeNormSet,ofGCT,fOverwrite=TRUE)
    ofGCT <- paste("../output/",ProbeName,"_Raw_Thread_",nthread,".gct",sep="")
    WriteGCT(lsfProbeSet[[ProbeName]],ofGCT)
  }
}
  
########################################
#           QNormgctPlate              #
########################################
QNormgctPlate <- function(lstPlateList,QTarget,Overwrite=FALSE){
  for(PlateName in lstPlateList){
    
    lstRawDataNames <- sub(".gct","",list.files("../output",pattern=paste(PlateName,"_[A-P][0-2][0-9].gct",sep="")))
##    lstRawDataNames <- lstRawDataNames[!is.na(lstRawDataNames)]
    if(length(lstRawDataNames)==0) return(NA)
    ifGCT <- paste("../output/",lstRawDataNames[1],".gct",sep="")
    tmpGCT <- read.delim(paste("../output/",ifGCT,sep=""),sep="\t",
                         header = TRUE,skip = 2,
                         stringsAsFactors = FALSE)
    ProbeName <- ParseDataName(ifGCT)["ShortModeID",]
    ProbeSet <- subset(tmpGCT, select=Gene:Description)

    ofGCT <- paste("../output/",PlateName,"_",ProbeName,".gct",sep="")
    ofGCTraw <- paste("../output/",PlateName,"_",ProbeName,"_Raw.gct",sep="")
    if (file.exists(ofGCT)) {
      if (Overwrite){
        file.remove(ofGCT)
        file.remove(ofGCTraw)
      }  else{
        next
      }
    }


    i <- 0
    n <- length(lstRawDataNames)
    cat ("\tProcessing \t",PlateName, "\n")
    for (RawDataName in lstRawDataNames){
      i <- i + 1 
      ifGCT <- paste("../output/",RawDataName,".gct",sep="")
      tmpGCT <- read.delim(ifGCT,sep="\t",
                         header = TRUE,skip = 2,
                         stringsAsFactors = FALSE)
      tmpDataSet<- data.frame(
                     tmpGCT[match(ProbeSet$Gene,tmpGCT$Gene),
                            RawDataName],
                     stringsAsFactors = FALSE)
      colnames(tmpDataSet) <- RawDataName
      ProbeSet <- cbind(ProbeSet,tmpDataSet)
    }

    tmpNA <- is.na(ProbeSet[,4:ncol(ProbeSet)])
    colBad <- colnames(tmpNA)[apply(tmpNA,2,sum)>50]
    colGood <- setdiff(colnames(ProbeSet[4:ncol(ProbeSet)]),colBad)
    MAnormProbe <- normalize.quantiles.use.target(as.matrix(ProbeSet[,colGood]),QTarget,copy=TRUE)
    normDataSet <- cbind(ProbeSet[,1:3],
                         MAnormProbe)
    colnames(normDataSet)[4:ncol(normDataSet)] <- colGood
     
    WriteGCT(normDataSet,ofGCT,fOverwrite=TRUE)
    WriteGCT(ProbeSet,ofGCTraw)
  }
}


########################################
#            QNormgct                  #
########################################
QNormgct <- function(lstRawDataNames,gctFileName){
  lstRawDataNames <- sub(".gct","",list.files("../output",pattern="_[A-P][0-2][0-9].gct"))
  gctFileName <- "QNorm_06112012"
  i <- 1
  n <- length(lstRawDataNames)
  for (RawDataName in lstRawDataNames){
    cat ("Processing \t", i, " of ",n, ": RawDataName","\n")
    ifGCT <- paste("../output/",RawDataName,".gct",sep="")
    tmpGCT <- read.delim(ifGCT,sep="\t",
                         header = TRUE,skip = 2,
                         stringsAsFactors = FALSE)
    if (length(grep("DUO52HI53LO",RawDataName))){
      if (!exists( "DataSet5253")) {
        DataSet5253 <- tmpGCT
      } else {
        tmpDataSet <- data.frame(
                                 tmpGCT[match(DataSet5253$Gene,tmpGCT$Gene),
                                        RawDataName],
                                 stringsAsFactors = FALSE)
        colnames(tmpDataSet) <- RawDataName
        DataSet5253 <- cbind(DataSet5253,tmpDataSet)
      }
    } else if (length(grep("DUO53HI52LO",RawDataName))) {
      if (!exists( "DataSet5352")) {
        DataSet5352 <- tmpGCT
      } else {
        tmpDataSet <- data.frame(
                                 tmpGCT[match(DataSet5352$Gene,tmpGCT$Gene),
                                        RawDataName],
                                 stringsAsFactors = FALSE)
        colnames(tmpDataSet) <- RawDataName
        DataSet5352 <- cbind(DataSet5352,tmpDataSet)
      }
    } else if (length(grep("DUO44HI45LO",RawDataName))) {
      if (!exists( "DataSet4445")) {
        DataSet4445 <- tmpGCT
      } else {
        tmpDataSet <- data.frame(
                                 tmpGCT[match(DataSet4445$Gene,tmpGCT$Gene),
                                        RawDataName],
                                 stringsAsFactors = FALSE)
        colnames(tmpDataSet) <- RawDataName
        DataSet4445 <- cbind(DataSet4445,tmpDataSet)
      }
    } else if (length(grep("DUO45HI44LO",RawDataName))
               | length(grep("DUO45HI44LO",RawDataName))) {
      if (!exists( "DataSet4544")) {
        DataSet4544 <- tmpGCT
      } else {
        tmpDataSet <- data.frame(
                                 tmpGCT[match(DataSet4544$Gene,tmpGCT$Gene),
                                        RawDataName],
                                 stringsAsFactors = FALSE)
        colnames(tmpDataSet) <- RawDataName
        DataSet4544 <- cbind(DataSet4544,tmpDataSet)
      }
    }
    i <- i+1
  }

  # Processing probe set 52/53
  oDataSet5253 <- DataSet5253
  tmpNA <- is.na(DataSet5253[,4:ncol(DataSet5253)])
  colBad <- colnames(tmpNA)[apply(tmpNA,2,sum)>50]
  colGood <- setdiff(colnames(DataSet5253[4:ncol(DataSet5253)]),colBad)
  oDataSet5253 <- cbind(DataSet5253[,1:3],
                        normalize.quantiles(as.matrix(DataSet5253[,colGood])))

  tmp <- as.matrix(DataSet5253[,colGood])
  tmpNA <- is.na(tmp)
  tmp <- normalize.quantiles(as.matrix(DataSet5253[,colGood]))
  
  ofGCT <- paste("../output/",gctFileName,"_5253.gct",sep="")
  if (file.exists(ofGCT)) file.remove(ofGCT)
  gctFile <- file(ofGCT,open = "a")
  cat(paste("#1.2","\n",sep=""),
      file = gctFile,append=TRUE)
  cat(paste(nrow(oDataSet5253),"\t",ncol(oDataSet5253)-3,"\n",sep=""),
      file = gctFile,append=TRUE)
  close(gctFile)
  write.table(oDataSet5253,file=ofGCT,append=TRUE,sep="\t")

  ofGCT <- paste("../output/",gctFileName,"_5253_Raw.gct",sep="")
  if (file.exists(ofGCT)) file.remove(ofGCT)
  gctFile <- file(ofGCT,open = "a")
  cat(paste("#1.2","\n",sep=""),
      file = gctFile,append=TRUE)
  cat(paste(nrow(DataSet5253),"\t",ncol(DataSet5253)-3,"\n",sep=""),
      file = gctFile,append=TRUE)
  close(gctFile)
  write.table(DataSet5253,file=ofGCT,append=TRUE,sep="\t")
  # Processing probe set 53/52
  oDataSet5352 <- DataSet5352
  tmpNA <- is.na(DataSet5352[,4:ncol(DataSet5352)])
  colBad <- colnames(tmpNA)[apply(tmpNA,2,sum)>50]
  colGood <- setdiff(colnames(DataSet5352[4:ncol(DataSet5352)]),colBad)
  oDataSet5352 <- cbind(DataSet5352[,1:3],
                        normalize.quantiles(as.matrix(DataSet5352[,colGood])))

  tmp <- as.matrix(DataSet5352[,colGood])
  tmpNA <- is.na(tmp)
  tmp <- normalize.quantiles(as.matrix(DataSet5352[,colGood]))
  
  ofGCT <- paste("../output/",gctFileName,"_5352.gct",sep="")
  if (file.exists(ofGCT)) file.remove(ofGCT)
  gctFile <- file(ofGCT,open = "a")
  cat(paste("#1.2","\n",sep=""),
      file = gctFile,append=TRUE)
  cat(paste(nrow(oDataSet5352),"\t",ncol(oDataSet5352)-3,"\n",sep=""),
      file = gctFile,append=TRUE)
  close(gctFile)
  write.table(oDataSet5352,file=ofGCT,append=TRUE,sep="\t")

  ofGCT <- paste("../output/",gctFileName,"_5352_Raw.gct",sep="")
  if (file.exists(ofGCT)) file.remove(ofGCT)
  gctFile <- file(ofGCT,open = "a")
  cat(paste("#1.2","\n",sep=""),
      file = gctFile,append=TRUE)
  cat(paste(nrow(DataSet5352),"\t",ncol(DataSet5352)-3,"\n",sep=""),
      file = gctFile,append=TRUE)
  close(gctFile)
  write.table(DataSet5352,file=ofGCT,append=TRUE,sep="\t")
  # Processing probe set 44/45
  oDataSet4445 <- DataSet4445
  tmpNA <- is.na(DataSet4445[,4:ncol(DataSet4445)])
  colBad <- colnames(tmpNA)[apply(tmpNA,2,sum)>50]
  colGood <- setdiff(colnames(DataSet4445[4:ncol(DataSet4445)]),colBad)
  oDataSet4445 <- cbind(DataSet4445[,1:3],
                        normalize.quantiles(as.matrix(DataSet4445[,colGood])))

  ofGCT <- paste("../output/",gctFileName,"_4445.gct",sep="")
  if (file.exists(ofGCT)) file.remove(ofGCT)
  gctFile <- file(ofGCT,open = "a")
  cat(paste("#1.2","\n",sep=""),
      file = gctFile,append=TRUE)
  cat(paste(nrow(oDataSet4445),"\t",ncol(oDataSet4445)-3,"\n",sep=""),
      file = gctFile,append=TRUE)
  close(gctFile)
  write.table(oDataSet4445,file=ofGCT,append=TRUE,sep="\t")

  ofGCT <- paste("../output/",gctFileName,"_4445_Raw.gct",sep="")
  if (file.exists(ofGCT)) file.remove(ofGCT)
  gctFile <- file(ofGCT,open = "a")
  cat(paste("#1.2","\n",sep=""),
      file = gctFile,append=TRUE)
  cat(paste(nrow(DataSet4445),"\t",ncol(DataSet4445)-3,"\n",sep=""),
      file = gctFile,append=TRUE)
  close(gctFile)
  write.table(DataSet4445,file=ofGCT,append=TRUE,sep="\t")
  # Processing probe set 45/44
  oDataSet4544 <- DataSet4544
  tmpNA <- is.na(DataSet4544[,4:ncol(DataSet4544)])
  colBad <- colnames(tmpNA)[apply(tmpNA,2,sum)>50]
  colGood <- setdiff(colnames(DataSet4544[4:ncol(DataSet4544)]),colBad)
  oDataSet4544 <- cbind(DataSet4544[,1:3],
                        normalize.quantiles(as.matrix(DataSet4544[,colGood])))

  ofGCT <- paste("../output/",gctFileName,"_4544.gct",sep="")
  if (file.exists(ofGCT)) file.remove(ofGCT)
  gctFile <- file(ofGCT,open = "a")
  cat(paste("#1.2","\n",sep=""),
      file = gctFile,append=TRUE)
  cat(paste(nrow(oDataSet4544),"\t",ncol(oDataSet4544)-3,"\n",sep=""),
      file = gctFile,append=TRUE)
  close(gctFile)
  write.table(oDataSet4544,file=ofGCT,append=TRUE,sep="\t")

  ofGCT <- paste("../output/",gctFileName,"_4544_Raw.gct",sep="")
  if (file.exists(ofGCT)) file.remove(ofGCT)
  gctFile <- file(ofGCT,open = "a")
  cat(paste("#1.2","\n",sep=""),
      file = gctFile,append=TRUE)
  cat(paste(nrow(DataSet4544),"\t",ncol(DataSet4544)-3,"\n",sep=""),
      file = gctFile,append=TRUE)
  close(gctFile)
  write.table(DataSet4544,file=ofGCT,append=TRUE,sep="\t")
}

########################################
#               FindDEG                #
########################################
FindDEG <- function(gctFileName,dfControlNames){
  ifGCT <- paste("../output/",gctFileName,".gct",sep="")
  tmpGCT <- read.delim(ifGCT,
                       header = TRUE,skip=2,
                       stringsAsFactors = FALSE)
  for (Sample in unique(dfControlNames$Sample)) {
    cat(paste(Sample,"\n"))
    for (Control in dfControlNames$Control[Sample == dfControlNames$Sample]) {
      cat(paste("\t",Control,"\n"))
    }
  }
  
  idxSample = c(as.character(lstOEDataNames))
  idxControls = dfControlNames$Control[
    (idxSample[1] == dfControlNames$Sample |
     idxSample[1] == dfControlNames$Sample)]
  ## c("KDA006_A375_96H_X1_B3_DUO52HI53LO_G19",
  ##   "KDA006_A375_96H_X1_B3_DUO52HI53LO_L04",
  ##   "KDA006_A375_96H_X1_B3_DUO52HI53LO_L07",
  ##   "KDA006_A375_96H_X1_B3_DUO52HI53LO_M10")
  ## tmpDataSet <- data.frame(

  ##                 Mean=apply(tmpGCT[,idxControls],1,mean),
  ##                 Sd=apply(tmpGCT[,idxControls],1,sd))
  ## DataSet <- cbind(tmpGCT,tmpDataSet)
  ## tmpTTest <- mat.or.vec(nrow(tmpDataSet),3)
  ## colnames(tmpTTest) <- c("t","p","p_adj")
  ## rownames(tmpTTest) <- tmpGCT$Gene
  ## for (i in 1:nrow(tmpDataSet)){
  ##   tmp <- t.test(tmpGCT[i,idxControls],mu=tmpGCT[i,idxSample])
  ##   tmpTTest[i,"t"] <- tmp$statistic
  ##   tmpTTest[i,"p"] <- tmp$p.value
  ## }
  ## tmpTTest[,"p_adj"] <- p.adjust(tmpTTest[,"p"])

  tmpLFC <- tmpGCT
  tmpLFC[,c(idxSample,idxControls)] <-
    log2(tmpGCT[,c(idxSample,idxControls)])

  tmpLFCDataSet <- data.frame(
              Mean=apply(tmpLFC[,idxControls],1,mean),
              Sd=apply(tmpLFC[,idxControls],1,sd))
  LFCDataSet <- cbind(tmpLFC,tmpLFCDataSet)
  tmpLFCTTest <- mat.or.vec(nrow(tmpLFCDataSet),3)
  colnames(tmpLFCTTest) <- c("t","p","p_adj")
  rownames(tmpLFCTTest) <- tmpGCT$Gene
  for (i in 1:nrow(tmpLFCDataSet)){
#    tmp <- t.test(tmpLFC[i,idxSample],y=tmpLFC[i,idxControl])
    tmp <- t.test(tmpLFC[i,idxControls],mu=tmpLFC[i,idxSample])
    tmpLFCTTest[i,"t"] <- tmp$statistic
    tmpLFCTTest[i,"p"] <- tmp$p.value
  }
  tmpLFCTTest[,"p_adj"] <- p.adjust(tmpLFCTTest[,"p"])

  idxDEG <- tmpLFCTTest[,"p_adj"] <= 0.1
  cbind(rownames(tmpLFCTTest)[idxDEG],
        LFCDataSet[idxDEG,"Description"],
        tmpLFCTTest[idxDEG,])
  
  tmpTTest <- data.frame(
                t=numeric(0),
                p=numeric(0),
                padj=numeric(0)
                )
}
## tmpDataSet <- data.frame(
##                    tmpGCT[match(DataSet$Gene,tmpGCT$Gene),
##                           RawDataName])
## colnames(tmpDataSet) <- RawDataName
## DataSet <- cbind(DataSet,tmpDataSet)

########################################
#                QN2LFC                #
########################################
QN2LFC <- function(PlateName,Overwrite=FALSE){
  ifPlateMap <- "../Annotation/plate_maps.csv" 
 ## ifPlateMap <- "../ljp001/total.csv"
  ProbeName <- ParseDataName(paste(PlateName,"_A01"))["ShortModeID",]
  LFCFileName <- paste(PlateName,"_",ProbeName,"_LFC",sep="")

  if (file.exists(LFCFileName)) {
    if (Overwrite)
      file.remove(LFCFileName)
    else {
      cat("\t\t",LFCFileName," exists. Skip. \n")
      return(NA)
    }
  }

  PlateMap <- read.csv(ifPlateMap, header = TRUE, stringsAsFactors = FALSE)
  Type <- substr(ParseDataName(paste(PlateName,"_A01"))["Plate",1],1,2)
  if(Type != "KD"& Type!="OE" & Type !="CP" & Type!="LJ"&&Type!="LI") return(NA)
##  KDNames <- FindKDNames(PlateName)
  QNData <- ReadQNData(PlateName)
  if (is.null(QNData)) {
    return (NA)
  }  
  if (Type=="KD") {
    lstControlWells <- FindKDControls(colnames(QNData)[4],PlateMap)
  } else if(Type=="OE"){
    lstControlWells <- FindOEControls(colnames(QNData)[4],PlateMap)
  } else if(Type=="CP"){
    lstControlWells <- FindDrugControls(colnames(QNData)[4],PlateMap)
  } else if(Type=="LJ"){
    lstControlWells <- FindDrugControlsbac(colnames(QNData)[4],PlateMap)
  } else if(Type =="LI"){
    if(substr(ParseDataName(paste(PlateName,"_A01"))["Plate",1],1,4) =="LITM"){
      lstControlWells <- FindKDControls(colnames(QNData)[4],PlateMap)}
  } 
  
  QN2LFC <- QNData
  ControlList <- na.omit(match(lstControlWells$Control,colnames(QNData)))
  if(length(ControlList)==0){
    cat("No Control Well available ",PlateName,"\n")
    return(NA)}
  QNControl <- QNData[ControlList]
  QN2LFC[,4:ncol(QN2LFC)] <- QNData[,4:ncol(QNData)]-
    apply(QNControl,1,mean,na.rm=TRUE)
  WriteGCT(QN2LFC,LFCFileName,fOverwrite = Overwrite)
}

########################################
#               FindKDNames            #
########################################
FindKDNames <- function(PlateName) {
  ifControlNames <- "../Annotation/controls.csv"
  ControlNames <- unique(read.csv(ifControlNames, header = TRUE, stringsAsFactors = FALSE)
                         $PerturbagenDesc)
  ifPlateMap <- "../Annotation/plate_maps.csv"
  PlateMap <- read.csv(ifPlateMap, header = TRUE, stringsAsFactors = FALSE)
  ifPlateList <- "../Annotation/plate_list.csv"
  PlateList <- read.csv(ifPlateList, header = TRUE,
                        stringsAsFactors = FALSE, row.names="PlateName")

  QNData <- ReadQNData(PlateName)
  if (is.null(QNData)) return (NA)
  
  lstQNWells <- as.list(ParseDataName(colnames(QNData)[4:ncol(QNData)])["Well",])
  lstWell <- paste(rep(LETTERS[1:16],each=24),
                    formatC(rep(1:24,times=16),width=2,flag="0"), sep="")
  
  this.PlateMap <- PlateMap[ParseDataName(colnames(QNData)[4])["Plate",1] ==
                            PlateMap$PerturbagenPlate,]
  rownames(this.PlateMap) <- this.PlateMap[,"WellName"]
  tmpControl <- data.frame(Control=rep(FALSE,nrow(this.PlateMap)))
  this.PlateMap <- cbind(this.PlateMap,tmpControl)

  for (ControlName in ControlNames) 
    this.PlateMap[ControlName==this.PlateMap$PerturbagenDesc,"Control"] <- TRUE
  tmpDataName <- data.frame(DataName=colnames(QNData)[4:ncol(QNData)],
                            stringsAsFactors = zFALSE)
  this.QNPlateMap <- cbind(this.PlateMap[unlist(lstQNWells),],tmpDataName)
  
  return(this.QNPlateMap)
}

########################################
#               ReadQNData             #
########################################
ReadQNData <- function(PlateName){
  ifQNData <- list.files("../output", pattern=paste(PlateName,"_[4-5][2-5][4-5][2-5].gct",sep=""))
  if (!1==length(ifQNData)) {
    cat("Error: there are", length(ifQNData), "matched QNorm gct file for", PlateName,".\n")
    cat ("\tWhich are", ifQNData,"\n")
    return (NULL)
  }
  return(ReadGCT(ifQNData))
}

########################################
#               ReadGCT                #
########################################
ReadGCT <- function(DataName){
  if (length(grep(".gct$",DataName))) ifGCT <- paste("../output/",DataName,sep="")
  else ifGCT <- paste("../output/",DataName,".gct",sep="")
  if (!file.exists(ifGCT)) {
    cat("Error: there is not GCT file named ", ifGCT,".\n")
    return (NA)
  }
  return (read.delim(paste("../output/",ifGCT, sep=""),
                     sep="\t",
                     header = TRUE,skip = 2,
                     stringsAsFactors = FALSE))
}

########################################
#              WriteGCT                #
########################################
WriteGCT <- function(GCTData,gctFileName,fOverwrite=FALSE){
  if (length(grep(".gct$",gctFileName)))
    ofGCT <- paste("../output/",gctFileName,sep="")
  else
    ofGCT <- paste("../output/",gctFileName,".gct",sep="")
  if (file.exists(ofGCT)) {
    if (fOverwrite)
      file.remove(ofGCT)
    else {
##      cat(ofGCT," exists. Skip. \n")
      return()
    }
  }

  gctFile <- file(ofGCT,open = "a")
  cat(paste("#1.2","\n",sep=""),
      file = gctFile,append=TRUE)
  cat(paste(nrow(GCTData),"\t",ncol(GCTData)-3,"\n",sep=""),
      file = gctFile,append=TRUE)
  close(gctFile)
  suppressWarnings(write.table(GCTData,file=ofGCT,append=TRUE,sep="\t"))
}

########################################
#             ParseDataName            #
########################################
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
  ShortModeID <- ParseDataName["DetectMode",]
  ShortModeID <- replace(ShortModeID,grep(44,ShortModeID),4445)
  ShortModeID <- replace(ShortModeID,grep(52,ShortModeID),5253)
  ParseDataName <- rbind(ParseDataName,ShortModeID)
  return (ParseDataName)
}
########################################
#            ModeShortName             #
########################################
ModeShortName <- function(DetectMode){
  ModeSName <- DetectMode
  for(i in 1:length(ModeSName)){
    if(grep(44,ModeSName[i])==1)
      ModeSName[i] <- 4445
    else
      ModeSName[i] <- 5253
  }
  return(ModeSName)
}

