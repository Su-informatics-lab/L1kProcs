rm(list=ls())
library(prada)
#library(lattice)
library(stats)
#library(cluster)
#library(sigclust)
#library(utils)
#library(lga)
#library(Rmpi)

setwd("C:/Users/tmhjxs20/Documents/Proposals/LINCS/L1K/Software")

infile <- "../Data/KDA006_A375_96H_X1_B3_DUO52HI53LO/KDA006_A375_96H_X1_B3_DUO52HI53LO_M19.lxb";
outfile <-  "../Output/KDA006_A375_96H_X1_B3_DUO52HI53LO_M19.csv";
outdatafile <-  "../Output/KDA006_A375_96H_X1_B3_DUO52HI53LO_M19.RData";
cat ("From", infile,"to",outfile,"\n")

dat<-readFCS(infile)

data <- dat@exprs

rawlfc <- cbind(data[,"RID"], log2(data[,"RP1"]))
colnames(rawlfc) <- c("RID","EXP")
lfc <- rawlfc
lfc <- rawlfc[!is.infinite(lfc[,"EXP"]),]

results <- mat.or.vec(500,7)
colnames(results) <- c("AnalyteN","n1","c1","n2","c2", "n_peak","n_tot")
results[,"AnalyteN"] <- 1:500

GMM <- function(inX, inP) {
  lambda <- 1.25/2
  return(lambda*pnorm(inX,mean=inP[1],sd=inP[2])
         + (1-lambda)*pnorm(inX,mean=inP[3],sd=inP[4]))
}

FGMM <- function(inP,inX,inY) sum((inY-GMM(inX,inP))^2)

lambda <- 1.25/2
for (i in 1:500){
  cat("Processing", i, "of 500... ")
  Analyte <- lfc[i == lfc[,"RID"],"EXP"]
  results[i,"n_tot"] <- length(Analyte)
  if ( i <= 10) {
    results[i,"n_peak"] <- 1
    results[i,"n1"] <- length(Analyte)
    results[i,"c1"] <- median(Analyte)
    results[i,"n2"] <- NA
    results[i,"c2"] <- NA
  } else {
    results[i,"n_peak"] <- 2
    HistAnalyte <- hist(Analyte,breaks=500,plot=FALSE)
    CDFAnalyte <- cumsum(HistAnalyte$counts/sum(HistAnalyte$counts))
    fit <- optim(c(2*median(Analyte)/3,sd(Analyte),4*median(Analyte)/3,sd(Analyte)),
                 FGMM,NULL,HistAnalyte$mids,CDFAnalyte,control=list(maxit=10000))
    results[i,"n1"] <- round(results[i,"n_tot"]*lambda)
    results[i,"n2"] <- results[i,"n_tot"]-results[i,"n1"]
    results[i,"c1"] <- fit$par[1]
    results[i,"c2"] <- fit$par[3]
  }
  cat("Done! \n")
}

save(results,dat,rawlfc,lfc,file=outdatafile,compress=TRUE)
write.csv(results,file=outfile)
