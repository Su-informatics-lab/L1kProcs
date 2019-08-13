PipeLine <- function(iStart,iEnd,Path,nthread,Overwrite){
  setwd(Path)
  suppressMessages(library(preprocessCore))
  suppressMessages(library(stats))
  suppressMessages(library(prada))
  source("Util.R")
  ## Files & Directories:
  ##   ../Annotation/bead_gene_map.csv
  ##   ../Annotation/HG_U133A.chip
  ##   ../Annotation/plate_maps.csv
  ##   ../Annotation/plate_list.csv
  ##   ../Output/dataName.csv
  ##   ../Output/dataName.gct

  ifBead2Gene <- "../Annotation/bead_gene_map.csv"
  ifU133A <- "../Annotation/HG_U133A.chip"
  ifPlateMap <- "../Annotation/plate_maps.csv"
  ifPlateList <- "../Annotation/plate_list.csv"
  Bead2Gene <- read.csv(ifBead2Gene, header = TRUE, stringsAsFactors = FALSE)
  PlateMap <- read.csv(ifPlateMap, header = TRUE, stringsAsFactors = FALSE)
  U133A <- read.delim(ifU133A, header = TRUE, stringsAsFactors = FALSE)
  rownames(U133A) <- U133A$pr_id


  if (!file.exists("lstPlates.RData")){
    cat("Please save all the plate names to lstPlates.Data first!\n")
    return(NA)
  }

  load(file = "lstPlates.RData")
  load(file="QTarget.RData")
  QTarget <- as.vector(QTarget)

  QNormgctPlate(lstPlates[iStart:iEnd],QTarget,Overwrite)

  for (PlateName in lstPlates[iStart:iEnd]){
    QN2LFC(PlateName, Overwrite)
  }  
}
