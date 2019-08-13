.First.lib <-function(pkg,lib){
  ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"),"Version")       
  ver <- as.character(ver)
  library.dynam("L1KAnno",pkg,lib)
}

