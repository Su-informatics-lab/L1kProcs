\name{libDEG}
\alias{libDEG}
\title{libDEG}
\description{
  differential expressed genes after shRNA genetic perturbagens
  based on LINCS data. It can be used in \code{egem} function of
  L1KProcs package.
}
\usage{
data(libDEG)
}
\format{
  A list with differential expressed genes after shRNA genetic
  perturbagens based on LINCS data of different cell lines and
  time for measurements. \cr
  Each item of libDEG is a list with five items:\cr
  up        list      positions of up regulated genes in the \code{genelist}\cr
  down      list      positions of down regulated genes in the \code{genelist}\cr
  mean.perm vector    mean values of random egem scores\cr
  sd.perm   vector    standard deviation values of random egem scores\cr
  PtDesc    vector    the target genes of the shRNA perturbagens\cr
}
\examples{
  data(libDEG)
  names(libDEG)
  names(libDEG[[1]])
}
\keyword{datasets}
