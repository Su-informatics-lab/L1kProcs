\name{libEGEM}
\alias{libEGEM}
\title{libEGEM}
\description{
  EGEM matrix based on LINCS data. This can be used for \code{csNMF} function of
  L1KProcs package.}
\usage{
data(libEGEM)
}
\format{
  A list with egem matrix information based on different cell lines and
  time for measurements. \cr
  Each item of libEGEM is a list with three items:\cr
  egem    matrix  egem matrix n genes X m compounds\cr
  CNames  vector  compound perturbagen description\cr
  PNames  vector  genetic pertrubagen deescription,like target genes.
}
\examples{
  data(libEGEM)
  names(libEGEM)
  str(libEGEM[[1]])
}
\keyword{datasets}
