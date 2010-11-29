\name{CheckParam}
\alias{CheckParam}

\title{Check of the Parameters}
\description{
  Subroutine called by InitParam.
  The procedure controls the validity of the correlation's parameters.
}
\usage{
CheckParam(corrmodel, namesparam, numparam)
} 
\arguments{
  \item{corrmodel}{String; the name of a correlation model.}
  \item{namesparam}{String; the names of the parameters.}
  \item{numparam}{Numeric; the number of the parameters.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://people.epfl.ch/simone.padoan}.}

\keyword{Composite}
