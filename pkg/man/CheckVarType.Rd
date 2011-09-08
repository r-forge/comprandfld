\name{CheckVarType}
\alias{CheckVarType}

\title{Check of the method for the computation of the estimates' variances}
\description{
  Subroutine called by InitParam.
  The procedure controls the method used to compute the estimates' variances.
}
\usage{
CheckVarType(type)
} 
\arguments{
  \item{type}{String; the method used to compute the estimates'
    variances. If \code{SubSamp} (the default) the estimates' variances
    are computed by the sub-sampling method, see \code{\link{FitComposite}}.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibg.it},
  \url{http://www.unibg.it/pers/?simone.padoan}.}

\keyword{Composite}
