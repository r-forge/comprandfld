\name{CheckVarType}
\alias{CheckVarType}

\title{Checking Variance Estimates Type}
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

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home}.}
\keyword{Composite}
