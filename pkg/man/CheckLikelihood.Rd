\name{CheckLikelihood}
\alias{CheckLikelihood}

\title{Check of the type of Composite-likelihood}
\description{
  Subroutine called by InitParam.
  The procedure controls the type of the composite-likelihood passed to the FitComposite procedure.
}
\usage{
CheckLikelihood(likelihood)
} 
\arguments{
  \item{likelihood}{String; the configuration of the composite
    likelihood. \code{Marginal} is the default.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibg.it},
  \url{http://www.unibg.it/pers/?simone.padoan}.}

\keyword{Composite}
