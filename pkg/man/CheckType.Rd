\name{CheckType}
\alias{CheckType}

\title{Check of the type of likelihood objects}
\description{
  Subroutine called by InitParam.
  The procedure controls the type of likelihood objects that form the composite-likelihood  .
}
\usage{
CheckType(type)
} 
\arguments{
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibg.it},
  \url{http://www.unibg.it/pers/?simone.padoan}.}

\keyword{Composite}
