\name{CheckModel}
\alias{CheckModel}

\title{Check of the type of Random Field}
\description{
  Subroutine called by InitParam.
  The procedure controls the type of random field passed to
  the fitting procedures.
}
\usage{
CheckModel(model)
} 
\arguments{
   \item{model}{String; the density associated to the likelihood
    objects. \code{Gaussian} is the default.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://people.epfl.ch/simone.padoan}.}

\keyword{Composite}
