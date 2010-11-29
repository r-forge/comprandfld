\name{CorrelationFct}
\alias{CorrelationFct}

\title{Computes the Correlation function}
\description{
  Subroutine called by Covariogram.
  The procedure computes the estimated correlation function for a given
  fitted model.
}
\usage{
CorrelationFct(corrmodel, lags, param)
} 
\arguments{
  \item{corrmodel}{Numeric; the id of the correlation model.}
  \item{lags}{A numeric vector of distances between points.}
  \item{param}{A numeric vector with the parameter values.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://people.epfl.ch/simone.padoan}.}

\keyword{Composite}
