\name{CorrelationFct}
\alias{CorrelationFct}

\title{Computation of the Correlation function}
\description{
  Subroutine called by Covariogram.
  The procedure computes the estimated correlation function for a given
  fitted model.
}
\usage{
CorrelationFct(corrmodel, lags, numpairs, param)
} 
\arguments{
  \item{corrmodel}{Numeric; the id of the correlation model.}
  \item{lags}{A numeric vector of distances between points.}
  \item{numpairs}{Numeric; the number of distances.}
  \item{param}{A numeric vector with the parameter values.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}.}

\keyword{Composite}
