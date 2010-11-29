\name{GevLogLik}
\alias{GevLogLik}
\title{Log-Likelihood of the Generalized Extreme Value
  Distribution}
\description{
  the function returns the log-likelihood value of the Generalized
  Extreme Value Distribution for a given set of data and parameters.
}
\usage{
GevLogLik(data, numdata, param)
} 
\arguments{
  \item{data}{A vector of extreme values.}
  \item{numdata}{The number of data observations.}
  \item{param}{The vector of GEV parameters (location, scale and shape).}
}
\value{
  The log-likelihood value is returned.}

\references{
  de Haan, L. and Ferreira, A. (2006) \emph{Extreme Value Theory An
    Introduction}. Springer Verlang, New York.
}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://people.epfl.ch/simone.padoan}}

\keyword{extremes}
