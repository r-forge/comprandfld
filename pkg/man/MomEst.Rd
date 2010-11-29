\name{MomEst}
\alias{MomEst}
\title{Moment Estimator of the Extreme Value Distribution's parameters}
\description{
  Using the moment estimator, the function returns the parameter estimates
  of the generalized extreme value distribution for a given dataset of extreme values.
}
\usage{
MomEst(data, n)
} 
\arguments{
  \item{data}{A vector of extreme values.}
  \item{n}{The number of observations.}
}

\value{
  The returned object is a list with:
  
    \item{location}{The location estimate.}
    \item{scale}{The scale estimate.}
    \item{shape}{The shape estimate.}
}

\references{
  de Haan, L. and Ferreira, A. (2006) \emph{Extreme Value Theory An
    Introduction}. Springer Verlang, New York.
  
}

\seealso{\code{\link{GevLogLik}}, \code{\link{FitGev}}}


\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://people.epfl.ch/simone.padoan}}


\keyword{extremes}
