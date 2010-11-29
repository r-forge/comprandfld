\name{LogNormDen}
\alias{LogNormDen}

\title{Computes the multivariate log-normal density}
\description{
  Subroutine called by the Likelihood procedure.
  The procedure compute the multivariate log-normal density
  for a given set of data and parameters.
}
\usage{
LogNormDen(stdata, detvarcov, ivarcov, numcoord, type)
} 
\arguments{
  \item{stdata}{A numeric vector (\eqn{d \times 1}{d x 1}) of data.}
  \item{detvarcov}{Numeric; the determinant of the variance-covariance
    matrix.}
  \item{ivarcov}{The inverse of the variance-covariance matrix (\eqn{d \times d}{d x d}).}
  \item{numcoord}{The number of point's coordinates.}
  \item{type}{The numeric id denoting the type of likelihood.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://people.epfl.ch/simone.padoan}.}

\keyword{Composite}
