\name{Likelihood}
\alias{Likelihood}

\title{Computation of the likelihood}
\description{
  Subroutine called by OptimLik.
  The procedure computes the likelihood function for a given
  set of data and parameters.
}
\usage{
Likelihood(corrmodel, data, fixed, grid, lags, model,
           namescorr, namesnuis, numcoord, numdata,
           numpairs, param, type)
} 
\arguments{
  \item{corrmodel}{String; the name of a correlation model.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{fixed}{A named list giving the values of the parameters that
    will be considered as known values.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{lags}{A numeric vector of distances between points.}
  \item{model}{String; the density associated to the likelihood
    objects. \code{Gaussian} is the default.}
  \item{namescorr}{String; the names of the correlation parameters.}
  \item{namesnuis}{String; the names of the nuisance parameters.}
  \item{numcoord}{Numeric; the number of coordinates;}
  \item{numdata}{Numeric; the number of data replications in time.}
  \item{numpairs}{Numeric; the number of pairwise points.}
  \item{param}{A numeric vector of parameters' values.}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}.}

\keyword{Composite}
