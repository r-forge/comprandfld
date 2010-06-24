\name{OptimLik}
\alias{OptimLik}

\title{Optimization of the log-likelihood}
\description{
  Subroutine called by FitComposite.
  The procedure estimates the model parameters by maximisation of the 
  log-likelihood.
}
\usage{
OptimLik(corrmodel, data, fixed, grid, hessian, lags, lower,
         model, namescorr, namesnuis, namesparam, numcoord,
         numdata, numpairs, optimizer, param, varest, type,
         upper)
} 
\arguments{
  \item{corrmodel}{Numeric; the id of the correlation model.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{fixed}{A numeric vector of parameters that
    will be considered as known values.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{hessian}{Logical; if \code{FALSE} (the default) the hessian
    matrix is not computed from the optimizer routine.}
  \item{lags}{A numeric vector of distances between points.}
  \item{lower}{A numeric vector with the lower bounds of the parameters' ranges.}
  \item{model}{Numeric; the id value of the density associated to the likelihood
    objects.}
  \item{namescorr}{String; the names of the correlation parameters.}
  \item{namesnuis}{String; the names of the nuisance parameters.}
  \item{namesparam}{String; the names of the parameters to be maximised.}
  \item{numcoord}{Numeric; the number of coordinates.}
  \item{numdata}{Numeric; the number of data replications in time.}
  \item{numpairs}{Numeric; the number of pairwise distances.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). 'Nelder-Mead' is the
    default.}
  \item{param}{A numeric vector of parameters' values.}
  \item{varest}{Logical; if \code{TRUE} the estimate' variances and
    standard errors are returned. 
    \code{FALSE} is the default.}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
  \item{upper}{A numeric vector with the upper bounds of the parameters' ranges.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}.}

\keyword{Composite}
