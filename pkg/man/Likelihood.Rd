\name{Likelihood}
\alias{Likelihood}

\title{Optimizes the Log Likelihood}
\description{
  Subroutine called by FitComposite.
  The procedure estimates the model parameters by maximization of the
  log-likelihood.
}
\usage{
Likelihood(corrmodel, data, fixed, flagcor, flagnuis, grid, lower,
           model, namescorr, namesnuis, namesparam, numcoord, numpairs,
           numparamcor, numrep, numtime, optimizer, param, setup,
           spacetime, varest, taper, type, upper)
}
\arguments{
  \item{corrmodel}{Numeric; the id of the correlation model.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{flagcor}{A numeric vector of flags denoting which correlation
    parameters have to be estimated.}
  \item{flagnuis}{A numeric verctor of flags denoting which nuisance
    parameters have to estimated.}
  \item{fixed}{A numeric vector of parameters that
    will be considered as known values.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{lower}{A numeric vector with the lower bounds of the parameters' ranges.}
  \item{model}{Numeric; the id value of the density associated to the likelihood
    objects.}
  \item{namescorr}{String; the names of the correlation parameters.}
  \item{namesnuis}{String; the names of the nuisance parameters.}
  \item{namesparam}{String; the names of the parameters to be maximised.}
  \item{numcoord}{Numeric; the number of coordinates.}
  \item{numpairs}{Numeric; the number of pairs.}
  \item{numparamcor}{Numeric; the number of the correlation parameters.}
  \item{numrep}{Numeric; the number of iid replications.}
  \item{numtime}{Numeric; the number of temporal observations.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). 'Nelder-Mead' is the
    default.}
  \item{param}{A numeric vector of parameters.}
  \item{setup}{A List of useful components for the estimation based on the maximum tapered likelihood.}
  \item{spacetime}{Logical; if the random field is spatial
    (\code{FALSE}) or spatio-temporal (\code{TRUE}).}
  \item{varest}{Logical; if \code{TRUE} the estimate' variances and
    standard errors are returned.
    \code{FALSE} is the default.}
  \item{taper}{String; the name of the taper correlation function.}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
  \item{upper}{A numeric vector with the upper bounds of the parameters' ranges.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{http://staff.deuv.cl/bevilacqua/}.}


\keyword{Composite}
