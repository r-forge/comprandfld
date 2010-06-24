\name{CompLikelihood}
\alias{CompLikelihood}

\title{Computation of the Composite-likelihood}
\description{
  Subroutine called by OptimCompLik.
  The procedure computes the composite-likelihood for a given set of
  data and parameters.
}
\usage{
CompLikelihood(coordx, coordy, corrmodel, data, dista, fixed, lags,
               model, namescorr, namesnuis, numcoord, numdata, param,
               type)
} 
\arguments{
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of points) assigning 2-dimensions of coordinates or a numeric vector assigning
    1-dimension of coordinates.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector otherwise it will be ignored.}
  \item{corrmodel}{Numeric; the id of the correlation model.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{dista}{A numeric value specifying the cut-off distance.}
  \item{fixed}{A numeric vector with the parameters that
    will be considered as known values.}
  \item{lags}{A numeric vector of distances between points.}
  \item{model}{Numeric; the id of the random field.}
  \item{namescorr}{String; the names of the correlation parameters.}
  \item{namesnuis}{String; the names of the nuisance parameters.}
  \item{numcoord}{Numeric; the number of coordinates.}
  \item{numdata}{Numeric; the number of the data in time.}
  \item{param}{A numeric vector with the parameter values.}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}.}

\keyword{Composite}
