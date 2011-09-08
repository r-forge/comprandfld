\name{CheckInput}
\alias{CheckInput}

\title{Checks the Inserted Input}
\description{
  Subroutine called by the fitting procedures.
  The procedure controls the input inserted into
  the fitting procedures.
}
\usage{
CheckInput(coordx, coordy, coordt, corrmodel, data, fixed,
           grid, likelihood, lonlat, margins, maxdist,
           maxtime, model, optimizer, replicates, start,
           taper, threshold, type, varest, vartype,
           weighted, weights, winconst)
}
\arguments{
   \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of points) assigning 2-dimensions of coordinates or a numeric vector assigning
    1-dimension of coordinates.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector otherwise it will be ignored.}
  \item{coordt}{A numeric vector assigning 1-dimension of
    temporal coordinates.}
  \item{corrmodel}{String; the name of a correlation model, for the
    description see \code{\link[RandomFields]{CovarianceFct}}.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{fixed}{A named list giving the values of the parameters that
    will be considered as known values. The listed parameters for a
    given correlation function will be not estimated, i.e. if
    \code{list(nugget=0)} the nugget effect is ignored.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{likelihood}{String; the configuration of the composite
    likelihood. \code{Marginal} is the default.}
  \item{lonlat}{Logical; if \code{FALSE} (the default), \code{coordx}
    and \code{coordy} are interpreted as Cartesian coordinates
    otherwise they are considered as longitude and latitude.}
  \item{margins}{String; the type of the marginal distribution of the
    max-stable field.}
  \item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial distance considered in the composite-likelihood
    computation.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal lag separation in the composite-likelihood.}
  \item{model}{String; the density associated to the likelihood
    objects. \code{Gaussian} is the default.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). 'Nelder-Mead' is the
    default.}
  \item{replicates}{Logical; if \code{FALSE} (the default) one spatial random field is
    considered, instead if \code{TRUE} the data are considered as iid replicates of a field.}
  \item{start}{A named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure. \code{NULL} is the default.}
  \item{taper}{String; the name of the tapered correlation
    function.}
  \item{threshold}{Numeric; a value indicating a threshold for the
    binary random field.}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
  \item{varest}{Logical; if \code{TRUE} the estimate' variances and
    standard errors are returned.
    \code{FALSE} is the default.}
  \item{vartype}{String; the type of estimation method for computing
  the estimate variances, see \code{\link{FitComposite}}.}
  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted. If \code{FALSE} (the
    default) the composite likelihood is not weighted.}
  \item{weights}{A numeric vector of weights.}
  \item{winconst}{Numeric; a positive real value -- if
    \code{vartype=Sub-Samp} -- that determines the window size in the sub-sampling
    estimates of the variances, see \code{\link{FitComposite}}.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibg.it},
  \url{http://www.unibg.it/pers/?simone.padoan}.}

\keyword{Composite}
