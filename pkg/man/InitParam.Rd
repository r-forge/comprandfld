\name{InitParam}
\alias{InitParam}

\title{Initializes the Parameters for the Fitting Procedures}
\description{
  Subroutine called by the fitting procedures.
  The procedure initializes the parameters for the fitting procedure.
}
\usage{
InitParam(coordx, coordy, coordt, corrmodel, data, fixed, grid,
          likelihood, lonlat, margins, maxdist, maxtime, model,
          parscale, paramrange, replicates, start, threshold,
          type, varest, vartype, weighted, winconst)
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
  \item{corrmodel}{String; the name of a correlation model.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{fixed}{A named list giving the values of the parameters that
    will be considered as known values.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{likelihood}{String; the configuration of the composite
    likelihood.}
  \item{lonlat}{Logical; if \code{FALSE} (the default), \code{coordx}
    and \code{coordy} are interpreted as Cartesian coordinates
    otherwise they are considered as longitude and latitude.}
  \item{margins}{String; the type of the marginal distribution of the
    max-stable field.}
  \item{model}{String; the density associated to the likelihood
    objects. \code{Gaussian} is the default.}
  \item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial distance considered in the composite-likelihood
    computation.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal lag considered in the composite-likelihood
    computation.}
  \item{parscale}{A numeric vector of scaling factor to improve the
    maximizing procedure, see \code{\link{optim}}.}
  \item{paramrange}{A numeric vector of parameters ranges, see
    \code{\link{optim}}.}
  \item{replicates}{Logical; if \code{FALSE} (the default) one spatial random field is
    considered, instead if \code{TRUE} the data are considered as iid replicates of a field.}
   \item{start}{A named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure.}
   \item{threshold}{Numeric; a value indicating a threshold for the
    binary random field.}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
  \item{varest}{Logical; if \code{TRUE} the estimates' variances and
    standard errors are returned.
    \code{FALSE} is the default.}
  \item{vartype}{String; the type of estimation method for computing
    the estimate variances, see the Section \bold{Details}.}
  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted, see \code{\link{FitComposite}}.}
  \item{winconst}{Numeric; a positive real value indicating the window size used from the
   sub-sampling method for the estimation of the parameters' variances.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibg.it},
  \url{http://www.unibg.it/pers/?simone.padoan}.}

\keyword{Composite}
