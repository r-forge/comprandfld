\name{WlsInit}
\alias{WlsInit}

\title{Computes Starting Values based on Weighted Least Squares}
\description{
  Subroutine called by FitComposite.
  The function returns opportune starting values for the
  composite-likelihood fitting procedure based on weigthed least squares.
}
\usage{
WlsInit(coordx, coordy, coordt, corrmodel, data, fcall, fixed, grid,
        likelihood, lonlat, margins, maxdist, maxtime, model, numblock,
        param, parscale, paramrange, replicates, start, threshold,
        type, varest, vartype, weighted, winconst, winstp)
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
    description.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{fcall}{String; "fitting" to call the fitting procedure and
    "simulation" to call the simulation procedure.}
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
  \item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial distance considered in the composite-likelihood
    computation.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal separation considered in the composite-likelihood
    computation.}
  \item{model}{String; the name of the model. Here the default is
    \code{NULL}.}
  \item{numblock}{Numeric; the observation size of the underlying random
    field. Only in case of max-stable random fields and in the
    simulation.}
  \item{param}{A numeric vector of parameter values required in the
    simulation procedure of random fields.}
  \item{parscale}{A numeric vector with scaling values for improving the
    maximisation routine.}
  \item{paramrange}{A numeric vector with the range of the parameter
    space.}
  \item{replicates}{Logical; if \code{FALSE} (the default) one spatial random field is
    considered, instead if \code{TRUE} the data are considered as iid replicates of a field.}
  \item{start}{A numeric vector with starting values.}
  \item{threshold}{Numeric; a value indicating a threshold for the
    binary random field.}
  \item{type}{String; the type of estimation method.}
  \item{varest}{Logical; if \code{TRUE} the estimates' variances and
    standard errors are returned.
    \code{FALSE} is the default.}
  \item{vartype}{String; the type of estimation method for computing
    the estimate variances, see the Section \bold{Details}.}
  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted, see \code{\link{FitComposite}}.}
  \item{winconst}{Numeric; a positive real value indicating the window size used from the
   sub-sampling method for the estimation of the parameters'
   variances.}.
 \item{winstp}{Numeric; a value in \eqn{(0,1]} for defining the window
    step. See \code{\link{FitComposite}}.}
}

\seealso{\code{\link{FitComposite}}, \code{\link{WLeastSquare}}.}

\author{Simone Padoan, \email{simone.padoan@stat.unipd.it},
  \url{http://homes.stat.unipd.it/padoan}.}

\keyword{LeastSquare}
