\name{InitParam}
\alias{InitParam}

\title{Initializes the Parameters for Estimation Procedures}
\description{
  Subroutine called by the fitting procedures.
  The procedure initializes the parameters for the fitting procedure.
}
\usage{
InitParam(coordx, coordy, coordt, corrmodel, data, distance, fcall,
          fixed, grid,likelihood, margins, maxdist, maxtime, model,
          numblock, param, parscale, paramrange, replicates, start,
          taper, tapsep, threshold, type, typereal, varest, vartype,
          weighted, winconst, winstp)
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
  \item{distance}{String; the name of the spatial distance. The default is \code{Eucl},
            the euclidean distance. See the Section \bold{Details}.}
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
  \item{margins}{String; the type of the marginal distribution of the
    max-stable field.}
  \item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial distance considered in the composite-likelihood
    computation.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal lag considered in the composite-likelihood
    computation.}
  \item{model}{String; the density associated to the likelihood
    objects. \code{Gaussian} is the default.}
  \item{numblock}{Numeric; the observation size of the underlying random
    field. Only in case of max-stable random fields and in the simulation.}
  \item{param}{A numeric vector of parameter values required in the
    simulation procedure of random fields.}
  \item{parscale}{A numeric vector of scaling factor to improve the
    maximizing procedure, see \code{\link{optim}}.}
  \item{paramrange}{A numeric vector of parameters ranges, see
    \code{\link{optim}}.}
  \item{replicates}{Logical; if \code{FALSE} (the default) one spatial random field is
    considered, instead if \code{TRUE} the data are considered as iid replicates of a field.}
   \item{start}{A named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure.}
  \item{taper}{String; the name of the type of covariance matrix.
     It can be   \code{Standard} (the default value) or
    \code{Tapering} for taperd covariance matrix.}
   \item{tapsep}{Numeric; an optional  value indicating the separabe parameter in the space time quasi taper (see \bold{Details}).}
   \item{threshold}{Numeric; a value indicating a threshold for the
    binary random field.}
  \item{type}{String; the type of likelihood objects. Temporary value set
    to be "WLeastSquare" (weigthed least-square) in order to compute the starting values.}
  \item{typereal}{String; the real type of likelihood objects. See \code{\link{FitComposite}}.}
  \item{varest}{Logical; if \code{TRUE} the estimates' variances and
    standard errors are returned.
    \code{FALSE} is the default.}
  \item{vartype}{String; the type of estimation method for computing
    the estimate variances, see the Section \bold{Details}.}
  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted, see \code{\link{FitComposite}}.}
  \item{winconst}{Numeric; a positive real value indicating the
    sub-window in the sub-sampling procedure. See \code{\link{FitComposite}}.}
  \item{winstp}{Numeric; a value in \eqn{(0,1]} for defining the window
    step. See \code{\link{FitComposite}}.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home}.}

\section{Acknowledgements}{
We thank the CRAN team for having helped us in finding a bug in the code. 
}

\keyword{Composite}