\name{WLeastSquare}
\alias{WLeastSquare}
\alias{print.WLS}

\title{Weighted Least Square Estimation of Random Fields}
\description{
  the function returns the parameters' estimates and the estimates' variances
  of a random field obtained by the weigthed least squares estimator.
}
\usage{
WLeastSquare(coordx, coordy, corrmodel, data, fixed=NULL,
             grid=FALSE, lonlat=FALSE, maxdist=NULL,
             optimizer='Nelder-Mead', numbins=NULL, start=NULL,
             time=FALSE, weighted=FALSE)
} 
\arguments{
 \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of points) assigning 2-dimensions of coordinates or a numeric vector assigning
    1-dimension of coordinates.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector otherwise it will be ignored.}
  \item{corrmodel}{String; the name of a correlation model, for the
    description (see \code{\link{FitComposite}}).}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations (see \code{\link{FitComposite}}).}
  \item{fixed}{A named list giving the values of the parameters that
    will be considered as known values. The listed parameters for a
    given correlation function will be not estimated, i.e. if
    \code{list(nugget=0)} the nugget effect is ignored.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{lonlat}{Logical; if \code{FALSE} (the default), \code{coordx}
    and \code{coordy} are interpreted as Cartesian coordinates
    otherwise they are considered as longitude and latitude.}
  \item{maxdist}{A numeric value denoting the maximum distance, see the Section \bold{Details}.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). 'Nelder-Mead' is the
    default.}
  \item{numbins}{A numeric value denoting the numbers of bins, see the Section \bold{Details}}
  \item{start}{A named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure. \code{NULL} is the default (see \code{\link{FitComposite}}).}
  \item{time}{Logical; if \code{FALSE} (the default) a spatial random field is
    considered (one temporal realisation), if \code{TRUE} a
    spatial-temporal random field is considered (see \code{\link{FitComposite}}).}
  \item{weighted}{Logical; if \code{TRUE} then the weighted least square
    estimator is considered. If \code{FALSE} (the
    default) then the classic least square is used.}
}

\details{
  Insert description of \code{maxdist} and \code{numbins}.
}

\value{
  The returned object is a list with:
  \itemize{
    \item \code{bins}, adjacent intervals of grouped distances;
    \item \code{coord}, the vector of coordinates;
    \item \code{convergence}, a string that denotes if convergence is
    reached;
    \item \code{corrmodel}, the correlation model;
    \item \code{data}, the vector or matrix of data;
    \item \code{fixed}, the vector of fixed parameters;
    \item \code{iterations}, the number of iteration used by the
    numerical routine;
    \item \code{message}, extra message passed from the numerical
    routines;
    \item \code{param}, the vector of  parameters' estimates;
    \item \code{variogram} the empirical variogram;
    }
}

\references{
  Cressie, N. A. C. (1993) Statistics for Spatial Data. New York: Wiley.

  Barry, J. T., Crowder, M. J. and Diggle, P. J. (1997)
  Parametric estimation of the variogram. Tech. Report, Dept Maths & Stats, Lancaster University.
}

\seealso{\code{\link{FitComposite}}, \code{\link{optim}}}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php};
  Moreno Bevilacqua, \email{moreno.bevilacqua@unive.it}.}

\examples{
library(RandomFields)
set.seed(2111)

# Set the coordinates of the sites:
x <- runif(100, 0, 20) 
y <- runif(100, 0, 20)


################################################################
###
### Example 1. Least square fitting of one
### spatial realisation of a Gaussian random field.
### Non weighted version (all weights equals to 1)
###
###############################################################

# Set the model's parameters:
corrmodel <- "stable"   
mean <- 0
variance <- 1
nugget <- 0
scale <- 10
power <- 1.5

# Simulation of the Gaussian random field in the specified points:
sim <- GaussRF(x=x, y=y, model=corrmodel, grid=FALSE,
              param=c(mean, variance, nugget, scale, power))

# Least square fitting of the random field:
fit <- WLeastSquare(x, y, corrmodel, sim)

# Results:
print(fit)


################################################################
###
### Example 1. Weighted least square fitting of one
### spatial realisation of a Gaussian random field.
### Weighted version.
###
###############################################################

# Set the model's parameters:
corrmodel <- "stable"   
mean <- 0
variance <- 1
nugget <- 0
scale <- 10
power <- 1.5

# Simulation of the Gaussian random field in the specified points:
sim <- GaussRF(x=x, y=y, model=corrmodel, grid=FALSE,
              param=c(mean, variance, nugget, scale, power))

# Least square fitting of the random field:
fit <- WLeastSquare(x, y, corrmodel, sim, weighted=TRUE)

# Results:
print(fit)


}

\keyword{LeastSquare}
