\name{EVariogram}
\alias{EVariogram}

\title{Empirical Variogram of Gaussian and Max-Stable Random Fields}
\description{
  The function returns empirical estimates of the variogram of a Gaussian
  or a max-stable random field.
}
\usage{
EVariogram(data, coordx, coordy, coordt=NULL, cloud=FALSE, grid=FALSE,
           gev=c(0,1,0), lonlat=FALSE, maxdist=NULL, maxtime=NULL,
           numbins=NULL, replicates=1, type='variogram')
}
\arguments{
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{n \times d}{n x d})-matrix
   (\eqn{n} iid spatial realisations) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid)
   or an (\eqn{d \times d \times n}{d x d x n})-array (\eqn{n} iid spatial realisations on regular grid) or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation) or an (\eqn{t \times d \times n }{t x d x n})-array
   (\eqn{n} iid spatial-temporal realisations) or or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid) or an (\eqn{d \times d \times t \times n }{d x d x t x n})-array
   (\eqn{n} iid spatial-temporal realisations on regular grid). See \code{\link{FitComposite}} for details.}
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of spatial sites) assigning 2-dimensions of spatial coordinates or a numeric \eqn{d}{d}-dimensional vector assigning
    1-dimension of spatial coordinates.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    spatial coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, the default is \code{NULL} then \code{coordx} is expected to
    be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
  \item{coordt}{A numeric vector assigning 1-dimension of
    temporal coordinates. Optional argument, the default is \code{NULL} then a spatial random field is expected.}
  \item{cloud}{Logical; if \code{TRUE} the variogram cloud is computed,
  otherwise if \code{FALSE} (the default) the empirical (binned)
  variogram is returned.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set of non-equispaced spatial sites.}
  \item{gev}{A numeric vector with the three GEV parameters;}
  \item{lonlat}{Logical; if \code{FALSE} (the default), \code{coordx}
    and \code{coordy} are interpreted as Cartesian coordinates
    otherwise they are considered as longitude and latitude.}
  \item{maxdist}{A numeric value denoting the spatial maximum distance,
   see the Section \bold{Details}.}
   \item{maxtime}{A numeric value denoting the temporal maximum distance,
   see the Section \bold{Details}.}
  \item{numbins}{A numeric value denoting the numbers of bins, see the
    Section \bold{Details}.}
  \item{replicates}{Numeric; a positive integer denoting the number of independent and identically distributed (iid)
    replications of a spatial or spatial-temporal random field. Optional argument, the default value is \eqn{1} then
    a single realisation is considered.}
  \item{type}{A String denoting the type of variogram. Three options
    are available, \code{variogram} which return the standard variogram
    and for extreme values the \code{madogram} or \code{Fmadogram}.}
}

\details{

  The procedure, for a spatial random field, uses the following
  definition of the sample variogram estimator:
  \deqn{\hat{\gamma}(h) = 0.5 \sum_{x_i, x_j \in N(h)} (Z(x_i) - Z(x_j))^2 / |N(h)|}
  where \eqn{N(h)} is the set of all the sample pairs whose distances fall into a tolerance region with size \eqn{h}
  (equispaced intervalls are considered).
  Observe, that in the literature often the above definition is termed semivariogram (see i.e. the first reference).
  Nevertheless, here this defition has been used in order to be consistent with the variogram defition used for the extremes
  (see i.e. the third reference).

  The procedure, for a max-stable random field, uses the following defition of the sample madogram estimator:
  \deqn{\hat{\nu}(h) = 0.5 \sum_{x_i, x_j \in N(h)} |Z(x_i) - Z(x_j)| / |N(h)|.}

  The procedure, for a spatial-temporal random field, uses the following
  definition of the sample variogram estimator:
  \deqn{\hat{\gamma}(h, u) = 0.5 \sum_{(x_i, l), (x_j, k) \in N(h, u)} (Z(x_i, l) - Z(x_j, k))^2 / |N(h, u)|}
  where \eqn{N(h, u)} is the set of all the sample pairs whose spatial distances fall into a tolerance region with size \eqn{h}
  and \eqn{|k - l| = u}. Note, that \eqn{Z(x_i, l)} is the observation at site \eqn{x_i} and time \eqn{l}.

  The  \code{numbins} parameter indicates the number of adjacent
  intervals to consider in order to grouped distances with which to
  compute the (weighted) lest squares.

  The \code{maxdist} parameter indicates the maximum spatial distance below which
  the shorter distances will be considered in the calculation of
  the (weigthed) least squares.

  The \code{maxtime} parameter indicates the maximum temporal distance below which
  the shorter distances will be considered in the calculation of
  the (weigthed) least squares.
}

\value{
  Returns an object of class \code{Variogram}.
  An object of class \code{Variogram} is a list containing
  at most the following components:

  \item{bins}{Adjacent intervals of grouped spatial distances if
  \code{cloud=FALSE}. Otherwise if \code{cloud=TRUE} all the spatial pairwise distances;}
  \item{bint}{Adjacent intervals of grouped temporal distances if
  \code{cloud=FALSE}. Otherwise if \code{cloud=TRUE} all the temporal pairwise distances;}
  \item{cloud}{If the variogram cloud is returned (\code{TRUE}) or the
  empirical variogram (\code{FALSE});}
  \item{centers}{The centers of the spatial bins;}
  \item{extcoeff}{The spatial extremal coefficient function. Available only if
  \code{type} is equal to \code{madogram} or \code{Fmadogram} (for the moment available only for a spatial random field);}
  \item{lenbins}{The number of pairs in each spatial bin;}
  \item{lenbinst}{The number of pairs in each spatial-temporal bin;}
  \item{lenbint}{The number of pairs in each temporal bin;}
  \item{srange}{The maximum and minimum spatial distances used for the calculation of the variogram;}
  \item{variograms}{The empirical spatial variogram;}
  \item{variogramst}{The empirical spatial-temporal variogram;}
  \item{variogramt}{The empirical temporal variogram;}
  \item{trange}{The maximum and minimum temporal distance used for the calculation of the variogram;}
  \item{type}{The type of estimated variogram: the standard variogram or
    the madogram.}
}

\references{

  Cressie, N. A. C. (1993) \emph{Statistics for Spatial Data}. New York: Wiley.

  Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modelling}.
  Spring Verlang, New York.

  Cooley, D., Naveau, P. and Poncet, P. (2006) \emph{Variograms for spatial
  max-stable random fields}. Dependence in Probability and
  Statistics, p. 373--390.
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibg.it},
  \url{http://people.epfl.ch/simone.padoan}.}

\examples{
library(RandomFields)
set.seed(2111)

# Set the coordinates of the sites:
x <- runif(100, 0, 20)
y <- runif(100, 0, 20)


################################################################
###
### Example 1. Empirical variogram of Gaussian random fields
### with stable correlation function.
### One single replication is considered.
###
###############################################################

# Set the model's parameters:
corrmodel <- "stable"
mean <- 0
variance <- 1
nugget <- 0
scale <- 10
power <- 1.2

# Simulation of a Gaussian random field in the specified points:
data <- GaussRF(x=x, y=y, model=corrmodel, grid=FALSE,
                param=c(mean, variance, nugget, scale, power))

# Empirical variogram estimation:
fit <- EVariogram(data, x, y)

# Results:
plot(fit$centers, fit$variograms, xlab='h', ylab=expression(gamma(h)),
     ylim=c(0, max(fit$variograms)), xlim=c(0, fit$srange[2]), pch=20)


################################################################
###
### Example 2. Empirical mandogram of max-stable random fields.
### Extremal Gaussian model with stable correlation function.
### n iid spatial replications are considered.
###
###############################################################

# Simulation of a max-stable random field in the specified points:
data <- MaxStableRF(x=x, y=y, model=corrmodel, grid=FALSE, maxstable="extr",
                    param=c(mean, variance, nugget, scale, power), n=40, pch='')
data <- t(data)
data <- Dist2Dist(data, to='sGumbel')

# Empirical madogram estimation:
fit <- EVariogram(data, x, y, type='madogram', replicates=40, cloud=TRUE)

# Results:
par(mfrow=c(1,2))
plot(fit$centers, fit$variograms, xlab='h', ylab=expression(nu(h)),
     ylim=c(0, max(fit$variograms)), xlim=c(0, fit$srange[2]), pch=20)
plot(fit$centers, fit$extcoeff, xlab='h', ylab=expression(theta(h)),
     ylim=c(1, 2), xlim=c(0, fit$srange[2]), pch=20)

}

\keyword{Variogram}
