\name{EVariogram}
\alias{EVariogram}

\title{Empirical Estimation of the Variogram for Gaussian and Max-Stable
Random Fields}
\description{
  the function returns an empirical estimate of the variogram for a Gaussian
  or a max-stable random field.
}
\usage{
EVariogram(coordx, coordy, data, cloud=FALSE, extcoeff=FALSE, grid=FALSE,
           gev=c(0,1,0), lonlat=FALSE, maxdist=NULL, numbins=NULL,
           replicates=FALSE, type='variogram')
} 
\arguments{
 \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of points) assigning 2-dimensions of coordinates or a numeric vector assigning
    1-dimension of coordinates.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector otherwise it will be ignored.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations (see
    \code{\link{FitComposite}}).}
  \item{cloud}{Logical; if \code{TRUE} the variogram cloud is computed,
  otherwise if \code{FALSE} (the default) the empirical (binned)
  variogram is returned.}
  \item{extcoeff}{Logical; if the extremal coefficient estimates need to
  be returned. Available only if \code{type} is equal to \code{madogram}
  or \code{Fmadogram}.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{gev}{A numeric vector with the three GEV parameters;}
  \item{lonlat}{Logical; if \code{FALSE} (the default), \code{coordx}
    and \code{coordy} are interpreted as Cartesian coordinates
    otherwise they are considered as longitude and latitude.}
  \item{maxdist}{A numeric value denoting the maximum distance, see the Section \bold{Details}.}
  \item{numbins}{A numeric value denoting the numbers of bins, see the
    Section \bold{Details}.}
  \item{replicates}{Logical; if \code{FALSE} (the default) one spatial random field is
    considered, instead if \code{TRUE} the data are considered as iid replicates of a field.}
  \item{type}{A String denoting the type of variogram. Three options
    are available, \code{variogram} which return the standard variogram
    and for extreme values the \code{madogram} or \code{Fmadogram}.}
}

\details{
  The  \code{numbins} parameter indicates the number of adjacent
  intervals to consider in order to grouped distances with which to
  compute the (weighted) lest squares.

  The \code{maxdist} parameter indicates the maximum distance below which
  the shorter distances will be considered in the calculation of
  the (weigthed) least squares.
}

\value{
  Returns an object of class \code{Variogram}.
  An object of class \code{Variogram} is a list containing
  at most the following components:

  \item{bins}{Adjacent intervals of grouped distances if
  \code{cloud=FALSE}. Otherwise the pairwise distances if
  \code{cloud=TRUE};}
  \item{cloud}{If the variogram cloud is returned (\code{TRUE}) or the
  empirical variogram (\code{FALSE});}
  \item{centers}{The centers of the bins;}
  \item{extremalcoeff}{The extremal coefficient estimates;}
  \item{lenbins}{The number of pairs in each bin;}
  \item{maxdist}{The maximum distance used for the calculation of the variogram;}
  \item{variogram}{The empirical variogram;}
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

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://people.epfl.ch/simone.padoan}.}

\examples{
library(RandomFields)
set.seed(2111)

# Set the coordinates of the sites:
x <- runif(100, 0, 20) 
y <- runif(100, 0, 20)


################################################################
###
### Example 1. Empirical variogram estimate of one
### spatial realisation of a Gaussian random field.
### 
###
###############################################################

# Set the model's parameters:
corrmodel <- "stable"   
mean <- 0
variance <- 1
nugget <- 0
scale <- 10
power <- 1.5

# Simulation of a Gaussian random field in the specified points:
sim <- GaussRF(x=x, y=y, model=corrmodel, grid=FALSE,
              param=c(mean, variance, nugget, scale, power))

# Empirical variogram estimation:
fit <- EVariogram(x, y, sim)

# Results:
plot(fit$centers, fit$variogram, xlab='h', ylab=expression(gamma(h)),
     ylim=c(0, max(fit$variogram)), xlim=c(0, fit$maxdist), pch=20)


################################################################
###
### Example 2. Empirical mandogram estimate of one
### spatial realisation of a max-stable random field.
### 
###
###############################################################

# Simulation of a max-stable random field in the specified points:
sim <- MaxStableRF(x=x, y=y, model=corrmodel, grid=FALSE, maxstable="extr",
                   param=c(mean, variance, nugget, scale, power), n=100)
sim <- t(sim)
sim <- Dist2Dist(sim, to='sGumbel')

# Empirical madogram estimation:
fit <- EVariogram(x, y, sim, type='madogram', replicates=TRUE,
                  extcoeff=TRUE, cloud=TRUE)

# Results:
par(mfrow=c(1,2))
plot(fit$centers, fit$variogram, xlab='h', ylab=expression(nu(h)),
     ylim=c(0, max(fit$variogram)), xlim=c(0, fit$maxdist), pch=20)
plot(fit$centers, fit$extremalcoeff, xlab='h', ylab=expression(theta(h)),
     ylim=c(1, 2), xlim=c(0, fit$maxdist), pch=20)

}

\keyword{Variogram}
