\name{Covariogram}
\alias{Covariogram}

\title{Computes (and Plots) the covariance function and the variogram}
\description{
  The procedure computes and plots the estimated covariance function
  and the variogram from a fitted model obtained fitting a random field
  with the composite-likelihood or using the weighted least square method.
}
\usage{
Covariogram(fitted, lags=NULL, answer.cov=FALSE, answer.vario=FALSE,
            answer.range=FALSE, show.cov=FALSE, show.vario=FALSE,
            show.range=FALSE, add.cov=FALSE, add.vario=FALSE,
            pract.range=95, vario, \dots)
} 
\arguments{
  \item{fitted}{The fitted object obtained from the
    \code{\link{FitComposite}} procedure.}
  \item{lags}{A numeric vector of distances.}
  \item{answer.cov}{Logical; if \code{TRUE} a vector with the
    estimated covariance function is returned; if \code{FALSE} (the
    default) the covariance is not returned.}
  \item{answer.vario}{Logical; if \code{TRUE} a vector with the
    estimated variogram is returned; if \code{FALSE} (the
    default) the variogram is not returned.}
  \item{answer.range}{Logical; if \code{TRUE} the estimated
    pratical range is returned; if \code{FALSE} (the
    default) the pratical range is not returned.}
  \item{show.cov}{Logical; if \code{TRUE} the estimated
    covariance function is plotted; if \code{FALSE} (the
    default) the covariance function is not plotted.}
  \item{show.vario}{Logical; if \code{TRUE} the estimated
    variogram is plotted; if \code{FALSE} (the
    default) the variogram is not plotted.}
  \item{show.range}{Logical; if \code{TRUE} the estimated
    pratical range is added on the plot; if \code{FALSE} (the
    default) the pratical range is not added.}
  \item{add.cov}{Logical; if \code{TRUE} the vector of the
    estimated covariance function is added on the current plot; if \code{FALSE} (the
    default) the covariance is not added.}
  \item{add.vario}{Logical; if \code{TRUE} the vector with the
    estimated variogram is added on the current plot; if \code{FALSE} (the
    default) the correlation is not added.}
  \item{pract.range}{Numeric; the percent of the sill to be reached.}
  \item{vario}{An object of type \code{Variogram} obtained from \code{\link{EVariogram}}.}
  \item{\dots}{other optional parameters which are passed to plot function.}
}

\value{
  The returned object is a list with:
  
    \item{covariance}{The vector of the estimated covariance function;}
    \item{variogram}{The vector of the estimated variogram function;}
    \item{pratical.range}{The estimated practial range.}
}

\references{
  Gaetan, C. and Guyon, X. (2010) Spatial Statistics and Modelling.
  \emph{Spring Verlang, New York}.
}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://people.epfl.ch/simone.padoan}.}

\examples{
library(RandomFields)
set.seed(2111)

# Set the coordinates of the points:
x <- runif(100, 0, 20) 
y <- runif(100, 0, 20)


################################################################
###
### Example 1. Plot of the estimated correlation function
### from a spatial realisation of a Gaussian random field.
###
###
###############################################################

# Set the model's parameters:
corrmodel <- "stable"   
mean <- 0
variance <- 3
nugget <- 1
scale <- 10
power <- 1.5

# Simulation of the Gaussian random field in the specified points:
sim <- GaussRF(x=x, y=y, model=corrmodel, grid=FALSE,
               param=c(mean, variance, nugget, scale, power))

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(x, y, corrmodel, sim)

# Plot of the Correlation function:
par(mfrow=c(1,2))
Covariogram(fit, show.cov=TRUE, show.range=TRUE, show.vario=TRUE)

}

\keyword{Composite}
