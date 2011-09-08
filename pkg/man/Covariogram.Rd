\name{Covariogram}
\alias{Covariogram}

\title{Computes covariance, variogram and extremal coefficient functions}
\description{
  The procedure computes and/or plots the covariance, the variogram
  or the extremal coefficient functions and the practical range estimated fitting a Gaussian or
  max-stable random field
  with the composite-likelihood or using the weighted least square
  method. Allows to add to the variogram or extremal coefficient plots
  the empirical estimates.
}
\usage{
Covariogram(fitted, lags=NULL, answer.cov=FALSE, answer.vario=FALSE,
            answer.extc=FALSE, answer.range=FALSE, show.cov=FALSE,
            show.vario=FALSE, show.extc=FALSE, show.range=FALSE,
            add.cov=FALSE, add.vario=FALSE, add.extc=FALSE,
            pract.range=95, vario, \dots)
}
\arguments{
  \item{fitted}{A fitted object obtained from the
    \code{\link{FitComposite}} or \code{\link{WLeastSquare}} procedures.}
  \item{lags}{A numeric vector of distances.}
  \item{answer.cov}{Logical; if \code{TRUE} a vector with the
    estimated covariance function is returned; if \code{FALSE} (the
    default) the covariance is not returned.}
  \item{answer.vario}{Logical; if \code{TRUE} a vector with the
    estimated variogram is returned; if \code{FALSE} (the
    default) the variogram is not returned.}
  \item{answer.extc}{Logical; if \code{TRUE} a vector with the
    estimated extremal coefficient is returned; if \code{FALSE} (the
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
  \item{show.extc}{Logical; if \code{TRUE} the estimated
    extremal coefficient is plotted; if \code{FALSE} (the
    default) the extremal coefficient is not plotted.}
  \item{show.range}{Logical; if \code{TRUE} the estimated
    pratical range is added on the plot; if \code{FALSE} (the
    default) the pratical range is not added.}
  \item{add.cov}{Logical; if \code{TRUE} the vector of the
    estimated covariance function is added on the current plot; if \code{FALSE} (the
    default) the covariance is not added.}
  \item{add.vario}{Logical; if \code{TRUE} the vector with the
    estimated variogram is added on the current plot; if \code{FALSE} (the
    default) the correlation is not added.}
  \item{add.extc}{Logical; if \code{TRUE} the vector with the
    estimated extremal coefficient is added on the current plot; if \code{FALSE} (the
    default) the correlation is not added.}
  \item{pract.range}{Numeric; the percent of the sill to be reached.}
  \item{vario}{A \code{Variogram} object obtained from the
    \code{\link{EVariogram}} procedure.}
  \item{\dots}{other optional parameters which are passed to plot functions.}
}

\value{
  The returned object is eventually a list with:

    \item{covariance}{The vector of the estimated covariance function;}
    \item{variogram}{The vector of the estimated variogram function;}
    \item{extrcoeff}{The vector of the estimated extremal coefficient function;}
    \item{pratical.range}{The estimated practial range.}
}

\references{
  Gaetan, C. and Guyon, X. (2010) Spatial Statistics and Modelling.
  \emph{Spring Verlang, New York}.
}

\seealso{\code{\link{FitComposite}}, \code{\link{WLeastSquare}}}

\author{Simone Padoan, \email{simone.padoan@unibg.it},
  \url{http://www.unibg.it/pers/?simone.padoan}.}

\examples{
library(RandomFields)
set.seed(32378)

# Set the coordinates of the points:
x <- runif(100, 0, 20)
y <- runif(100, 0, 20)


################################################################
###
### Example 1. Plot of the covariance and variogram functions.
### One spatial realisation of a Gaussian random field.
###
###
###############################################################

# Set the model's parameters:
corrmodel <- "stable"
mean <- 0
variance <- 1
nugget <- 0
scale <- 10
power <- 0.8

# Simulation of the Gaussian random field in the specified points:
data <- GaussRF(x=x, y=y, model=corrmodel, grid=FALSE,
                param=c(mean, variance, nugget, scale, power))

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(data, x, y, corrmodel=corrmodel,
                    type='Difference')
vario <- EVariogram(data, x, y)
# Plot of the covariance and variogram functions:
par(mfrow=c(1,2))
Covariogram(fit, show.cov=TRUE, show.range=TRUE,
            show.vario=TRUE, vario=vario)


################################################################
##
### Example 2. Plot of the covariance and extremal coefficient
### functions.
### n idd spatial realisations of a max-stable random field.
###
###
###############################################################

set.seed(179922)
# Simulation of a max-stable random field in the specified points:
data <- MaxStableRF(x=x, y=y, model=corrmodel, grid=FALSE, maxstable="extr",
                    param=c(mean, variance, nugget, scale, power), n=20, pch='')
data <- t(data)

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(data, x, y, corrmodel=corrmodel, model='ExtGauss',
                    replicates=20, varest=TRUE, vartype='Sampling',
                    margins="Frechet")

data <- Dist2Dist(data, to='sGumbel')

# Empirical madogram estimation:
vario <- EVariogram(data, x, y, type='madogram', replicates=20)

# Plot of the correlation and extremal coefficient functions:
par(mfrow=c(1,2))
Covariogram(fit, show.cov=TRUE, show.range=TRUE,
            show.extc=TRUE, vario=vario, pract.range=84)

}

\keyword{Composite}
