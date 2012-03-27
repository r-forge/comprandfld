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
Covariogram(fitted, lags=NULL, lagt=NULL, answer.cov=FALSE, answer.vario=FALSE,
            answer.extc=FALSE, answer.range=FALSE, fix.lags=NULL, fix.lagt=NULL,
            show.cov=FALSE, show.vario=FALSE, show.extc=FALSE, show.range=FALSE,
            add.cov=FALSE, add.vario=FALSE, add.extc=FALSE, pract.range=95,
            vario, \dots)
}
\arguments{
  \item{fitted}{A fitted object obtained from the
    \code{\link{FitComposite}} or \code{\link{WLeastSquare}} procedures.}
  \item{lags}{A numeric vector of distances.}
  \item{lagt}{A numeric vector of temporal separations.}
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
  \item{fix.lags}{Integer; a positive value denoting the spatial lag to
    consider for the plot of the temporal profile.}
  \item{fix.lagt}{Integer; a positive value denoting the temporal lag to
    consider for the plot of the spatial profile.}
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

\author{Simone Padoan, \email{simone.padoan@stat.unipd.it},
  \url{http://homes.stat.unipd.it/padoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@unibg.it},
   \url{www.unibg.it/pers/?moreno.bevilacqua}}

\examples{
library(CompRandFld)
library(RandomFields)
library(scatterplot3d)
set.seed(31231)

# Set the coordinates of the points:
x <- runif(100, 0, 10)
y <- runif(100, 0, 10)


################################################################
###
### Example 1. Plot of covariance and variogram functions
### estimated from a Gaussian random field with exponent 
### correlation. One spatial replication is simulated.
###
###
###############################################################

# Set the model's parameters:
corrmodel <- "exponential"
mean <- 0
sill <- 1
nugget <- 0
scale <- 3

# Simulation of the Gaussian random field:
data <- RFsim(x, y, corrmodel=corrmodel, param=list(mean=mean,
              sill=sill, nugget=nugget, scale=scale))$data

# Maximum composite-likelihood fitting of the Gaussian random field:
fit <- FitComposite(data, x, y, corrmodel=corrmodel,
                    type='Difference')

# Empirical estimation of the madogram:
vario <- EVariogram(data, x, y)

# Plot of covariance and variogram functions:
par(mfrow=c(1,2))
Covariogram(fit, show.cov=TRUE, show.range=TRUE,
            show.vario=TRUE, vario=vario)


################################################################
##
### Example 2. Plot of covariance and extremal coefficient
### functions estimated from a max-stable random field with
### exponential correlation. n idd spatial replications are
### simulated.
###
###############################################################

set.seed(1126)
# Simulation of the max-stable random field:
data <- RFsim(x, y, corrmodel=corrmodel, model="ExtGauss", replicates=20,
              param=list(mean=mean,sill=sill,nugget=nugget,scale=scale))$data

# Maximum composite-likelihood fitting of the max-stable random field:
fit <- FitComposite(data, x, y, corrmodel=corrmodel, model='ExtGauss',
                    replicates=20, varest=TRUE, vartype='Sampling',
                    margins="Frechet")

data <- Dist2Dist(data, to='sGumbel')

# Empirical estimation of the madogram:
vario <- EVariogram(data, x, y, type='madogram', replicates=20)

# Plot of correlation and extremal coefficient functions:
par(mfrow=c(1,2))
Covariogram(fit, show.cov=TRUE, show.range=TRUE, show.extc=TRUE,
            vario=vario, pract.range=84)


################################################################
###
### Example 3. Plot of covariance and variogram functions
### estimated from a Gaussian spatio-temporal random field with
### double-exp correlation.
### One spatio-temporal replication is simulated.
###
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(20, 0, 1)
y <- runif(20, 0, 1)
# Define the temporal sequence:
time <- seq(0, 50, 1)

# Simulation of the spatio-temporal Gaussian random field:
data <- RFsim(x, y, time, corrmodel="exp_exp",param=list(mean=mean,
              nugget=nugget,scale_s=0.5,scale_t=1,sill=sill))$data

# Maximum composite-likelihood fitting of the space-time Gaussian random field:
fit <- FitComposite(data, x, y, time, corrmodel="exp_exp", maxtime=5,
                    likelihood="Marginal",type="Pairwise", fixed=list(
                    nugget=nugget, mean=mean),start=list(scale_s=0.2,
                    scale_t=1, sill=sill))

# Empirical estimation of spatio-temporal covariance:
vario <- EVariogram(data, x, y, time, maxtime=10)

# Plot of the fitted space-time covariace
Covariogram(fit,show.cov=TRUE)

# Plot of the fitted space-time variogram
Covariogram(fit,vario=vario,show.vario=TRUE)

# Plot of covariance, variogram and spatio and temporal profiles:
Covariogram(fit,vario=vario,fix.lagt=1,fix.lags=1,show.vario=TRUE)

################################################################
###
### Example 4. Plot of parametric and empirical lorelograms
### estimated from a Binary Gaussian random fields with 
### exponential correlation. One spatial replication is
### simulated.
###
###############################################################

set.seed(1240)

# Define the spatial-coordinates of the points:
x <- seq(0,3, 0.1)
y <- seq(0,3, 0.1)

# Simulation of the Binary Gaussian random field:
data <- RFsim(x, y, corrmodel=corrmodel, model="BinaryGauss",
              threshold=0,param=list(nugget=nugget,mean=mean,
              scale=.6,sill=0.8))$data

# Maximum composite-likelihood fitting of the Binary Gaussian random field:
fit <- FitComposite(data, x, y, corrmodel=corrmodel, model="BinaryGauss",
                    maxdist=0.8, likelihood="Marginal", type="Pairwise",
                    start=list(mean=mean,scale=0.1,sill=0.1))

# Empirical estimation of the lorelogram:
vario <- EVariogram(data, x, y, type="lorelogram", maxdist=2)

# Plot of fitted and empirical lorelograms:
Covariogram(fit, vario=vario, show.vario=TRUE, lags=seq(0.1,2,0.1))
}

\keyword{Composite}
