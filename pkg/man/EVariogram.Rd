\name{EVariogram}
\alias{EVariogram}

\title{Empirical Variogram(variants) of Gaussian, Binary and Max-Stable Fields}
\description{
  The function returns an empirical estimate of the variogram (or its variants) for Gaussian,
  Binary and max-stable random field.
}
\usage{
EVariogram(data, coordx, coordy, coordt=NULL, cloud=FALSE, distance='Eucl',
           grid=FALSE, gev=c(0,1,0), maxdist=NULL, maxtime=NULL, numbins=NULL,
           replicates=1, type='variogram')
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
  \item{distance}{String; the name of the spatial distance. The default
    is \code{Eucl}, the euclidean distance. See the Section
    \bold{Details} of \code{\link{FitComposite}}.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set of non-equispaced spatial sites.}
  \item{gev}{A numeric vector with the three GEV parameters;}
  \item{maxdist}{A numeric value denoting the spatial maximum distance,
   see the Section \bold{Details}.}
   \item{maxtime}{A numeric value denoting the temporal maximum distance,
   see the Section \bold{Details}.}
  \item{numbins}{A numeric value denoting the numbers of bins, see the
    Section \bold{Details}.}
  \item{replicates}{Numeric; a positive integer denoting the number of independent and identically distributed (iid)
    replications of a spatial or spatial-temporal random field. Optional argument, the default value is \eqn{1} then
    a single realisation is considered.}
  \item{type}{A String denoting the type of variogram. Four options
    are available: \code{variogram}, \code{madogram}, \code{Fmadogram}
    and \code{lorelogram}. It is returned respectively,
    the standard variogram with the first (Gaussian responses),
    the madogram with the second and third (extreme values), the
    lorelogram with the fourth (Binary data).}
}

\details{

  We briefly report the definitions of variogram used in this function.

  In the case of a spatial Gaussian random field
  the sample \code{variogram} estimator is defined by
  \deqn{\hat{\gamma}(h) = 0.5 \sum_{x_i, x_j \in N(h)} (Z(x_i) - Z(x_j))^2 / |N(h)|}
  where \eqn{N(h)} is the set of all the sample pairs whose distances fall into a tolerance region with size \eqn{h}
  (equispaced intervalls are considered).
  Observe, that in the literature often the above definition is termed semivariogram (see e.g. the first reference).
  Nevertheless, here this defition has been used in order to be consistent with the variogram defition used for the extremes
  (see e.g. the third reference).

  In the case of a spatial max-stable random field, the sample \code{madogram}
  estimator is defined similarly to the Gaussian case by
  \deqn{\hat{\nu}(h) = 0.5 \sum_{x_i, x_j \in N(h)} |Z(x_i) - Z(x_j)| /
    |N(h)|.}

  In the case of a spatial binary random field, the sample \code{lorelogram}
  estimator (the analogue of the correlation) is defined by
  \deqn{\hat{L}(h) = (N_{11}(h) N_{00}(h) )/ (N_{01}(h) N_{10}(h)).}
  where \eqn{N_{11}(h)} is the number of pairs who are both equal to
  \eqn{1}{1} and that falls in the bin
  \eqn{h}{h}. Similarly are defined the other quantities.

  In the case of a spatio-temporal Gaussian random field the sample
  \code{variogram} estimator is defined by
  \deqn{\hat{\gamma}(h, u) = 0.5 \sum_{(x_i, l), (x_j, k) \in N(h, u)} (Z(x_i, l) - Z(x_j, k))^2 / |N(h, u)|}
  where \eqn{N(h, u)}{N(h, u)} is the set of all the sample pairs whose
  spatial distances fall into a tolerance region with size \eqn{h}{h}
  and \eqn{|k - l| = u}{\|k-l\|=u}. Note, that \eqn{Z(x_i, l)}{Z(x_i,l)}
  is the observation at site \eqn{x_i}{x_i} and time
  \eqn{l}{l}. Taking this in mind and given the above definition of
  lorelogram, the spatio-temporal extention is straightforward.

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
  \item{distance}{The type of spatial distance;}
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
  Cooley, D., Naveau, P. and Poncet, P. (2006) \emph{Variograms for spatial
  max-stable random fields}. Dependence in Probability and
Statistics, p. 373--390.

Cressie, N. A. C. (1993) \emph{Statistics for Spatial Data}. New York: Wiley.

  Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modelling}.
  Spring Verlang, New York.

  Heagerty, P. J., and Zeger, S. L. (1998). Lorelogram: A Regression
  Approach to Exploring Dependence in Longitudinal Categorical
  Responses.
  \emph{Journal of the American Statistical Association},
  \bold{93}(441), 150--162

  Smith, R. L. (1990) Max-Stable Processes and Spatial Extremes.
  \emph{Unpublished manuscript}, University of North California.
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{http://staff.deuv.cl/bevilacqua/}.}


 \examples{
library(CompRandFld)
library(RandomFields)
set.seed(514)

# Set the coordinates of the sites:
x <- runif(200, 0, 10)
y <- runif(200, 0, 10)


################################################################
###
### Example 1. Empirical estimation of the variogram from a
### Gaussian random field with exponential correlation.
### One spatial replication is simulated.
###
###
###############################################################

# Set the model's parameters:
corrmodel <- "exponential"
mean <- 0
sill <- 1
nugget <- 0
scale <- 3

# Simulation of the spatial Gaussian random field:
data <- RFsim(x, y, corrmodel=corrmodel, param=list(mean=mean,
              sill=sill, nugget=nugget, scale=scale))$data

# Empirical spatial variogram estimation:
fit <- EVariogram(data, x, y)

# Results:
plot(fit$centers, fit$variograms, xlab='h', ylab=expression(gamma(h)),
     ylim=c(0, max(fit$variograms)), xlim=c(0, fit$srange[2]), pch=20,
     main="variogram")


################################################################
###
### Example 2. Empirical estimation of the variogram from a
### spatio-temporal Gaussian random fields with Gneiting
### correlation function.
### One spatio-temporal replication is simulated
###
###############################################################

set.seed(331)
# Define the temporal sequence:
times <- seq(1,10,1)

# Simulation of a spatio-temporal Gaussian random field:
data <- RFsim(x, y, times, corrmodel="gneiting",
              param=list(mean=0,scale_s=0.4,scale_t=1,sill=sill,
              nugget=0,power_s=1,power_t=1,sep=0.5))$data

# Empirical spatio-temporal variogram estimation:
fit <- EVariogram(data, x, y, times, maxtime=5,maxdist=4)

# Results: Marginal spatial empirical  variogram
par(mfrow=c(2,2), mai=c(.5,.5,.3,.3), mgp=c(1.4,.5, 0))
plot(fit$centers, fit$variograms, xlab='h', ylab=expression(gamma(h)),
     ylim=c(0, max(fit$variograms)), xlim=c(0, max(fit$centers)),
     pch=20,main="Marginal spatial Variogram",cex.axis=.8)

# Results: Marginal temporal empirical  variogram
plot(fit$bint, fit$variogramt, xlab='t', ylab=expression(gamma(t)),
     ylim=c(0, max(fit$variograms)),xlim=c(0,max(fit$bint)),
     pch=20,main="Marginal temporal Variogram",cex.axis=.8)

# Building space-time variogram
st.vario <- matrix(fit$variogramst,length(fit$centers),length(fit$bint))
st.vario <- cbind(c(0,fit$variograms), rbind(fit$variogramt,st.vario))

# Results: 3d Spatio-temporal variogram
require(scatterplot3d)
st.grid <- expand.grid(c(0,fit$centers),c(0,fit$bint))
scatterplot3d(st.grid[,1], st.grid[,2], c(st.vario),
              highlight.3d=TRUE, xlab="h",ylab="t",
              zlab=expression(gamma(h,t)), pch=20,
              main="Space-time variogram",cex.axis=.7,
              mar=c(2,2,2,2), mgp=c(0,0,0),
              cex.lab=.7)

# A smoothed version
par(mai=c(.2,.2,.2,.2),mgp=c(1,.3, 0))
persp(c(0,fit$centers), c(0,fit$bint), st.vario,
      xlab="h", ylab="u", zlab=expression(gamma(h,u)),
      ltheta=90, shade=0.75, ticktype="detailed", phi=30,
      theta=30,main="Space-time variogram",cex.axis=.8,
      cex.lab=.8)


################################################################
###
### Example 3. Empirical estimation of the madogram from a
### max-stable random field (Extremal Gaussian model) with
### exponential correlation.
### n iid spatial replications are simulated.
###
###############################################################

set.seed(7273)
# Simulation of the max-stable random field:
data <- RFsim(x, y, corrmodel=corrmodel, model="ExtGauss",
              param=list(mean=mean, sill=sill, nugget=nugget,
              scale=scale), replicates=40)$data
# Tranform data from from common unit Frechet to standard Gumbel margins:
data <- Dist2Dist(data, to='sGumbel')

# Empirical madogram estimation:
fit <- EVariogram(data, x, y, type='madogram', replicates=40, cloud=FALSE)

# Results:
par(mfrow=c(1,2), mai=c(.6,.6,.3,.3), mgp=c(1.6,.6, 0))
plot(fit$centers, fit$variograms, xlab='h', ylab=expression(nu(h)),
     ylim=c(0, max(fit$variograms)), xlim=c(0, fit$srange[2]), pch=20,
     main="madogram")
plot(fit$centers, fit$extcoeff, xlab='h', ylab=expression(theta(h)),
     ylim=c(1, 2), xlim=c(0, fit$srange[2]), pch=20,
     main="extremal coefficient")


################################################################
###
### Example 4. Empirical estimation of the lorelogram from a
### spatio-temporal Binary-Gaussian random field with double
### exponential correlation.
### One spatio-temporal replication is simulated.
###
###############################################################

set.seed(5351)
# Define the temporal sequence:
times <- seq(1,10,1)

# Simulation of othe space-time Binary Gaussian random field:

data <- RFsim(x, y, times, corrmodel="exp_exp", model="BinaryGauss",
        threshold=0, param=list(nugget=nugget,mean=mean,scale_t=1,
        scale_s=1,sill=1))$data

# Empirical lorelogram estimation:
fit <- EVariogram(data, x, y, times, type="lorelogram", maxtime=6,maxdist=4)

# Building space-time lorelogram
st.vario <- matrix(fit$variogramst,length(fit$centers),length(fit$bint))
tempvario <- rbind(fit$variogramt,st.vario)
st.vario <- cbind(c(NA,fit$variograms),tempvario)

par(mfrow=c(2,2))

# Results:
# Marginal spatial empirical lorelogram
plot(fit$centers, fit$variograms, xlab='h', ylab="lorelogram",
     xlim=c(0, max(fit$centers)), pch=20,main="Marginal spatial lorelgram")

# Marginal temporal empirical  lorelogram
plot(fit$bint, fit$variogramt, xlab='t', ylab="lorelogram",
     xlim=c(0,max(fit$bint)), pch=20,main="Marginal temporal lorelogram")

# 3d Spatio-temporal lorelogram:
require(scatterplot3d)
st.grid <- expand.grid(c(0,fit$centers), c(0,fit$bint))

scatterplot3d(st.grid[,1], st.grid[,2], c(st.vario),
              highlight.3d=TRUE, xlab = "h", ylab = "t",
              zlab="Lorelogram",pch=20, main="Space-time lorelogram",
              cex.axis=.7, mar=c(2,2,2,2), mgp=c(0,0,0),cex.lab=.7)

# smoothed version:
par(mai=c(.2,.2,.2,.2),mgp=c(1.2,.5, 0))
persp(c(0,fit$centers),c(0,fit$bint),st.vario,
      xlab = "h", ylab = "t", zlab ="Lorelogram",
      ltheta=90, shade=0.75, ticktype="detailed",
      phi=30,theta=30, main="Space-time lorelogram",
      cex.axis=.8, cex.lab=.8)
}

\keyword{Variogram}
