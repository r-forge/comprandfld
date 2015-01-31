\name{WLeastSquare}
\alias{WLeastSquare}
\alias{print.WLS}

\title{WLS of Gaussian and Max-Stable Random Fields}
\description{
  the function returns the parameters' estimates and the estimates' variances
  of a random field obtained by the weigthed least squares estimator.
}
\usage{
WLeastSquare(data, coordx, coordy=NULL, coordt=NULL, corrmodel,
             distance="Eucl", fixed=NULL, grid=FALSE, maxdist=NULL,
             maxtime=NULL,  model='Gaussian', optimizer='Nelder-Mead',
             numbins=NULL, replicates=1, start=NULL, weighted=FALSE)
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
    \code{d} is the number of spatial sites) giving 2-dimensions of spatial coordinates or a numeric \eqn{d}{d}-dimensional vector giving
    1-dimension of spatial coordinates.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, the default is \code{NULL} then \code{coordx} is expected to
    be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
  \item{coordt}{A numeric vector giving 1-dimension of
    temporal coordinates. Optional argument, the default is \code{NULL} then a spatial random field is expected.}
  \item{corrmodel}{String; the name of a correlation model, for the
    description (see \code{\link{FitComposite}}).}
  \item{distance}{String; the name of the spatial distance. The default is \code{Eucl},
                  the euclidean distance. See the Section \bold{Details}
		  of \code{\link{FitComposite}}.}
  \item{fixed}{A named list giving the values of the parameters that
    will be considered as known values. The listed parameters for a
    given correlation function will be not estimated, i.e. if
    \code{list(nugget=0)} the nugget effect is ignored.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{maxdist}{A numeric value denoting the maximum distance, see
    \bold{Details} and \code{\link{FitComposite}}.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal lag considered in the composite-likelihood
    computation (see \code{\link{FitComposite}}.}
  \item{model}{String; the type of random field. \code{Gaussian} is the default, see
    \code{\link{FitComposite}} for the different types.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). 'Nelder-Mead' is the
    default.}
  \item{numbins}{A numeric value denoting the numbers of bins, see the
    Section \bold{Details}}
   \item{replicates}{Numeric; a positive integer denoting the number of independent and identically distributed (iid)
    replications of a spatial or spatial-temporal random field. Optional argument, the default value is \eqn{1} then
    a single realisation is considered.}
  \item{start}{A named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure. \code{NULL} is the default (see \code{\link{FitComposite}}).}
  \item{weighted}{Logical; if \code{TRUE} then the weighted least square
    estimator is considered. If \code{FALSE} (the
    default) then the classic least square is used.}
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
   Returns an object of class \code{WLS}.
  An object of class \code{WLS} is a list containing
  at most the following components:

  \item{bins}{Adjacent intervals of grouped distances;}
  \item{bint}{Adjacent intervals of grouped temporal separations}
  \item{centers}{The centers of the bins;}
  \item{coordx}{The vector or matrix of spatial coordinates;}
  \item{coordy}{The vector of spatial coordinates;}
  \item{coordt}{The vector of temporal coordinates;}
  \item{convergence}{A string that denotes if convergence is reached;}
  \item{corrmodel}{The correlation model;}
  \item{data}{The vector or matrix of data;}
  \item{distance}{The type of spatial distance;}
  \item{fixed}{The vector of fixed parameters;}
  \item{iterations}{The number of iteration used by the numerical routine;}
  \item{message}{Extra message passed from the numerical routines;}
  \item{model}{The type of random fields;}
  \item{numcoord}{The number of spatial coordinates;}
  \item{numrep}{The number of the iid replicatations of the random field;}
  \item{numtime}{The number the temporal realisations of the random field;}
  \item{param}{The vector of  parameters' estimates;}
  \item{srange}{The minimum and maximum spatial distance;}
  \item{trange}{The minimum and maximum temporal separations;}
  \item{variograms}{The empirical spatial variogram;}
  \item{variogramt}{The empirical temporal variogram;}
  \item{variogramst}{The empirical spatial-temporal variogram;}
  \item{weighted}{A logical value indicating if its the weighted
    method;}
  \item{wls}{The value of the least squares at the minimum.}

}

\references{

Padoan, S. A. and Bevilacqua, M. (2015). Analysis of Random Fields Using CompRandFld.
  \emph{Journal of Statistical Software}, \bold{63}(9), 1--27.

  Barry, J. T., Crowder, M. J. and Diggle, P. J. (1997)
  \emph{Parametric estimation of the variogram}. Tech. Report, Dept Maths &
  Stats, Lancaster University.

  Cressie, N. A. C. (1993) \emph{Statistics for Spatial Data}. New York:
  Wiley.

  Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modelling}.
  Spring Verlang, New York.

  Smith, R. L. (1990) Max-Stable Processes and Spatial Extremes.
  \emph{Unpublished manuscript}, University of North California.
}

\seealso{\code{\link{FitComposite}}, \code{\link{optim}}}

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home}.}


\examples{
library(CompRandFld)
library(RandomFields)
set.seed(2111)

# Set the coordinates of the sites:
x <- runif(100, 0, 10)
y <- runif(100, 0, 10)


################################################################
###
### Example 1. Least square fitting of a Gaussian random field
### with exponential correlation.
### One spatial replication is simulated.
### Unweighted version (all weights equals to 1).
###
###############################################################

# Set the model's parameters:
corrmodel <- "exponential"
mean <- 0
sill <- 1
nugget <- 0
scale <- 2

# Simulation of the Gaussian random field:
data <- RFsim(x, y, corrmodel=corrmodel, param=list(mean=mean,
              sill=sill, nugget=nugget, scale=scale))$data

fix<-list(nugget=0)
ini<-list(scale=scale,sill=sill)
# Least square fitting of the random field:
fit <- WLeastSquare(data, x, y, corrmodel=corrmodel,fixed=fix,start=ini)

# Results:
print(fit)


################################################################
###
### Example 2. Least square fitting of a max-stable random field
### (Extremal Gaussian model) with exponential correlation
### n iid spatial replications.
### Unweighted version (all weights equals to 1).
###
###############################################################


# Simulation of the max-stable random field:
data <- RFsim(x, y, corrmodel=corrmodel, model="ExtGauss",
              param=list(mean=mean, sill=sill, nugget=nugget,
              scale=scale), replicates=40)$data

# Least square fitting of the random field:
fit <- WLeastSquare(data, x, y, corrmodel=corrmodel, model="ExtGauss",
                    replicates=40)

# Results:
print(fit)

################################################################
###
### Example 3. Least square fitting of a spatio-temporal
### Gaussian random field with double exponential correlation.
### One replication is simulated.
### Weighted version (all weights equals to 1).
###
###############################################################

# Define the temporal sequence:
#time <- seq(1, 25, 1)

# Simulation of the Gaussian random field:
#data <- RFsim(x, y, time, corrmodel="exp_exp", param=list(mean=mean,
#             scale_s=scale,scale_t=1,sill=sill,nugget=nugget))$data

#fix<-list(nugget=nugget)
#ini<-list(scale_s=scale,scale_t=1,sill=1)
# Weighted least square estimation:
#fit <- WLeastSquare(data, x, y, time, corrmodel="exp_exp", maxdist=5,
#                    maxtime=5,fixed=fix,start=ini)

# Results
#print(fit)
}

\keyword{LeastSquare}
