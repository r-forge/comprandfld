\name{FitComposite}
\alias{FitComposite}
\alias{print.FitComposite}

\title{Maximum Composite-likelihood Fitting of Random Fields}
\description{
  Maximum composite-likelihood fitting for random fields.
  The function returns the parameters' estimates and the estimates'
  variances of random fields obtained by maximisation of the
  composite-likelihood and allows to fix any of the parameters.
}
\usage{
FitComposite(coordx, coordy=NULL, corrmodel, data, fixed=NULL,
             grid=FALSE, hessian=FALSE, likelihood='Marginal',
             lonlat=FALSE, model='Gaussian', optimizer='Nelder-Mead',
             start=NULL, varest=FALSE, time=FALSE, type='Pairwise',
             weighted=FALSE, weights=NULL)
} 
\arguments{
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of points) assigning 2-dimensions of coordinates or a numeric vector assigning
    1-dimension of coordinates.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector otherwise it will be ignored.}
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations (see
    \bold{Details}).}
  \item{fixed}{A named list giving the values of the parameters that
    will be considered as known values. The listed parameters for a
    given correlation function will be not estimated, i.e. if
    \code{list(nugget=0)} the nugget effect is ignored.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{hessian}{Logical; if \code{FALSE} (the default) the hessian
  matrix is not computed from the optimizer routine.}
  \item{likelihood}{String; the configuration of the composite
    likelihood. \code{Marginal} is the default, see the Section
    \bold{Details}.}
  \item{lonlat}{Logical; if \code{FALSE} (the default), \code{coordx}
    and \code{coordy} are interpreted as Cartesian coordinates
    otherwise they are considered as longitude and latitude.}
  \item{model}{String; the density associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section
    \bold{Details}.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). 'Nelder-Mead' is the default.}
  \item{start}{A named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure. \code{NULL} is the default (see \bold{Details}).}
  \item{varest}{Logical; if \code{TRUE} the estimate' variances and
    standard errors are returned. 
    \code{FALSE} is the default.}
  \item{time}{Logical; if \code{FALSE} (the default) a spatial random field is
    considered (one temporal realisation), if \code{TRUE} a
    spatial-temporal random field is considered, see the Section \bold{Details}.}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods (see \bold{Details}).}
  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted, see the Section \bold{Details}. If \code{FALSE} (the
    default) the composite likelihood is not weighted.}
  \item{weights}{A numeric vector of weights (see \bold{Details}).}
}

\details{
  The \code{corrmodel} parameter allows to select a specific correlation
  function for the random field. The implemented correlation models
  are:
  \enumerate{
    \item \code{cauchy};
    \item \code{exponential};
    \item \code{gauss} (Gaussian);
    \item \code{gencauchy} (generalised Cauchy);
    \item \code{stable} (or powered exponential);
    \item \code{whittlematern} (Whittle-Matern).
  }
  See for more details \code{\link[RandomFields]{CovarianceFct}}.

  With the \code{data} parameter:
  \itemize{
    \item If a numeric vector, the data are interpreted as one spatial
    realisation; 
    \item If a numeric (\eqn{n \times d}{n x d})-matrix, the columns
    represent the data observed at different points and
    the rows represent the data for different time steps.
    \item If a numeric (\eqn{d \times d \times n}{d x d x n})-matrix the data are
    observed at (\eqn{d \times d}{d x d}) points for
    \code{n} time steps.
  } 
  
  The \code{likelihood} parameter represents the composite-likelihood
  configurations. The settings alternatives are:
  \enumerate{
    \item \code{Conditional}, the composite-likelihood is formed by
    conditionals likelihoods (not implemented yet);
    \item \code{Marginal}, the composite-likelihood is formed by
  marginals likelihoods; 
    \item \code{Full}, the composite-likelihood turns out to be the standard likelihood;
  }

  The \code{model} paramter represents the density function underlying
  the definition of the likelihoods which form the composite-likelihood.
  The settings alternatives are:
  \itemize{
    \item \code{Gaussian}, the Gaussian density.
   }   
 
   The \code{start} parameter allows to specify starting values.
   If \code{start} is omitted the routine is computing the
  starting values using the weighted moment estimator.

  The \code{time} parameter allows to specify the type of random field.
  If \code{FALSE} a spatial random field is considered, if \code{TRUE}
  a spatial-temporal random field is used. For the moment
  the case of i.i.d. time replications is implemented. Soon will be possible to
  specify also dependence structure for the temporal component.
  

  The \code{type} parameter represents the type of likelihood used in the
  composite-likelihood definition. The settings alternatives are:
  \enumerate{
    \item If the composite is formed by \code{Marginal} likelihoods:
    \itemize{
      \item If each likelihood is obtained by the \code{Gaussian}
      density then with:
      \itemize{
      \item \code{Pairwise}, the composite-likelihood is defined by
      the pairwise likelihoods;
      \item \code{Difference}, the composite-likelihood is defined by
      likelihoods which are obtained as difference of the pairwise likelihoods.}
    }
    
    \item If the \code{Full} likelihood is considered:
    \itemize{
      \item If the likelihood is obtained by the \code{Gaussian}
      density then with:
      \itemize{
      \item \code{Standard}, the likelihood used is the standard version;
      \item \code{Restricted}, the likelihood used is the restricted version.}
    }
  }

  The \code{weighted} parameter specifies if the likelihoods forming the
  composite-likelihood must be weighted. If \code{TRUE} the weights are
  selected by opportune procedures that improve the efficient of the
  maximum composite-likelihood estimator (not implemented yet). If
  \code{FALSE} the efficient improvement procedure is not used.

  The \code{weights} parameter allows to weight the composite-likelihood
  by weights insert by the users. These do not imply any gain in
  efficiency of the maximum composite-likelihood estimator but still be
  a reasonable setting.
}

\value{
  The returned object is a list with:
  \itemize{
    \item \code{clic}, the composite information criterion, if the full likelihood is considered then it conicide with the Akaike information criterion;
    \item \code{coord}, the vector of coordinates;
    \item \code{convergence}, a string that denotes if convergence is reached;
    \item \code{corrmodel}, the correlation model;
    \item \code{data}, the vector or matrix of data;
    \item \code{fixed}, the vector of fixed parameters;
    \item \code{iterations}, the number of iteration used by the numerical routine;
    \item \code{likelihood}, the configuration of the composite likelihood;
    \item \code{logCompLik}, the value of the log composite-likelihood
    at the maximum;
    \item \code{lonlat}, the type of coordinates;
    \item \code{message}, extra message passed from the numerical routines;
    \item \code{model}, the density associated to the likelihood objects;
    \item \code{param}, the vector of  parameters' estimates;
    \item \code{stderr}, the vector of standard errors;
    \item \code{sensmat}, the sensitivity matrix;
    \item \code{varcov}, the matrix of the variance-covariance of the estimates;
    \item \code{varimat}, the variability matrix;
    \item \code{type}, the type of the likelihood objects.
    }
}

\references{
  Harville, D. A. (1977)
  Maximum Likelihood Approaches to Variance Component Estimation and to Related Problems.
  \emph{Journal of the American Statistical Association}, \bold{72}, 320--338.
  
  Varin, C. and Vidoni, P. (2005)
  A Note on Composite Likelihood Inference and Model Selection.
  \emph{Biometrika}, \bold{92}, 519--528.

  Varin, C. (2008)
  On Composite Marginal Likelihoods.
  \emph{Advances in Statistical Analysis}, \bold{92}, 1--28.

  Padoan, S. A. Ribatet, M and Sisson, S. A. (2010)
  Likelihood-Based Inference for Max-Stable Processes.
  \emph{Journal of the American Statistical Association, Theory & Methods}, \bold{105}, 263--277.
}

\seealso{\code{\link[RandomFields]{CovarianceFct}}, \code{\link{WLeastSquare}}, \code{\link{optim}}}


\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php};
  Moreno Bevilacqua, \email{moreno.bevilacqua@unive.it}.}

\examples{
library(RandomFields)
set.seed(2111)

# Set the coordinates of the points:
x <- runif(100, 0, 20) 
y <- runif(100, 0, 20)


################################################################
###
### Example 1. Maximum composite-likelihood fitting of one
### spatial realisation of a Gaussian random field.
### Composite-likelihood setting: pairwise marginal likelihoods.
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

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(x, y, corrmodel, sim)

# Results:
print(fit)


################################################################
###
### Example 2. Maximum composite-likelihood fitting of one
### spatial realisation of a Gaussian random field.
### Composite-likelihood setting: difference likelihoods.
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

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(x, y, corrmodel, sim, type='Difference')

# Results:
print(fit)


################################################################
###
### Example 3. Maximum likelihood fitting of one
### spatial realisation of a Gaussian random field.
### Likelihood setting: restricted likelihood.
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

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(x, y, corrmodel, sim, likelihood='Full',
                    type='Restricted')

# Results:
print(fit)



################################################################
###
### Example 4. Maximum composite-likelihood fitting of n i.i.d.
### spatial realisations of a Gaussian random field.
### Composite-likelihood setting: difference likelihoods.
###
###############################################################

numdata <- 30

# Simulation of the Gaussian random field in the specified points:
sim <- GaussRF(x=x, y=y, model=corrmodel, grid=FALSE, n = numdata,
               param=c(mean, variance, nugget, scale, power), pch='')

sim <- t(sim)

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(x, y, corrmodel, sim, time=TRUE, type='Difference')


# Results:
print(fit)


################################################################
###
### Example 5. Maximum composite-likelihood fitting of one
### spatial realisations of a Gaussian random field on a
### regular grid.
### Composite-likelihood setting: difference likelihoods.
###
###############################################################

step <- 1 
x <- seq(0, 20, step) 
y <- seq(0, 20, step)

# Simulation of the Gaussian random field in the specified points:
sim <- GaussRF(x=x, y=y, model=corrmodel, grid=TRUE,
               param=c(mean, variance, nugget, scale, power))

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(x, y, corrmodel, sim, grid=TRUE, type='Difference')

# Results:
print(fit)


################################################################
###
### Example 6. Maximum composite-likelihood fitting of n i.i.d.
### spatial realisations of a Gaussian random field on a
### regular grid.
### Composite-likelihood setting: difference likelihoods.
###
###############################################################

numdata <- 5

# Simulation of the Gaussian random field in the specified points:
sim <- GaussRF(x=x, y=y, model=corrmodel, grid=TRUE, n = numdata,
               param=c(mean, variance, nugget, scale, power), pch='')

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(x, y, corrmodel, sim, grid=TRUE, time=TRUE,
                    type='Difference')

# Results:
print(fit)


}

\keyword{Composite}
