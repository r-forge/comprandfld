\name{FitComposite}
\alias{FitComposite}
\alias{print.FitComposite}

\title{Maximum Weighted Composite-likelihood Fitting of Gaussian Random Fields}
\description{
  Maximum weighted composite-likelihood fitting for Gaussian and max-stable random fields.
  The function returns the model parameters' estimates and the estimates'
  variances by weighted maximisation of the composite-likelihood and allows to
  fix any of the parameters.
}
\usage{
FitComposite(coordx, coordy=NULL, corrmodel, data, fixed=NULL, grid=FALSE,
             likelihood='Marginal', lonlat=FALSE, model='Gaussian',
             optimizer='Nelder-Mead', replicates=FALSE, start=NULL,
             type='Pairwise', varest=FALSE, vartype='SubSamp',
             weighted=FALSE, weights=NULL, winconst=NULL)
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
  \item{likelihood}{String; the configuration of the composite
    likelihood. \code{Marginal} is the default, see the Section
    \bold{Details}.}
  \item{lonlat}{Logical; if \code{FALSE} (the default), \code{coordx}
    and \code{coordy} are interpreted as Cartesian coordinates
    otherwise they are considered as longitude and latitude.}
  \item{model}{String; the type of random field and therefore the densities associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section
    \bold{Details}.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). 'Nelder-Mead' is the
    default.}
  \item{replicates}{Logical; if \code{FALSE} (the default) one spatial random field is
    considered, instead if \code{TRUE} the data are considered as iid replicates of a field.}
  \item{start}{A named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure. \code{NULL} is the default (see \bold{Details}).}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods (see \bold{Details}).}
  \item{varest}{Logical; if \code{TRUE} the estimate' variances and
    standard errors are returned.
    \code{FALSE} is the default.}
  \item{vartype}{String; (\code{Sub-Samp} the default) the type of
    method used for computing
  the estimates' variances, see the Section \bold{Details}.}
  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted, see the Section \bold{Details}. If \code{FALSE} (the
    default) the composite likelihood is not weighted.}
  \item{weights}{A numeric vector of weights (see \bold{Details}).}
  \item{winconst}{Numeric; a positive real value -- if
    \code{vartype=Sub-Samp} -- that determines the window size in the sub-sampling
    estimates of the variances (see \bold{Details}).}
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
    \item If it is a numeric vector, the data are interpreted as a spatial or temporal
    realisation; 
    \item If it is a numeric (\eqn{n \times d}{n x d})-matrix, the
  \code{d} columns represent the different points where the data
  observed and the \code{n} rows represent the available iid replications.
    \item If it is a numeric (\eqn{d \times d \times n}{d x d x
      n})-matrix then the data are considered as observations observed at
    (\eqn{d \times d}{d x d}) points with eventually \code{n} available replications.
  } 
  
  The \code{likelihood} parameter represents the composite-likelihood
  configurations. The settings alternatives are:
  \enumerate{
    \item \code{Conditional}, the composite-likelihood is formed by
    conditionals likelihoods;
    \item \code{Marginal}, the composite-likelihood is formed by
  marginals likelihoods; 
    \item \code{Full}, the composite-likelihood turns out to be the standard likelihood;
  }

  The \code{model} paramter represents the type of random field
  considered, an example is the Gaussian random field. Therefore, this also
  determines the associated density functions that will define the
  likelihood objects and then the composite-likelihood.
  The available options:
  \itemize{
    \item \code{Gaussian}, for a Gaussian random field.
  }

  The \code{replicates} parameter specifies if iid replications of the
  random field are available.
  If \code{FALSE}, then a realisation of spatial random field is considered. Instead if \code{TRUE}
  the \code{data} are considered as iid replications of a spatial random field. For the moment
  the option of iid replications is implemented only for a spatial
  random field. The alternative is to work with a realisation of a
  spatial-temporal random field.
 
   The \code{start} parameter allows to specify starting values.
   If \code{start} is omitted the routine is computing the
   starting values using the weighted moment estimator.

   The \code{type} parameter represents the type of likelihood used in the
  composite-likelihood definition. The possible alternatives are listed
  in the following scheme.
  \enumerate{
    \item If a Gaussian random field is considered
    (\code{model=Gaussian}):
    \itemize{
      \item If the composite is formed by marginal likelihoods (\code{likelihood=Marginal}):
      \itemize{
	\item \code{Pairwise}, the composite-likelihood is defined by
	the pairwise likelihoods;
	\item \code{Difference}, the composite-likelihood is defined by
	likelihoods 
	which are obtained as difference of the pairwise likelihoods.}
      \item If the composite is formed by conditional likelihoods
      (\code{likelihood=Conditional})
      \itemize{
	\item \code{Pairwise}, the composite-likelihood is defined by
	the pairwise conditional likelihoods.}
      \item If the composite is formed by a full likelihood
      (\code{likelihood=Full}):
      \itemize{
	\item \code{Standard}, the likelihood that will be considered is the standard
	multivariate version;
	\item \code{Restricted}, the likelihood that will be considered
        is the restricted version (see \bold{References}).}
      }
      }

  The \code{vartype} parameter -- if the \code{varest} is
  \code{TRUE} -- specifies the method used to compute the estimates' variances. The
  default \code{Sub-Samp} uses the Sub-Sampling method to estimate the
  variability matrix in the Godambe matrix. Other options are:
  \code{Theoretical} where for the variability matrix it is used the exact
  expression (for the moment it has been implemented only for the 
  \code{Difference} likelihood, see the \code{type} field) and \code{Sampling} where the variability matrix
  in this case is estimated by the sample contro-part (available only for iid.
  replicates of a random field, see the \code{replicates} parameter).

  The \code{weighted} parameter specifies if the likelihoods forming the
  composite-likelihood must be weighted. If \code{TRUE} the weights are
  selected by opportune procedures that improve the efficient of the
  maximum composite-likelihood estimator (not implemented yet). If
  \code{FALSE} the efficient improvement procedure is not used.

  The \code{weights} parameter allows to weight the composite-likelihood
  by weights insert by the users. These do not imply any gain in
  efficiency of the maximum composite-likelihood estimator but still be
  a reasonable setting (not implemented yet!).
}

\value{
  Returns an object of class \code{FitComposite}.
  An object of class \code{FitComposite} is a list containing
  at most the following components:
  
  \item{clic}{The composite information criterion, if the full
  likelihood is considered then it conicide with the Akaike information
  criterion;}
  \item{coord}{The vector of coordinates;}
  \item{convergence}{A string that denotes if convergence is reached;}
  \item{corrmodel}{The correlation model;}
  \item{data}{The vector or matrix of data;}
  \item{fixed}{The vector of fixed parameters;}
  \item{iterations}{The number of iteration used by the numerical routine;}
  \item{likelihood}{The configuration of the composite likelihood;}
  \item{logCompLik}{The value of the log composite-likelihood at the maximum;}
  \item{lonlat}{The type of coordinates;}
  \item{message}{Extra message passed from the numerical routines;}
  \item{model}{The density associated to the likelihood objects;}
  \item{param}{The vector of  parameters' estimates;}
  \item{stderr}{The vector of standard errors;}
  \item{sensmat}{The sensitivity matrix;}
  \item{varcov}{The matrix of the variance-covariance of the estimates;}
  \item{varimat}{The variability matrix;}
  \item{type}{The type of the likelihood objects.}
}


\references{
  
  Maximum Restricted Likelihood Estimator:
  
  Harville, D. A. (1977)
  Maximum Likelihood Approaches to Variance Component Estimation and to Related Problems.
  \emph{Journal of the American Statistical Association}, \bold{72},
  320--338.
   
  Composite-likelihood:
  
  Varin, C. and Vidoni, P. (2005)
  A Note on Composite Likelihood Inference and Model Selection.
  \emph{Biometrika}, \bold{92}, 519--528.

  Varin, C. (2008)
  On Composite Marginal Likelihoods.
  \emph{Advances in Statistical Analysis}, \bold{92}, 1--28.

  Weighted Composite-likelihood for max-stable random fields:
  
  Padoan, S. A. (2008). \emph{Computational
    Methods for Complex Problems in Extreme Value Theory}. PhD. Thesis,
  Department of Statistics, University of Padua.
  
  Padoan, S. A. Ribatet, M. and Sisson, S. A. (2010)
  Likelihood-Based Inference for Max-Stable Processes.
  \emph{Journal of the American Statistical Association, Theory & Methods}, \bold{105}, 263--277.

  Weighted Composite-likelihood for Gaussian random fields:

  Bevilacqua, M. Gaetan, C., Mateu, J. and Porcu, E. (2011)
  Estimating space and space-time covariance functions for large data
  sets: a weighted composite likelihood approach.
  \emph{Journal of the American Statistical Association}, to appear.

  Spatial Statistics:

  Schlather, M. (1999) \emph{An introduction to positive definite
    functions and to unconditional simulation of random fields}.
  Technical report ST 99--10, Dept. of Maths and Statistics, Lancaster University

  Wackernagel, H. (1998) \emph{Multivariate Geostatistics}. Berlin:
  Springer, 2nd edition.

  Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modelling}.
  Spring Verlang, New York.

  Spatial Extremes:

  Smith, R. L. (1990) Max-Stable Processes and Spatial Extremes.
  \emph{Unpublished manuscript}, University of North California.

  Schlather, M. (2002) Models for Stationary Max-Stable Random
  Fields. \emph{Extremes}, \bold{5}, 33--44.

  de Haan, L., and Pereira, T. T. (2006) Spatial Extremes: Models for
  the Stationary Case. \emph{The Annals of Statistics}, \bold{34},
  146--168.

  Kabluchko, Z., Schlather, M., and de Haan, L. (2009) Stationary max-stable fields
associated to negative definite functions. \emph{The Annals of
  Probability}, \bold{37}, 2042--2065.

  Kabluchko, Z. (2010) Extremes of Independent Gaussian
  Processes. \emph{Extremes}, to appear.

  Sub-sampling estimation:

  Heagerty, P. J. and Lumley T. (2000) Window Subsampling of Estimating
  Functions with Application to Regression Models.
  \emph{Journal of the American Statistical Association, Theory &
    Methods}, \bold{95}, 197--211.

  Lee, Y. D. and Lahiri S. N. (2002) Variogram Fitting by Spatial
  Subsampling.
  \emph{Journal of the Royal Statistical Society. Series B},
  \bold{64}, 837--854.
}

\seealso{\code{\link[RandomFields]{CovarianceFct}}, \code{\link{WLeastSquare}}, \code{\link{optim}}}


\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://people.epfl.ch/simone.padoan};
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

}

\keyword{Composite}
