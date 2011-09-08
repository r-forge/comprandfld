\name{FitComposite}
\alias{FitComposite}
\alias{print.FitComposite}

\title{Max-Likelihood Fitting of Gaussian, Binary and Max-Stable Random Fields}
\description{
  Maximum weighted composite-likelihood fitting of Gaussian, binary and max-stable random fields.
  The function returns the model parameters' estimates and the estimates'
  variances by weighted maximisation of the composite-likelihood and allows to
  fix any of the parameters.
}
\usage{
FitComposite(data, coordx, coordy=NULL, coordt=NULL, corrmodel,
             fixed=NULL, grid=FALSE, likelihood='Marginal',
             lonlat=FALSE, margins='Gev', maxdist=NULL, maxtime=NULL,
             model='Gaussian', optimizer='Nelder-Mead', replicates=1,
             start=NULL, taper=NULL, threshold=NULL, type='Pairwise',
             varest=FALSE, vartype='SubSamp', weighted=FALSE,
             weights=NULL, winconst=NULL)
}
\arguments{
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{n \times d}{n x d})-matrix
   (\eqn{n} iid spatial realisations) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid)
   or an (\eqn{d \times d \times n}{d x d x n})-array (\eqn{n} iid spatial realisations on regular grid) or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation) or an (\eqn{t \times d \times n }{t x d x n})-array
   (\eqn{n} iid spatial-temporal realisations) or or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid) or an (\eqn{d \times d \times t \times n }{d x d x t x n})-array
   (\eqn{n} iid spatial-temporal realisations on regular grid).
   For the description see the Section \bold{Details}.}
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of spatial sites) assigning 2-dimensions of spatial coordinates or a numeric \eqn{d}{d}-dimensional vector assigning
    1-dimension of spatial coordinates.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    spatial coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, the default is \code{NULL} then \code{coordx} is expected to
    be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
  \item{coordt}{A numeric vector assigning 1-dimension of
    temporal coordinates. At the moment implemented only for the
    Gaussian case. Optional argument, the default is \code{NULL}
    then a spatial random field is expected.}
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
  \item{fixed}{An optional named list giving the values of the parameters that
    will be considered as known values. The listed parameters for a
    given correlation function will be not estimated, i.e. if
    \code{list(nugget=0)} the nugget effect is ignored.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set of non-equispaced spatial sites (irregular grid).}
  \item{likelihood}{String; the configuration of the composite
    likelihood. \code{Marginal} is the default, see the Section
    \bold{Details}.}
  \item{lonlat}{Logical; if \code{FALSE} (the default), \code{coordx}
    and \code{coordy} are interpreted as Cartesian coordinates
    otherwise they are considered as longitude and latitude.}
  \item{margins}{String; the type of the marginal distribution of the
    max-stable field. \code{Gev} is the default, see the Section
    \bold{Details}.}
  \item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial distance considered in the composite-likelihood
    computation. See the Section \bold{Details} for more information.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal separation considered in the composite-likelihood
    computation (see \bold{Details}).}
  \item{model}{String; the type of random field and therefore the densities associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section
    \bold{Details}.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). 'Nelder-Mead' is the
    default.}
  \item{replicates}{Numeric; a positive integer denoting the number of independent and identically distributed (iid)
    replications of a spatial or spatial-temporal random field. Optional argument, the default value is \eqn{1} then
    a single realisation is considered.}
  \item{start}{An optional named list with the initial values of the
    parameters that are used by the numerical routines in maximization
    procedure. \code{NULL} is the default (see \bold{Details}).}
  \item{taper}{String; the name of the taper correlation
    function, see the Section \bold{Details}. Optional
    parameter considered only if the \code{type} parameter is equal to
    \code{Tapering}.}
  \item{threshold}{Numeric; a value indicating a threshold for the
    binary random field. Optional in the case that \code{model}
    is \code{BinaryGauss}, see the Section \bold{Details}.}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods (see \bold{Details}).}
  \item{varest}{Logical; if \code{TRUE} the estimates' variances and
    standard errors are returned.
    \code{FALSE} is the default.}
  \item{vartype}{String; (\code{SubSamp} the default) the type of
    method used for computing
  the estimates' variances, see the Section \bold{Details}.}
  \item{weighted}{Logical; if \code{TRUE} the likelihood objects are
    weighted, see the Section \bold{Details}. If \code{FALSE} (the
    default) the composite likelihood is not weighted.}
  \item{weights}{A numeric vector of weights (see \bold{Details}). Optional argument, the default value is \code{NULL}
   then no weights are considered}
  \item{winconst}{Numeric; a positive real value indicating the window size used from the
   sub-sampling method for the estimation of the parameters' variances. This window size
   is considered only if \code{vartype=SubSamp}. Optional argument, the default value is \eqn{1}
   (see \bold{Details}).}
}

\details{
  With \code{data}, \code{coordx}, \code{coordy}, \code{coordt}, \code{grid} and \code{replicates} parameters:
  \itemize{
    \item If \code{data} is a numeric \eqn{d}{d}-dimensional vector, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors (or \code{coordx} is (\eqn{d \times 2}{d x 2})-matrix and \code{coordy=NULL}),
    \code{coordt=NULL}, \code{grid=FALSE} and \code{replicates=1}, then the data are interpreted as a single spatial
    realisation observed on \eqn{d}{d} spatial sites;
    \item If \code{data} is a numeric (\eqn{n \times d}{n x d})-matrix, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors (or \code{coordx} is (\eqn{d \times 2}{d x 2})-matrix and \code{coordy=NULL}),
    \code{coordt=NULL}, \code{grid=FALSE} and \code{replicates=n}, then the data are interpreted as \code{n} iid replications of a spatial
    random field observed on \code{d} spatial sites.
    \item If \code{data} is a numeric (\eqn{d \times d}{d x d})-matrix, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors, \code{coordt=NULL}, \code{grid=TRUE} and \code{replicates=1}, then the data are interpreted as
    a single spatial random field realisation observed on \code{d} equispaced spatial sites (named regular grid).
    \item If \code{data} is a numeric (\eqn{d \times d \time n}{d x d x n})-array, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors, \code{coordt=NULL}, \code{grid=TRUE} and \code{replicates=n}, then the data are interpreted as
    \code{n} iid realisations of a spatial random field observed on \code{d} equispaced spatial sites.
    \item If \code{data} is a numeric (\eqn{t \times d}{t x d})-matrix, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors (or \code{coordx} is (\eqn{d \times 2}{d x 2})-matrix and \code{coordy=NULL}),
    \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, \code{grid=FALSE} and \code{replicates=1}, then the data
    are interpreted as a single spatial-temporal realisation of a random field observed on \code{d} spatial sites and for \code{t} times.
    \item If \code{data} is a numeric (\eqn{t \times d \times n}{t x d x n})-array, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors (or \code{coordx} is (\eqn{d \times 2}{d x 2})-matrix and \code{coordy=NULL}),
    \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, \code{grid=FALSE} and \code{replicates=n}, then the data
    are interpreted as \code{n} iid realisations of a spatial-temporal random field observed on \code{d} spatial sites and for \code{t} times.
    \item If \code{data} is a numeric (\eqn{d \times d \times t}{d x d x t})-array, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors, \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, \code{grid=TRUE} and
    \code{replicates=1}, then the data are interpreted as a single spatial-temporal realisation of a random field observed on
    \code{d} equispaced spatial sites and for \code{t} times.
    \item If \code{data} is a numeric (\eqn{d \times d \times t \times n}{d x d x t x n})-array, \code{coordx} and \code{coordy} are two
    numeric \eqn{d}{d}-dimensional vectors, \code{coordt} is a numeric \eqn{t}{t}-dimensional vector, \code{grid=TRUE} and
    \code{replicates=n}, then the data are interpreted as \code{n} iid realisation of a spatial-temporal random field observed on
    \code{d} equispaced spatial sites and for \code{t} times.
  }

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

  The \code{likelihood} parameter represents the composite-likelihood
  configurations. The settings alternatives are:
  \enumerate{
    \item \code{Conditional}, the composite-likelihood is formed by
    conditionals likelihoods;
    \item \code{Marginal}, the composite-likelihood is formed by
  marginals likelihoods;
    \item \code{Full}, the composite-likelihood turns out to be the standard likelihood;
  }

  The \code{margins} parameter concerns only max-stable fields and indicates how the margins are
  considered. The options are \code{Gev} or \code{Frechet}, where in the
  former case the marginals are supposed generalized
  extreme value distributed and in the latter case unit Frechet
  distributed.

  The \code{maxdist} parameter set the maximum
  spatial distance below which pairs of sites with inferior distances
  are considered in the composite-likelihood. This can be
  inferior of the effective maximum spatial distance. \bold{Note} that
  this corresponds to use a weighted composite-likelihood with binary
  weights. Pairs with distance less than \code{maxdist} have weight 1
  and are included in the likelihood computation, instead those with
  greater distance have weight 0 and then excluded.
  The default
  is \code{NULL}, in this case the effective maximum spatial distance
  between sites is considered.

  The same arguments of \code{maxdist} are valid for \code{maxtime} but
  here the weigthed composite-likelihood regards the case of
  spatial-temporal field. At the moment is
  implemented only for Gaussian random fields. The default
  is \code{NULL}, in this case the effective maximum temporal lag
  between pairs of observations is considered.

  The \code{model} paramter indicates the type of random field
  considered, for instance \code{model=Gaussian} denotes a Gaussian random field.
  Accordingly, this also determines the analytical expression of the finite dimensional distribution associated with the random field.
  The available options are:
  \itemize{
    \item \code{Gaussian}, for a Gaussian random field (see
    i.e. Wackernagel, H. 1998);
    \item \code{BinaryGauss}, for a Binary random field (see Heagerty
    and Lele 1998)
    \item \code{BrowResn}, for a Brown-Resnick max-stable random field (see Kabluchko, Z. et al. 2009);
    \item \code{ExtGauss}, for a Extremal Gaussian max-stable random field (see Schlather, M. 2002);
    \item \code{ExtT}, for a Extremal \code{t} max-stable random field (see Davison, A. C. et al. 2011);
  }
  Note, that only for the \code{Gaussian} case the estimation procedure is implemented for spatial and spatial-temporal
  random fields.

   The \code{start} parameter allows to specify starting values.
   If \code{start} is omitted the routine is computing the
   starting values using the weighted moment estimator.

   The \code{taper} parameter, optional in case that
   \code{type=Tapering}, indicates the type of taper
   correlation model. At the moment the implemented models are:
   \enumerate{
     \item Wendland type 1;
     \item Wendland type 2;
     \item Wendland type 3.
   }
   For more details e.g. Kaufman et al. (2008) (see
   \bold{References}).
   \bold{ATTENTION}! Insert also a value for the
   \code{maxdist} parameter in order to set the upper bound of the compact support; if
   no \code{maxdist} is inserted then the effective maximum spatial
   distance is considered and therefore the tapering approach does not
   have any effect.

   The \code{threshold} parameter indicates the value (common for all
   the spatial sites) above which the values of the underlying Gaussian latent process
   are considered sucesses events (values below are instead
   failures). See e.g. Heagerty and Lele (1998) for more details.

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
	\item \code{Standard}, the objective function is the classical
        multivariate likelihood;
	\item \code{Restricted}, the objective function is the
        restricted version of the full likelihood (e.g. Harville 1977, see \bold{References});
	\item \code{Tapering}, the objective function is the tapered
        version of the full likelihood (e.g. Kaufman et al. 2008, see \bold{References}). At
        the moment is implemented only for spatial Gaussian random fields.
      }
      }
      }

  The \code{vartype} parameter specifies the method used to compute the estimates' variances. In particular
  for estimating the variability matrix \code{J} in the Godambe expression matrix.
  This parameter is considered if \code{varest=TRUE}. The options are:
  \itemize{
    \item \code{SubSamp} (the default), indicates the Sub-Sampling method;
    \item \code{Theoretical}, denotes that the variability matrix is computed using the exact
    analytical expression. At the moment is implemented for spatial
    Gaussian random fields when using difference composite likelihood,
    namely \code{Model='Gaussian'} and \code{likelihood='Difference'};
    \item \code{Sampling}, indicates that the variability matrix
    is estimated by the sample contro-part (available only for \eqn{n}{n} iid
    replications of the random field, i.e. \code{replicates=n});
  }

  The \code{weighted} parameter specifies if the likelihoods forming the
  composite-likelihood must be weighted. If \code{TRUE} the weights are
  selected by opportune procedures that improve the efficient of the
  maximum composite-likelihood estimator (not implemented yet). If
  \code{FALSE} the efficient improvement procedure is not used.

  The \code{weights} parameter allows to weight the composite-likelihood
  by weights insert by the users. These do not imply any gain in
  efficiency of the maximum composite-likelihood estimator but still be
  a reasonable setting (not implemented yet).

  The \code{winconst} parameter set the size of the sub spatial window into which
  observations are sampled. The window is moved covering the entire spatial
  domain. Observations belonging to disjoint or overlapping windows are
  considered indipendent. If you do not insert a window size the rule
  described by Lee and Lahiri (2002) is followed to determine
  the dimension.
  For the case of spatial-temporal Gaussian random fields,
  you can insert the size of a temporal window, if it is not inserted the rule of
  Carlstein (1986) is followed.

}

\value{
  Returns an object of class \code{FitComposite}.
  An object of class \code{FitComposite} is a list containing
  at most the following components:

  \item{clic}{The composite information criterion, if the full
  likelihood is considered then it coincides with the Akaike information
  criterion;}
  \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates;}
  \item{convergence}{A string that denotes if convergence is reached;}
  \item{corrmodel}{The correlation model;}
  \item{data}{The vector or matrix or array of data;}
  \item{fixed}{The vector of fixed parameters;}
  \item{iterations}{The number of iteration used by the numerical routine;}
  \item{likelihood}{The configuration of the composite likelihood;}
  \item{logCompLik}{The value of the log composite-likelihood at the maximum;}
  \item{lonlat}{The type of coordinates;}
  \item{message}{Extra message passed from the numerical routines;}
  \item{model}{The density associated to the likelihood objects;}
  \item{numcoord}{The number of spatial coordinates;}
  \item{numrep}{The number of the iid replicatations of the random field;}
  \item{numtime}{The number the temporal realisations of the random field;}
  \item{param}{The vector of  parameters' estimates;}
  \item{srange}{The minimum and maximum spatial distance;}
  \item{stderr}{The vector of standard errors;}
  \item{sensmat}{The sensitivity matrix;}
  \item{varcov}{The matrix of the variance-covariance of the estimates;}
  \item{varimat}{The variability matrix;}
  \item{trange}{The minimum and maximum temporal separation;}
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

  Patrick, J. H. and Subhash, R. L. (1998)
  A Composite Likelihood Approach to Binary Spatial Data.
  \emph{Journal of the American Statistical Association, Theory & Methods}, \bold{93}, 1099--1111.

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

  Kaufman, C. G., Schervish, M. J. and Nychka, D. W. (2008)
  Covariance Tapering for Likelihood-Based Estimation in Large Spatial Dataset.
  \emph{Journal of the American Statistical Association}, \bold{103}, 1545--1555.

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
  Processes. \emph{Extremes}, \bold{14}, 285--310.

  Davison, A. C., Padoan, S. A., and Ribatet, M. (2011) Statistical Modelling of Spatial Extremes.
  \emph{Statistical Science}, to appear.

  Sub-sampling estimation:

  Heagerty, P. J. and Lumley T. (2000) Window Subsampling of Estimating
  Functions with Application to Regression Models.
  \emph{Journal of the American Statistical Association, Theory &
    Methods}, \bold{95}, 197--211.

  Lee, Y. D. and Lahiri S. N. (2002) Variogram Fitting by Spatial
  Subsampling.
  \emph{Journal of the Royal Statistical Society. Series B},
  \bold{64}, 837--854.

  Carlstein, E. (1986) The Use of Subseries Values for Estimating the
  Variance.
  \emph{The Annals of Statistics},
  \bold{14}, 1171--1179.
}

\seealso{\code{\link[RandomFields]{CovarianceFct}},
  \code{\link[RandomFields]{GaussRF}},
  \code{\link[RandomFields]{MaxStableRF}}, \code{\link{WLeastSquare}}, \code{\link{optim}}}


\author{Simone Padoan, \email{simone.padoan@unibg.it},
  \url{http://www.unibg.it/pers/?simone.padoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@unibg.it}.}

\examples{
library(RandomFields)
library(spam)
set.seed(2111)

# Set the coordinates of the points:
x <- runif(100, 0, 20)
y <- runif(100, 0, 20)


################################################################
###
### Example 1. Maximum composite-likelihood fitting of
### Gaussian random fields with Stable correlation model
### One spatial replication.
### Composite setting: difference marginal likelihoods.
###
###############################################################

# Set the model's parameters:
corrmodel <- "stable"
mean <- 0
variance <- 1
nugget <- 0
scale <- 10
power <- 1.2

# Simulation of the Gaussian random field in the specified points:
data <- GaussRF(x=x, y=y, model=corrmodel, grid=FALSE,
                param=c(mean, variance, nugget, scale, power))

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(data, x, y, corrmodel=corrmodel,
                    type='Difference', varest=TRUE)

# Results:
print(fit)

################################################################
###
### Example 2. Maximum likelihood fitting of
### Gaussian random fields with Stable correlation model.
### One spatial replication.
### Likelihood setting: standard full likelihood.
###
###############################################################


# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(data, x, y, corrmodel=corrmodel,
                    likelihood='Full', type='Standard', varest=TRUE)

# Results:
print(fit)

################################################################
###
### Example 3. Maximum likelihood fitting of
### Gaussian random fields with Stable correlation model.
### One spatial replication.
### Likelihood setting: tapered full likelihood.
###
###############################################################

# Note that this example is really time consuming!
# However faster than the standard full likelihood!

# Set the coordinates of the points:
#x <- runif(1000, 0, 40)
#y <- runif(1000, 0, 40)

# Simulation of the Gaussian random field in the specified points:
#data <- GaussRF(x=x, y=y, model=corrmodel, grid=FALSE,
#                param=c(mean, variance, nugget, scale, power))

# Maximum composite-likelihood fitting of the random field:
#fit <- FitComposite(data, x, y, corrmodel=corrmodel,
#                    likelihood='Full', type='Tapering',
#                    taper="Wendland1", maxdist=15)

# Results:
#print(fit)

################################################################
###
### Example 4. Maximum composite-likelihood fitting of
### max-stable random fields. Extremal Gaussian model with
### Stable correlation function. n iid spatial replications.
### Composite-likelihood setting: pairwise marginal likelihoods.
###
###############################################################

# Simulation of a max-stable random field in the specified points:
data <- MaxStableRF(x=x, y=y, model=corrmodel, grid=FALSE, maxstable="extr",
                   param=c(mean, variance, nugget, scale, power),
                   n=20,pch='')
data <- t(data)

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(data, x, y, corrmodel=corrmodel, model='ExtGauss',
                    replicates=20, varest=TRUE, vartype='Sampling',
                    margins="Frechet")

# Results:
print(fit)

}

\keyword{Composite}
