\name{FitComposite}
\alias{FitComposite}
\alias{print.FitComposite}

\title{Max-Likelihood-Based Fitting of Gaussian, Binary and Max-Stable Fields}
\description{
  Maximum weighted composite-likelihood fitting of spatio-temporal
  Gaussian, binary and spatial max-stable random fields.
  For the spatio-temporal Gaussian random field, (restricted)
  maximum likelihood and tapered likelihood fitting is also avalable.
  The function returns the model parameters' estimates and the estimates'
  variances and allows to fix any of the parameters.
}
\usage{
FitComposite(data, coordx, coordy=NULL, coordt=NULL, corrmodel,
             distance='Eucl', fixed=NULL, grid=FALSE,
             likelihood='Marginal', margins='Gev', maxdist=NULL,
             maxtime=NULL, model='Gaussian', optimizer='Nelder-Mead',
             replicates=1, start=NULL, taper=NULL, tapsep=NULL,
             threshold=NULL,type='Pairwise', varest=FALSE,
             vartype='SubSamp', weighted=FALSE, winconst, winstp)
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
  \item{distance}{String; the name of the spatial distance. The default is \code{Eucl}, the euclidean distance.
    See the Section \bold{Details}.}
  \item{fixed}{An optional named list giving the values of the parameters that
    will be considered as known values. The listed parameters for a
    given correlation function will be not estimated, i.e. if
    \code{list(nugget=0)} the nugget effect is ignored.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set of non-equispaced spatial sites (irregular grid).}
  \item{likelihood}{String; the configuration of the composite
    likelihood. \code{Marginal} is the default, see the Section
    \bold{Details}.}
  \item{margins}{String; the type of the marginal distribution of the
    max-stable field. \code{Gev} is the default, see the Section
    \bold{Details}.}
  \item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial distance considered in the composite or tapered likelihood
    computation. See the Section \bold{Details} for more information.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal separation considered in the composite or tapered likelihood
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
  \item{taper}{String; the name of the type of covariance matrix.
     It can be \code{Standard} (the default value) or
    \code{Tapering} for taperd covariance matrix.}
   \item{tapsep}{Numeric; an optional  value indicating the separabe parameter in the space time quasi taper (see \bold{Details}).}
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
  \item{winconst}{Numeric; a positive value for computing the sub-window
    size where observations are sampled in the sub-sampling procedure (if \code{vartype=SubSamp}).
    For increasing \code{winconst}
    increasing sub-window sizes are obtained. Optional argument, the
   default is \eqn{1}. See \bold{Details} for more information.}
 \item{winstp}{Numeric; a value in \eqn{(0,1]} for computing the
   sub-window step (in the sub-sampling procedure). This value denote the proportion
   of the sub-window size. Optional argument, the default is
   \eqn{0.5}. See \bold{Details} for more information.}
}

\details{
  Note, that the standard likelihood may be seen as particular case of the
  composite likelihood. In this respect \code{FitComposite} provides maximum (restricted) likelihood
  fitting. Only composite likelihood estimation based on pairs are considered. Specifically marginal pairwise,
  conditional pairwise and difference pairwise. Covariance tapering is considered only for Gaussian random fields.

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
    \item If \code{data} is a numeric (\eqn{d \times d \times n}{d x d x n})-array, \code{coordx} and \code{coordy} are two
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
  function for the random field. (See   \code{\link{Covmatrix}} ).

  The \code{distance} parameter allows to consider differents kinds of spatial distances.
  The settings alternatives are:
  \enumerate{
    \item \code{Eucl}, the euclidean distance (default value);
    \item \code{Chor}, the chordal distance;
    \item \code{Geod}, the geodesic distance;
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

    In the case of tapering likelihood \code{maxdist} and \code{maxtime} describes the spatial and temporal compact support
    of the taper model (see \code{\link{Covmatrix}}). If they are not specified then the maximum spatial and temporal distances are considered.
    In the case of space time quasi taper the \code{tapsep} parameter allows to specify the spatio temporal compact support (see \code{\link{Covmatrix}}).


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
    \item \code{ExtGauss}, for an Extremal Gaussian max-stable random
  field (known also as Schlather model) (see Schlather, M. 2002);
    \item \code{ExtT}, for an Extremal \code{t} max-stable random field (see Davison, A. C. et al. 2012);
  }
  Note, that only for the \code{Gaussian} case the estimation procedure is implemented for spatial and spatial-temporal
  random fields.

   The \code{start} parameter allows to specify starting values.
   If \code{start} is omitted the routine is computing the
   starting values using the weighted moment estimator.

   The \code{taper} parameter, optional in case that
   \code{type=Tapering}, indicates the type of taper
   correlation model. (See \code{\link{Covmatrix}})

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
        version of the full likelihood (e.g. Kaufman et al. 2008, see \bold{References}).
      }
      }
      }


  The \code{varest} parameter specifies if the standard error estimation of the estimated parameters must be computed.
  For Gaussian random field and standard (restricted) likelihood estimation, standard  errrors  are computed as square root of the diagonal elements of the Fisher
  Information matrix (asymptotic covariance matrix of the estimates under increasing domain).
  For Gaussian random field and tapered and composite  likelihood estimation, standard errors estimate are computed as square root of the diagonal elements of the
  Godambe Information matrix. (asymptotic covariance matrix of the estimates under increasing domain  (see Shaby, B. and D. Ruppert (2012) for tapering and
  Bevilacqua et. al. (2012) , Bevilacqua and Gaetan (2013) for weighted composite likelihood)).
  The \code{vartype} parameter specifies the method used to compute the estimates' variances in the composite likelihood case. In particular
  for estimating the variability matrix \code{J} in the Godambe expression matrix.
  This parameter is considered if \code{varest=TRUE}. The options are:
  \itemize{
    \item \code{SubSamp} (the default), indicates the Sub-Sampling method;
    \item \code{Sampling}, indicates that the variability matrix
    is estimated by the sample contro-part (available only for \eqn{n}{n} iid
    replications of the random field, i.e. \code{replicates=n});
  }

  The \code{weighted} parameter specifies if the likelihoods forming the
  composite-likelihood must be weighted. If \code{TRUE} the weights are
  selected by opportune procedures that improve the efficient of the
  maximum composite-likelihood estimator (not implemented yet). If
  \code{FALSE} the efficient improvement procedure is not used.

  For computing the standard errors by the sub-sampling procedure,
  \code{winconst} and \code{winstp} parameters represent respectively a positive constant used to
  determine the sub-window size and the the step with which the
  sub-window moves.

  In the spatial case (subset of \eqn{R^2}{R^2}), the domain is seen as
  a rectangle \eqn{B \times H}{BxH}, therefore the size of the
  sub-window side \eqn{b}{b} is given by \eqn{b=winconst \times
  \sqrt(B)}{b=winconst x sqrt(B)} (similar is of \eqn{h}{h}).
  For a complete description see Lee and Lahiri (2002).
  By default \code{winconst} is set \eqn{B / (2 \times \sqrt(B))}{B / (2
  x sqrt(B))}.
  The \code{winstp} parameter is used to determine the sub-window step. The latter is given by the
  proportion of the sub-window size, so that when \code{winstp=1} there
  is not overlapping between contiguous sub-windows.
  In the spatial case by default \code{winstp=0.5}.
  The sub-window is moved
  by successive steps in order to cover the entire spatial
  domain.
  Observations, that fall in disjoint or overlapping windows
  are considered indipendent samples.

  In the spatio-temporal case the subsampling is meant only in time as
  described by Li et al. (2007). Thus, \code{winconst} represents
  the lenght of the temporal sub-window. By default the size of the
  sub-window is computed following the rule established in Li et al. (2007).
  By default \code{winstp} is the time step.

  Observe that in the spatio-temporal case, the returned values by
  \code{srange}
  and \code{trange}, represent respectively the minimum and maximum
  of the marginal spatial distances and those of the temporal
  separations. Thus, the minimum being not the overall (i.e. considering
  the spatio-temporal coordinates) is not zero, as one could be expect and
  the latter can be easily added by the user.
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
  \item{distance}{The type of spatial distance;}
  \item{fixed}{The vector of fixed parameters;}
  \item{iterations}{The number of iteration used by the numerical routine;}
  \item{likelihood}{The configuration of the composite likelihood;}
  \item{logCompLik}{The value of the log composite-likelihood at the maximum;}
  \item{message}{Extra message passed from the numerical routines;}
  \item{model}{The density associated to the likelihood objects;}
   \item{nozero}{In the case of tapered likelihood  the percentage of non zero values in the covariance matrix.
  Otherwise is NULL.}
  \item{numcoord}{The number of spatial coordinates;}
  \item{numrep}{The number of the iid replicatations of the random field;}
  \item{numtime}{The number the temporal realisations of the random field;}
  \item{param}{The vector of  parameters' estimates;}
  \item{srange}{The minimum and maximum spatial distance (see \bold{Details}). The maximum
  is \code{maxdist}, if inserted, rather the effective maximum distance;}
  \item{stderr}{The vector of standard errors;}
  \item{sensmat}{The sensitivity matrix;}
  \item{varcov}{The matrix of the variance-covariance of the estimates;}
  \item{varimat}{The variability matrix;}
  \item{vartype}{The method used to compute the variance of the estimates;}
  \item{trange}{The minimum and maximum temporal separation (see \bold{Details}). The maximum
    is \code{maxtime}, if inserted, rather then the effective maximum
  separation;}
  \item{threshold}{The threshold used in the binary random field.}
  \item{type}{The type of the likelihood objects.}
  \item{winconst}{The constant use to compute the window size in the
  sub-sampling procedure;}
  \item{winstp}{The step used for moving the window in the sub-sampling procedure}
}


\references{

  Maximum Restricted Likelihood Estimator:

  Harville, D. A. (1977)
  Maximum Likelihood Approaches to Variance Component Estimation and to Related Problems.
  \emph{Journal of the American Statistical Association}, \bold{72},
  320--338.

  Tapered likelihood:

    Kaufman, C. G., Schervish, M. J. and Nychka, D. W. (2008)
  Covariance Tapering for Likelihood-Based Estimation in Large Spatial Dataset.
  \emph{Journal of the American Statistical Association}, \bold{103},
  1545--1555.

    Shaby, B. and D. Ruppert (2012).
  Tapered covariance: Bayesian estimation and asymptotics. \emph{J. Comp. Graph.
Stat.}, \bold{21-2}, 433--452.


  Composite-likelihood:

  Varin, C., Reid, N. and Firth, D. (2011). An Overview of Composite
  Likelihood Methods. \emph{Statistica Sinica}, \bold{21}, 5--42.

  Varin, C. and Vidoni, P. (2005)
  A Note on Composite Likelihood Inference and Model Selection.
  \emph{Biometrika}, \bold{92}, 519--528.


  Weighted Composite-likelihood for binary random fields:

   Patrick, J. H. and Subhash, R. L. (1998)
  A Composite Likelihood Approach to Binary Spatial Data.
  \emph{Journal of the American Statistical Association, Theory & Methods}, \bold{93}, 1099--1111.


  Weighted Composite-likelihood for max-stable random fields:

  Davison, A. C. and Gholamrezaee, M. M. (2012)
  Geostatistics of extremes.
  \emph{Proceedings of the Royal Society of London, series A}, \bold{468}, 581--608.

  Padoan, S. A. (2008). \emph{Computational
  Methods for Complex Problems in Extreme Value Theory}. PhD. Thesis,
  Department of Statistics, University of Padua.

  Padoan, S. A. Ribatet, M. and Sisson, S. A. (2010)
  Likelihood-Based Inference for Max-Stable Processes.
  \emph{Journal of the American Statistical Association, Theory & Methods}, \bold{105}, 263--277.

  Weighted Composite-likelihood for Gaussian random fields:

    Bevilacqua, M. Gaetan, C., Mateu, J. and Porcu, E. (2012)
  Estimating space and space-time covariance functions for large data
  sets: a weighted composite likelihood approach.
  \emph{Journal of the American Statistical Association, Theory &
    Methods}, \bold{107}, 268--280.

  Bevilacqua, M. Gaetan, C. (2013)
  On composite likelihood inference based on pairs for
spatial Gaussian random fields
  \emph{Techical Report, Department of Statistics, de Valparaiso University}.

  Spatial Extremes:

  Davison, A. C., Padoan, S. A., and Ribatet, M. (2012) Statistical
  Modelling of Spatial Extremes, with discussion.
  \emph{Statistical Science}, \bold{27}, 161--186.

  de Haan, L., and Pereira, T. T. (2006) Spatial Extremes: Models for
  the Stationary Case. \emph{The Annals of Statistics}, \bold{34},
  146--168.

  Kabluchko, Z. (2010) Extremes of Independent Gaussian
  Processes. \emph{Extremes}, \bold{14}, 285--310.

  Kabluchko, Z., Schlather, M., and de Haan, L. (2009) Stationary max-stable fields
  associated to negative definite functions. \emph{The Annals of
    Probability}, \bold{37}, 2042--2065.

  Schlather, M. (2002) Models for Stationary Max-Stable Random
  Fields. \emph{Extremes}, \bold{5}, 33--44.

  Smith, R. L. (1990) Max-Stable Processes and Spatial Extremes.
  \emph{Unpublished manuscript}, University of North California.

  Sub-sampling estimation:

  Carlstein, E. (1986) The Use of Subseries Values for Estimating the
  Variance.
  \emph{The Annals of Statistics},
  \bold{14}, 1171--1179.

  Heagerty, P. J. and Lumley T. (2000) Window Subsampling of Estimating
  Functions with Application to Regression Models.
  \emph{Journal of the American Statistical Association, Theory &
    Methods}, \bold{95}, 197--211.

  Lee, Y. D. and Lahiri S. N. (2002) Variogram Fitting by Spatial
  Subsampling.
  \emph{Journal of the Royal Statistical Society. Series B},
  \bold{64}, 837--854.

  Li, B., Genton, M. G. and Sherman, M. (2007). A nonparametric assessment of properties of space-time
  covariance functions. \emph{Journal of the American Statistical
    Association}, \bold{102}, 736--744
}

\seealso{
  \code{\link{Covmatrix}},
  \code{\link[RandomFields]{GaussRF}},
  \code{\link[RandomFields]{MaxStableRF}}, \code{\link{WLeastSquare}}, \code{\link{optim}}}


\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{http://staff.deuv.cl/bevilacqua/}.}

\examples{
library(CompRandFld)
library(RandomFields)
library(spam)
set.seed(3132)

###############################################################
############ Examples of spatial random fields ################
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(100, 0, 10)
y <- runif(100, 0, 10)

# Set the covariance model's parameters:
corrmodel <- "exponential"
mean <- 0
sill <- 1
nugget <- 0
scale <- 1.5
param<-list(mean=mean,sill=sill,nugget=nugget,scale=scale)
coords=cbind(x,y)
# Simulation of the spatial Gaussian random field:
data <- RFsim(coordx=coords, corrmodel=corrmodel, param=param)$data

# Fixed parameters
fixed<-list(mean=mean,nugget=nugget)

# Starting value for the estimated parameters
start<-list(scale=scale,sill=sill)


################################################################
###
### Example 1. Maximum likelihood fitting of
### Gaussian random fields with exponential correlation.
### One spatial replication.
### Likelihood setting: composite with
### marginal pairwise likelihood objects.
###
###############################################################


# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(data, coordx=coords, corrmodel=corrmodel, maxdist=5,
                    likelihood="Marginal",type="Pairwise",start=start,fixed=fixed)

# Results:
print(fit)

################################################################
###
### Example 2. Maximum likelihood fitting of
### Gaussian random fields with exponential correlation.
### One spatial replication.
### Likelihood setting: standard full likelihood.
###
###############################################################

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(data, coordx=coords, corrmodel=corrmodel,likelihood="Full",
                    type="Standard",varest=TRUE,start=start,fixed=fixed)
# Results:
print(fit)

################################################################
###
### Example 3. Maximum likelihood fitting of
### Gaussian random fields with exponetial correlation.
### One spatial replication.
### Likelihood setting: tapered full likelihood.
###
###############################################################

# Maximum tapered likelihood fitting of the random field:
fit <- FitComposite(data, coordx=coords, corrmodel=corrmodel,likelihood="Full",
                    type="Tapering",taper="Wendland1",maxdist=4,
                    varest=TRUE,start=start,fixed=fixed)

# Results:
print(fit)


################################################################
###
### Example 4. Maximum composite-likelihood fitting of
### max-stable random fields. Extremal Gaussian model with
### exponential correlation. n iid spatial replications.
### Likelihood setting: composite with marginal pairwise
### likelihood objects.
###
###############################################################

# Simulation of a max-stable random field in the specified points:
data <- RFsim(x, y, corrmodel=corrmodel, model="ExtGauss", replicates=30,
              param=list(mean=mean,sill=sill,nugget=nugget,scale=scale))$data

# Maximum composite-likelihood fitting of the random field:
fit <- FitComposite(data, x, y, corrmodel=corrmodel, model="ExtGauss",
                    replicates=30, varest=TRUE, vartype="Sampling",
                    margins="Frechet",start=list(sill=sill,scale=scale))

# Results:
print(fit)

################################################################
###
### Example 5. Maximum likelihood fitting of
### Binary-Gaussian random fields with exponential correlation.
### One spatial replication.
### Likelihood setting: composite with marginal pairwise
### likelihood objects.
###
###############################################################

#set.seed(3128)

#x <- runif(200, 0, 10)
#y <- runif(200, 0, 10)

# Simulation of the spatial Binary-Gaussian random field:
#data <- RFsim(coordx=coords, corrmodel=corrmodel, model="BinaryGauss",
#              threshold=0, param=list(mean=mean,sill=.8,
#              nugget=nugget,scale=scale))$data

# Maximum composite-likelihood fitting of the random field:
#fit <- FitComposite(data, coordx=coords, corrmodel=corrmodel, threshold=0,
#                    model="BinaryGauss", fixed=list(nugget=nugget,
#                    mean=0),start=list(scale=.1,sill=.1))

# Results:
#print(fit)


###############################################################
######### Examples of spatio-temporal random fields ###########
###############################################################

# Define the temporal sequence:
#time <- seq(1, 80, 1)

# Define the spatial-coordinates of the points:
#x <- runif(10, 0, 10)
#y <- runif(10, 0, 10)
#coords=cbind(x,y)

# Set the covariance model's parameters:
#corrmodel="exp_exp"
#scale_s=0.5
#scale_t=1
#sill=1
#nugget=0
#mean=0

#param<-list(mean=0,scale_s=1,scale_t=1,sill=sill,nugget=nugget)

# Simulation of the spatial-temporal Gaussian random field:
#data <- RFsim(coordx=coords,coordt=time,corrmodel=corrmodel,
#              param=param)$data

# Fixed parameters
#fixed<-list(mean=mean,nugget=nugget)

# Starting value for the estimated parameters
#start<-list(scale_s=scale_s,scale_t=scale_t,sill=sill)

################################################################
###
### Example 6. Maximum likelihood fitting of
### Gaussian random field with double-exponential correlation.
### One spatio-temporal replication.
### Likelihood setting: composite with conditional pairwise
### likelihood objects.
###
###############################################################

# Maximum composite-likelihood fitting of the random field:
#fit <- FitComposite(data=data,coordx=coords,coordt=time,corrmodel="exp_exp",
#                    maxtime=2,maxdist=1,likelihood="Marginal",type="Pairwise",
#                    start=start,fixed=fixed)

# Results:
#print(fit)

}
\keyword{Composite}
