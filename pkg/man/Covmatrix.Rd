\name{Covmatrix}
\alias{Covmatrix}

\title{Spatio-temporal (tapered) Covariance Matrix}
\description{
  The function computes the (tapered) covariance matrix for
  a spatial (temporal or spatio-temporal) covariance model
  and a set of spatial (temporal or spatio-temporal) points.
}
\usage{
Covmatrix(coordx, coordy=NULL, coordt=NULL, corrmodel,
          distance="Eucl", grid=FALSE, iskrig=FALSE,
          maxdist=NULL, maxtime=NULL, param ,
          taper=NULL, tapsep=NULL, type="Standard")
}

\arguments{
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of spatial sites) giving 2-dimensions of spatial coordinates or a numeric \eqn{d}{d}-dimensional vector giving
    1-dimension of spatial coordinates.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, the default is \code{NULL} then \code{coordx} is expected to
    be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
  \item{coordt}{A numeric vector giving 1-dimension of
    temporal coordinates. At the moment implemented only for the
    Gaussian case. Optional argument, the default is \code{NULL}
    then a spatial random field is expected.}
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
  \item{distance}{String; the name of the spatial distance. The default
  is \code{Eucl}, the euclidean distance. See \code{\link{FitComposite}}.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set
    of non-equispaced spatial sites (irregular grid). See
    \code{\link{FitComposite}}.}
  \item{iskrig}{Logical: the default value is \code{FALSE}. It is  \code{TRUE} if the function is called by  the function \code{\link{Kri}}.}
   \item{maxdist}{Numeric; an optional positive value indicating the
  marginal spatial compact support. See \code{\link{FitComposite}}.}
  \item{maxtime}{Numeric; an optional positive value indicating the
  marginal temporal compact support. See \code{\link{FitComposite}}.}
  \item{param}{A list of parameter values required  for the correlation
  model. See \code{\link{FitComposite}} and \code{\link{CorrelationParam}}.}
    \item{taper}{String; the name of the taper correlation
    function if type is  \code{Tapering}, see the Section \bold{Details}.}
   \item{tapsep}{Numeric; an optional  value indicating the separabe parameter in the space-time quasi taper (see \bold{Details}).}
   \item{type}{String; the type of covariance matrix
   \code{Standard} (the default) or \code{Tapering} for tapered covariance matrix}
   }


\value{
  Returns an object of class \code{CovMat}.
  An object of class \code{CovMat} is a list containing
  at most the following components:

  \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates;}
  \item{covmatrix}{The covariance matrix if \code{type} is\code{Standard}. An object of class spam if \code{type} is \code{Tapering}}
  \item{corrmodel}{String: the correlation model;}
  \item{distance}{String: the type of spatial distance;}
  \item{grid}{Logical:\code{TRUE} if the spatial data are in a regular grid, otherwise  \code{FALSE};}
  \item{nozero}{In the case of tapered matrix  the percentage of non zero values in the covariance matrix.
  Otherwise is NULL.}
  \item{maxdist}{Numeric: the marginal spatial compact support if \code{type} is \code{Tapering};}
  \item{maxtime}{Numeric: the marginal temporal compact support if \code{type} is \code{Tapering}; }
  \item{namescorr}{String: The names of the correlation parameters; }
  \item{numcoord}{Numeric: the number of spatial coordinates;}
  \item{numtime}{Numeric: the number the temporal  coordinates;}
  \item{param}{Numeric: The covariance parameters;}
  \item{tapmod}{String: the taper model if \code{type} is \code{Tapering}. Otherwise is NULL.}
  \item{spacetime}{\code{TRUE} if spatio-temporal and \code{FALSE} if spatial covariance model;}

  In the space-time case  \code{covmatrix} is the covariance matrix of the random vector \deqn{Z(s_1,t_1),Z(s_1,t_2),..Z(s_n,t_1),..,Z(s_n,t_m)}
  for \code{n} spatial locatione sites and \code{m} temporal instants.
}

\details{
The parameter \code{param} is a list including all the parameters of a
covariance function model.\cr
In particular, the covariance models share the following paramaters: the
\code{sill} that
represents the common variance
of the random field,
the \code{nugget} that represents the local variation (white
noise) at the origin.
For each correlation model you can check the list of the specific parameters using \code{\link{CorrelationParam}}.\cr
Here there is the list of all the implemented space and space-time
correlation models. The list of space-time correlation functions
includes separable
and non-separable models.

\itemize{
  \item Purerly spatial correlation models:
  \enumerate{
    \item \code{cauchy }
    \deqn{R(h) = \left(1+h^2\right)^{-\beta}}{R(h) = (1+h^2)^(-\beta)}
    The parameter \eqn{\beta} is positive. It is a special case of the  \code{gencauchy} model.

    \item \code{exponential}
    \deqn{R(h)  =e^{-h},  \quad h\ge0}{R(h) = exp(-h)}
    This model is a special case of the \code{whittle} and the \code{stable} model.

    \item \code{gauss}
    \deqn{R(h) = e^{-h^2}}{R(h)=exp(-h^2)}
    This model is a special case of the \code{stable} model.

    \item \code{gencauchy} (generalised \code{cauchy})\cr

    \deqn{R(h) = ( 1+h^\alpha )^{-\frac{\beta}{\alpha}}}
    The parameter \eqn{\alpha} is in (0,2], and \eqn{\beta}
    is positive. \cr


  \item \code{spherical}
    \deqn{R(h) = (1- 1.5 h+0.5 h^3)
      1_{[0,1]}(h)}
    This isotropic covariance function is valid only for dimensions
    less than or equal to 3.


    \item \code{stable}
    \deqn{R(h) = e^{-h^\alpha},  \quad h\ge0}{R(h)=exp(-h ^\alpha)}
    The parameter \eqn{\alpha} is in \eqn{(0,2]}.


    \item \code{wave}
    \deqn{R(h)=\frac{\sin h}h, \quad h>0 \qquad \hbox{and } R(0)=1}{
      R(h)=sin(h)/h if h>0 and C(0)=1}
    This isotropic covariance function is valid only for dimensions less
    than or equal to 3.

    \item \code{matern}
    \deqn{R(h) = 2^{1-\nu} \Gamma(\nu)^{-1} x^\nu
      K_\nu(h)}
    The parameter \eqn{\nu} is positive.
    \cr
    This is the model of choice if the smoothness of a random field is to
    be parametrised: if \eqn{\nu > m} then the
    graph is \eqn{m} times differentiable.
  }

  \item Spatio-temporal correlation models:

  \itemize{
    \item Non-separable models:
    \enumerate{
      \item \code{gneiting} (non-separabel space time model)
    \deqn{R(h, u) = \frac{e^{ \frac{-h^\nu} { (1+u^\lambda)^{0.5 \gamma \nu }}}}  {1+u^\lambda}}
    The parameters \eqn{\nu} and \eqn{\lambda} take values
    in \eqn{[0,2]}; the parameter \eqn{\gamma} take values in \eqn{[0,1]}.
    For \eqn{\gamma=0} it is a separable model.

	\item \code{gneiting_GC} (non-separabel space time model with great circle distances)
    \deqn{R(h, u) = \frac{ e^{  \frac{-u^\lambda }{1+h^\nu)^{0.5 \gamma \lambda}} }} { 1+h^\nu} }


    \item \code{iacocesare} (non-separabel space time model)
    \deqn{R(h, u) = (1+h^\nu+u^\lambda)^{-\delta}}
    The parameters \eqn{\nu} and \eqn{\lambda} take values
    in \eqn{[1,2]}; the parameters \eqn{\delta} must be greater
    than or equal to half the space-time dimension.


       \item \code{porcu} (non-separabel space time model)
    \deqn{R(h, u) = (0.5 (1+h^\nu)^\gamma +0.5 (1+u^\lambda)^\gamma)^{-\gamma^{-1}}}
    The parameters \eqn{\nu} and \eqn{\lambda} take values
    in \eqn{[0,2]}; the paramete \eqn{\gamma} take values in \eqn{[0,1]}.
    The limit of the correlation model  as \eqn{\gamma} tends to zero  leads to a  separable model.

       \item \code{porcu2} (non-separabel space time model)
    \deqn{R(h, u) =\frac{ e^{ -h^\nu ( 1+u^\lambda)^{0.5 \gamma \nu}}} { (1+u^\lambda)^{1.5}}}
    The parameters \eqn{\nu} and \eqn{\lambda} take values
    in \eqn{[0,2]}; the parameter \eqn{\gamma } take values in \eqn{[0,1]}.
    For \eqn{\gamma=0} it is a separable model.
      }
      \item Separable models.\cr

      Space-time separable correlation models are easly obtained as the product of a spatial and a temporal correlation model, that is
 \deqn{R(h,u)=R(h) R(u)}
 Several combinations are possible:
    \enumerate{
    \item  \code{exp_exp}: spatial exponential model and  temporal  exponential model
     \item  \code{exp_cauchy}: spatial exponential model and  temporal cauchy  model
      \item  \code{matern_cauchy}: spatial matern model and  temporal cauchy  model
        \item  \code{stable_stable}: spatial stabel model and temporal stable  model
  }

  Note that some models are nested. (The \code{exp_exp} with the  \code{stable_stable} for instance.)
  }

  \item Spatial taper function models.\cr

  For spatial covariance tapering the tapered correlation functions are:
    \enumerate{
        \item \code{Wendland1}
    \deqn{R(h) = (1-h)^2 (1+0.5 h) 1_{[0,1]}(h)}

       \item \code{Wendland2}
    \deqn{R(h) = (1-h)^4 (1+4 h) 1_{[0,1]}(h)}

       \item \code{Wendland3}
    \deqn{R(h) = (1-h)^6 (1+6 h + 35 h^2 /3) 1_{[0,1]}(h)}
      }
      \item Spatio-tempora tapered correlation models.\cr

      For space-time covariance tapering likelihood the taper functions
      are obtained
      as the product of a spatial and a temporal taper (Separable taper).
  Several combinations are possible:
    \itemize{\item  \code{Wendlandi_Wendlandj}: spatial \code{Wendlandi}
      taper and temporal \code{Wendlandj} taper with i,j=1,2,3.}

    \item Space-time non separable adaptive-taper with dynamically space-time compact support is:
      \itemize{\item \code{qt_time} and \code{qt_space}. In The case of \code{qt_time} the space-time quasi taper is:
      \deqn{T(h,u) = (arg)^{-6} (1+7 x) (1-x)^7 1_{[0,\frac{maxtime}{arg}]}(u)}
      \deqn{arg=(1+\frac{h}{maxdist} )^\beta, x=u \frac{ arg}{maxtime} }
      }
      where \eqn{0<=\beta<=1} is a fixed parameter of separability (\code{tapsep}), \code{maxtime } the fixed temporal compact support and \code{maxdist} the
      fixed spatial scale parameter.
      The adaptive-taper \code{qt_space} is the same taper but changing the time with the space.

}

 \bold{Remarks}:\cr

   Let \code{R(h)} be a  spatial correlation model given in standard notation.
   Then the covariance model
   applied with arbitrary variance and scale equals to:
  \deqn{C(h)=sill +nugget ,   \quad if \quad h=0 }    \cr
  \deqn{C(h)=sill * R( \frac{h}{scale},...)  ,  \quad if \quad h>0 }
  Similarly if \code{R(h,u)} is a spatio-temporal correlation model given in standard notation,
   then the covariance model is:
  \deqn{C(h,u)=sill +nugget ,  \quad if \quad  h=0,u=0 }    \cr
  \deqn{C(h,u)=sill * R( \frac{h}{scale_s},\frac{u}{scale_t},...)  ,  \quad if \quad h>0 \quad or \quad u>0 }


  Here \sQuote{...} stands for additional parameters.






    Let \code{R(h)} be a spatial taper  given in standard notation.
   Then the taper function applied with an arbitrary compact support (maxdist) equals to:
  \deqn{T(h)= R( \frac{h}{maxdist})}
  Similarly if \code{R(h,u)} is a spatio-temporal taper given in standard notation, then the taper
  function applied with  arbitrary compact supports (maxdist, maxtime)  equals to:

  \deqn{T(h,u)= R( \frac{h}{maxdist},\frac{u}{maxtime})}


  Then the tapered covariance matrix is obtained as:
       \deqn{C_{tap}(h,u)= T(h,u)C(h,u) }


}


\references{

Bevilacqua, M., Mateu, J., Porcu, E., Zhang, H. and Zini, A. (2010). Weighted composite
 likelihood-based tests for space-time separability of covariance functions. \emph{Statistics and Computing}, \bold{20(3)},
283-293.

Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modelling}.
  Spring Verlang, New York.

  Gneiting, T. (2002). Nonseparable, stationary covariance functions for
  space-time data.
  \emph{Journal of the American Statistical Association}, \bold{97},
  590--600.

  Gneiting, T., Genton, M. G. and Guttorp, P. (2007).
  \emph{Geostatistical space-time models, stationarity, separability and full
  symmetry}.
  In Finkenstadt, B., Held, L. and Isham, V. (eds.),
  Statistical Methods for Spatio-Temporal Systems, Chapman & Hall/CRC,
  Boca Raton, pp. 151-175


  Schlather, M. (1999) \emph{An introduction to positive definite
    functions and to unconditional simulation of random fields}.
  Technical report ST 99--10, Dept. of Maths and Statistics, Lancaster University

}

\seealso{
  \code{\link{Kri}},
  \code{\link{RFsim}},
   \code{\link{FitComposite}}}



\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
   \url{http://staff.deuv.cl/bevilacqua/}.}

\examples{
library(CompRandFld)
library(spam)

################################################################
###
### Example 1. Covariance matrix associated to
### a Matern correlation model
###
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(500, 0, 2)
y <- runif(500, 0, 2)

matrix1 <- Covmatrix(x, y, corrmodel="matern", param=list(smooth=0.5,
                    sill=1,scale=0.2))
dim(matrix1$covmatrix)

################################################################
###
### Example 2. Tapered Covariance matrix associated to
### a Matern correlation model
###
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(500, 0, 2)
y <- runif(500, 0, 2)

matrix2 <- Covmatrix(x, y, corrmodel="matern", param=list(smooth=0.5,
                    sill=1,scale=0.2),maxdist=0.3,taper="Wendland1",
                    type="Tapering")
# Tapered covariance matrix
as.matrix(matrix2$covmatrix)[1:15,1:15]

# Percentage of no zero values in the tapered matrix
matrix2$nozero

################################################################
###
### Example 3. Covariance matrix associated to
### a space-time double exponential correlation model
###
###############################################################

# Define the temporal-coordinates:
times <- seq(1, 5, 1)

# Define correlation model
corrmodel="exp_exp"

# Define covariance parameters
param=list(scale_s=0.3,scale_t=0.5,sill=1)

# Simulation of a spatial Gaussian random field:
matrix3 <- Covmatrix(x, y, times, corrmodel=corrmodel,
                     param=param)

dim(matrix3$covmatrix)

################################################################
###
### Example 4. Tapered Covariance matrix associated to
### a space-time double exponential correlation model
###
###############################################################

# Simulation of a spatial Gaussian random field:
matrix4 <- Covmatrix(x, y, times, corrmodel=corrmodel, param=param, maxdist=0.3,
                     maxtime=2,taper="Wendland2_Wendland2",type="Tapering")

# Tapered space time covariance matrix
as.matrix(matrix4$covmatrix)[1:10,1:10]


}
\keyword{Simulation}
