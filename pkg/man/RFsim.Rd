\name{RFsim}
\alias{RFsim}
\alias{print.RFsim}

\title{Simulation of Gaussian, Binary and Max-stable Random Fields}
\description{
  Simulation of spatial and spatio-temporal Gaussian, binary and max-stable random fields.
  The function returns one or more replications of a random field for a given covariance model and covariance parameters.
}
\usage{
RFsim(coordx, coordy=NULL, coordt=NULL, corrmodel, distance="Eucl",
      grid=FALSE, model='Gaussian', numblock=NULL, param,
      replicates=1, threshold=NULL)
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
    is \code{Eucl}, the euclidean distance. See the Section
    \bold{Details}  of \code{\link{FitComposite}}.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as spatial or spatial-temporal realisations on a set
    of non-equispaced spatial sites (irregular grid).}
  \item{model}{String; the type of random field and therefore the densities associated to the likelihood
    objects. \code{Gaussian} is the default, see the Section
    \bold{Details}.}
  \item{numblock}{Numeric; the observation size of the underlying random
    field. Only in case of max-stable random fields.}
  \item{param}{A list of parameter values required in the
    simulation procedure of random fields, see \bold{Examples}.}
  \item{replicates}{Numeric; a positive integer denoting the number of independent and identically distributed (iid)
    replications of a spatial or spatial-temporal random field. Optional argument, the default value is \eqn{1} then
    a single realisation is considered.}
  \item{threshold}{Numeric; a value indicating a threshold for the
    binary random field. Optional in the case that \code{model}
    is \code{BinaryGauss}, see the Section \bold{Details}.}
}

\details{
  Note that this function is also interfaced to the \bold{R} package \bold{RandomFields},
  using fast routines therein developed for the simulation of random fields, see
  for example
  \code{\link[RandomFields]{GaussRF}},
  \code{\link[RandomFields]{MaxStableRF}}, ect.
}

\value{
  Returns an object of class \code{RFsim}.
  An object of class \code{RFsim} is a list containing
  at most the following components:

  \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates;}
  \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates;}
  \item{corrmodel}{The correlation model; see \code{\link{Covmatrix}}.}
  \item{data}{The vector or matrix or array of data, see
    \code{\link{FitComposite}};}
  \item{distance}{The type of spatial distance;}
  \item{model}{The type of random field, see \code{\link{FitComposite}}.}
  \item{numcoord}{The number of spatial coordinates;}
  \item{numtime}{The number the temporal realisations of the random field;}
  \item{param}{The vector of  parameters' estimates;}
  \item{randseed}{The seed used for the random simulation;}
  \item{replicates}{The number of the iid replicatations of the random field;}
  \item{spacetime}{\code{TRUE} if spatio-temporal and \code{FALSE} if
  spatial random field;}
  \item{threshold}{The threshold for deriving the binary random field.}
}

\seealso{\code{\link{Covmatrix}},
         \code{\link[RandomFields]{GaussRF}},
         \code{\link[RandomFields]{MaxStableRF}}}

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{http://staff.deuv.cl/bevilacqua/}.}


\examples{
library(CompRandFld)
library(RandomFields)
library(mapproj)
library(fields)

################################################################
###
### Example 1. Simulation of a Gaussian random field.
### Gaussian random fields with Whittle-Matern correlation.
### One spatial replication.
###
###
###############################################################

# Define the spatial-coordinates of the points:
x <- runif(500, 0, 2)
y <- runif(500, 0, 2)

set.seed(261)
# Simulation of a spatial Gaussian random field:
data <- RFsim(x, y, corrmodel="matern", param=list(smooth=0.5,
             mean=0,sill=1,scale=0.2,nugget=0))$data

################################################################
###
### Example 2. Simulation of a binary random field based on
### the latent Gaussian random field with exponential correlation.
### One spatial replication on a regular grid
###
###
###############################################################

# Define the spatial-coordinates of the points:
x <- seq(0, 1, 0.02)
y <- seq(0, 1, 0.02)

set.seed(251)

# Simulation of a spatial binary random field:
sim <- RFsim(x, y, corrmodel="exponential", grid=TRUE,
             model="BinaryGauss", threshold=0,
             param=list(nugget=0,mean=0,scale=.1,sill=1))

image.plot(x,y,sim$data,col=terrain.colors(100))

################################################################
###
### Example 3. Simulation of a max-stable random
### extremal-t type with exponential correlation.
### One spatial replication on a regular grid
###
###
###############################################################

set.seed(341)

# Simulation of a spatial binary random field:
sim <- RFsim(x, y, corrmodel="exponential", grid=TRUE, model="ExtT",
             numblock=500, param=list(nugget=0,mean=0,scale=.1,
             sill=1,df=5))

image.plot(x,y,log(sim$data))


################################################################
###
### Example 4. Simulation of a Gaussian random field.
### with double exponential correlation.
### One spatio-temporal replication.
###
###
###############################################################

# Define the spatial-coordinates of the points:
x <- seq(0, 1, 0.1)
y <- seq(0, 1, 0.1)
# Define the temporal-coordinates:
times <- seq(1, 10, 1)

# Simulation of a spatial Gaussian random field:
sim <- RFsim(x, y, times, corrmodel="exp_exp", grid=TRUE,
             param=list(nugget=0,mean=0,scale_s=0.3,
             scale_t=0.5,sill=1))

# For a simple animation movie
# library(animation)
# saveMovie(for (i in 1:length(times)) image(x,y,sim$data[,,i]),
#           height=300,width=500,interval=1,
#           outdir = "c:/",movie.name = "spacetime_animation.mpg")




################################################################
###
### Example 5. Simulation of a Gaussian random field
### with  exponential correlation on a portion of  the earth surface
### One spatial replication.
###
###
###############################################################

lon_region<-c(-60,60)
lat_region<-c(-60,60)

lon<-seq(min(lon_region),max(lon_region),2)
lat<-seq(min(lat_region),max(lat_region),2)

data<-RFsim(coordx=lon,coordy=lat,corrmodel="exponential",
         distance="Geod",grid=TRUE,param=list(nugget=0,mean=0
         ,scale=8000,sill=1))$data

image.plot(lon,lat,data,xlab="Longitude",ylab="Latitude")
map(database="world",xlim=lon_region,ylim=lat_region,add=TRUE)
}
\keyword{Simulation}
