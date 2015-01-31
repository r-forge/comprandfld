\name{Kri}
\alias{Kri}

\title{Spatial and spatio temporal simple and ordinary (tapered) kriging}
\description{
  The function computes simple or ordinary (tapered) kriging, in
  addition, for a set of unknown spatial location sites and temporal
  instants and a given space or space-time covariance model, it computes
  the Kriging variance.
}
\usage{
Kri(data, coordx, coordy=NULL, coordt=NULL, corrmodel,distance="Eucl",
    grid=FALSE, loc, maxdist=NULL, maxtime=NULL, param, taper=NULL,
    tapsep=NULL, time=NULL, type="Standard",type_krig="Simple")
}
\arguments{
  \item{data}{A \eqn{d}{d}-dimensional vector (a single spatial realisation) or a (\eqn{d \times d}{d x d})-matrix (a single spatial realisation on regular grid) or a
   (\eqn{t \times d}{t x d})-matrix (a single spatial-temporal realisation)  or an (\eqn{d \times d \times t \times n }{d x d x t})-array
   (a single spatial-temporal realisation on regular grid) giving the data used for prediction.}
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of spatial sites) giving 2-dimensions of spatial coordinates or a numeric \eqn{d}{d}-dimensional vector giving
    1-dimension of spatial coordinates used for prediction.}
  \item{coordy}{A numeric vector giving 1-dimension of
    spatial coordinates used for prediction; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector or \code{grid=TRUE} otherwise it will be ignored. Optional argument, the default is \code{NULL} then \code{coordx} is expected to
    be numeric a (\eqn{d \times 2}{d x 2})-matrix.}
  \item{coordt}{A numeric vector giving 1-dimension of
    temporal coordinates used for prediction; the default is \code{NULL}
    then a spatial random field is expected.}
  \item{corrmodel}{String; the name of a correlation model, for the
    description see the Section \bold{Details}.}
   \item{distance}{String; the name of the spatial distance. The default
     is \code{Eucl}, the euclidean distance. See the Section
    \bold{Details} of \code{\link{FitComposite}}.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data used for prediction
    are interpreted as spatial or spatial-temporal realisations on a set
    of non-equispaced spatial sites (irregular grid).}
  \item{loc}{A numeric (\eqn{n \times 2}{n x 2})-matrix (where
    \code{n} is the number of spatial sites) giving 2-dimensions of spatial coordinates to be predicted.}
  \item{maxdist}{Numeric; an optional positive value indicating the maximum
    spatial compact support in the case of covariance tapering kriging.}
  \item{maxtime}{Numeric; an optional positive value indicating the maximum
    temporal compact support in the case of covariance tapering kriging.}
   \item{param}{A list of parameter values required  for the correlation model.See the Section \bold{Details}.}
   \item{taper}{String; the name of the taper correlation
    function, see the Section \bold{Details}.}
    \item{tapsep}{Numeric; an optional  value indicating the separabe
      parameter in the space time quasi taper (see \bold{Details}).}
    \item{time}{A numeric (\eqn{m \times 1}{m x 1}) vector (where
    \code{m} is the number of temporal instants) giving the temporal instants to be predicted; the default is \code{NULL}
    then only spatial prediction is performed.}
   \item{type}{String; if \code{standard}   then standard kriging is performed;if \code{Tapering}   then  kriging with covariance tapering is performed.}
   \item{type_krig}{String; the type of kriging. If \code{Simple} (the
    default) then simple kriging is performed.
    If \code{ordinary}  then ordinary kriging is performed. (See the Section \bold{Details}).}
}

\details{
  For a spatial or spatio-temporal dataset, given a set of locations and
  temporal istants and a correlation model
\code{corrmodel} with some fixed parameters, the function computes
  simple or ordinary kriging, for the specified spatial locations
  \code{loc} and temporal instants \code{time},
  providing also the respective standard error.
  For the choice of the spatial or spatio temporal correlation model see details in \code{\link{Covmatrix}} function.
  The parameter \code{param} specifies the covariance parameters, see
  \code{\link{CorrelationParam}} and \code{\link{Covmatrix}} for details.
The \code{type_krig} parameter indicates the type of kriging. In the
case of simple kriging, the known mean can be specified by the parameter
  \code{mean} within list \code{param} (See examples).
In addition, it is possible to perform kriging based on covariance
tapering for simple kriging (Furrer et. al, 2008).
In this case, space or space-time tapered function and spatial or spatio- temporal compact support
 must be specified. For the choice of a space or space-time tapered function see \code{\link{Covmatrix}}.
 When performing kriging with covariance tapering,
 sparse matrix algorithms are exploited using the package \code{spam}.}

\value{
  Returns an object of class \code{Kg}.
  An object of class \code{Kg} is a list containing
  at most the following components:
  \item{coordx}{A \eqn{d}{d}-dimensional vector of spatial coordinates used for prediction;}
  \item{coordy}{A \eqn{d}{d}-dimensional vector of spatial coordinates used for prediction;}
  \item{coordt}{A \eqn{t}{t}-dimensional vector of temporal coordinates used for prediction;}
  \item{corrmodel}{String: the correlation model;}
  \item{covmatrix}{The covariance matrix if \code{type} is \code{Standard}. An object of class spam if \code{type} is \code{Tapering}   }
  \item{data}{The vector or matrix or array of data used for prediction}
  \item{distance}{String: the type of spatial distance;}
  \item{grid}{\code{TRUE} if the spatial data used for prediction are observed in a regular grid, otherwise  \code{FALSE};}
  \item{loc}{A (\eqn{n \times 2}{n x 2})-matrix of spatial locations to be predicted.}
  \item{nozero}{In the case of tapered simple kriging the percentage of non zero values in the covariance matrix.
  Otherwise is NULL.}
  \item{numcoord}{Numeric:he number \eqn{d}{d} of spatial coordinates used for prediction;}
  \item{numloc}{Numeric: the number \eqn{n}{n} of spatial coordinates to be predicted;}
  \item{numtime}{Numeric: the number \eqn{d}{d}  of the temporal instants used for prediction;}
  \item{numt}{Numeric: the number \eqn{m}{m}  of the temporal instants to be predicted;}
  \item{param}{Numeric: The covariance parameters;}
  \item{pred}{A (\eqn{m \times n}{n x m})-matrix of spatio or spatio temporal kriging prediction;}
  \item{spacetime}{\code{TRUE} if spatio-temporal kriging and \code{FALSE} if
  spatial kriging;}
   \item{tapmod}{String: the taper model if \code{type} is \code{Tapering}. Otherwise is NULL.}
  \item{time}{A \eqn{m}{m}-dimensional vector of temporal coordinates to be predicted;}
  \item{type}{String:  the type of kriging (Standard or Tapering).}
  \item{type_krig}{String: the type of kriging (Simple or Ordinary).}
  \item{varpred}{A (\eqn{m \times n}{n x m})-matrix of spatio or spatio temporal variance kriging prediction;}

}


\references{
Padoan, S. A. and Bevilacqua, M. (2015). Analysis of Random Fields Using CompRandFld.
  \emph{Journal of Statistical Software}, \bold{63}(9), 1--27.
  
  Gaetan, C. and Guyon, X. (2010) \emph{Spatial Statistics and Modelling}.
  Spring Verlang, New York.
  Furrer R.,  Genton, M.G. and  Nychka D. (2006). \emph{Covariance
  Tapering for Interpolation of Large Spatial Datasets}.
  Journal of Computational and Graphical Statistics, \bold{15-3}, 502--523.
}

\seealso{\code{\link{Covmatrix}}}

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home}.}

\examples{

library(CompRandFld)
library(fields)

################################################################
############### Example of Spatial kriging  ####################
################################################################

# Define the spatial-coordinates of the points:
x <- runif(50, 0, 1)
y <- runif(50, 0, 1)

# Set the model's parameters:
corrmodel <- "exponential"
mean<-0
sill<-1
nugget<-0
scale<-0.5
param<-list(mean=mean,sill=sill,nugget=nugget,scale=scale)

# spatial matrix location sites
coords<-cbind(x,y)

# Simulation of the spatial Gaussian random field:
set.seed(3132)
data <- RFsim(coordx=coords, corrmodel=corrmodel,
              param=param)$data
start<-list(scale=scale,sill=sill)
fixed<-list(mean=mean,nugget=nugget)
# Maximum likelihood fitting :
fit <- FitComposite(data, coordx=coords, corrmodel=corrmodel,
                    likelihood='Full', type='Standard',
                    start=start,fixed=fixed)

# locations to predict
xx<-seq(0,1,0.02)
loc_to_pred<-as.matrix(expand.grid(xx,xx))

################################################################
###
### Example 1. Spatial simple kriging of n sites of a
### Gaussian random fields with exponential correlation.
###
###############################################################

pr<-Kri(loc=loc_to_pred,coordx=coords,corrmodel=corrmodel,
       param= as.list(c(fit$param,fit$fixed)), data=data)

################################################################
###
### Example 2. Spatial tapered simple kriging of n sites of a
### Gaussian random fields with exponential correlation.
###
###############################################################

##pr_tap<-Kri(loc=loc_to_pred,coordx=coords,corrmodel=corrmodel,data=data,
##       param= as.list(c(fit$param,fit$fixed)),type="Tapering",
##       maxdist=0.15,taper="Wendland1")

##colour <- rainbow(100)

##par(mfrow=c(2,2))
# simple kriging map prediction
##image.plot(xx, xx, matrix(pr$pred,ncol=length(xx)),col=colour,
     ##      xlab="",ylab="",main="Simple Kriging")

# simple kriging map prediction variance
##image.plot(xx, xx, matrix(pr$varpred,ncol=length(xx)),col=colour,
   ##        xlab="",ylab="",main="Std error")

# simple tapered kriging map prediction
##image.plot(xx, xx, matrix(pr_tap$pred,ncol=length(xx)),col=colour,
  ##         xlab="",ylab="",main="Simple Tapered Kriging")

# simple taperd kriging map prediction variance
##image.plot(xx, xx, matrix(pr_tap$varpred,ncol=length(xx)),col=colour,
 ##          xlab="",ylab="",main="Std error")


################################################################
########### Examples of Spatio-temporal kriging  ###############
################################################################

# Define the spatial-coordinates of the points:
x <- runif(15, 0, 1)
y <- runif(15, 0, 1)
coords<-cbind(x,y)
times<-1:7

#  Define the times to predict
times_to_pred<-8:10

# Define model correlation and associated parameters
corrmodel<-"exp_exp"
param<-list(nugget=0,mean=1,scale_s=1,scale_t=2,sill=2)

# Simulation of the space time Gaussian random field:
set.seed(31)
data<-RFsim(coordx=coords,coordt=times,corrmodel=corrmodel,
           param=param)$data

# Maximum likelihood fitting of the space time random field:
start <- list(scale_s=1,scale_t=2,sill=2)
fixed <- list(mean=1,nugget=0)
fit <- FitComposite(data, coordx=coords, coordt=times,
                    corrmodel=corrmodel, likelihood='Marginal',
                    type='Pairwise',start=start,fixed=fixed,
                    maxdist=0.5,maxtime=3)



################################################################
###
### Example 3. Spatio temporal simple kriging of n locations
### sites and m temporal instants for a Gaussian random fields
### with estimated double exponential correlation.
###
###############################################################

param<-as.list(c(fit$param,fit$fixed))

pr<-Kri(loc=loc_to_pred,time=times_to_pred,coordx=coords,coordt=times,
       corrmodel=corrmodel, param=param,data=data)

par(mfrow=c(3,2))

colour <- rainbow(100)

for(i in 1:3){
image.plot(xx, xx, matrix(pr$pred[i,],ncol=length(xx)),col=colour,
           main = paste("Kriging Time=" , i),ylab="")
image.plot(xx, xx, matrix(pr$varpred[i,],ncol=length(xx)),col=colour,
           main = paste("Std error Time=" , i),ylab="")
}



################################################################
###
### Example 4. Spatio temporal tapered simple kriging of n locations
### sites and m temporal instants for a Gaussian random fields
### with estimated double exponential correlation.
###
###############################################################

#pr_tap<-Kri(loc=loc_to_pred,time=times_to_pred,coordx=coords,coordt=times,
#       corrmodel=corrmodel, param=param,type="Tapering",maxdist=0.4,maxtime=4,
#       taper="Wendland2_Wendland2",data=data)

#par(mfrow=c(3,2))

#for(i in 1:3){
#image.plot(xx, xx, matrix(pr_tap$pred[i,],ncol=length(xx)),col=colour,
#           main = paste("Tapered Kriging Time=" , i),ylab="")
#image.plot(xx, xx, matrix(pr_tap$varpred[i,],ncol=length(xx)),col=colour,
#           main = paste("Tapered Std error Time=" , i),ylab="")
#}
}

\keyword{Composite}
