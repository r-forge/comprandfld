\name{CompLikelihood}
\alias{CompLikelihood}

\title{Optimizes the Composite log-likelihood}
\description{
  Subroutine called by FitComposite.
  The procedure estimates the model parameters by maximisation of the
  composite log-likelihood.
}
\usage{
CompLikelihood(coordx, coordy, corrmodel, data, distance, 
               flagcorr, flagnuis, fixed, grid, likelihood, 
               lower, model, namescorr, namesnuis, namesparam,
               numparam, numparamcorr, optimizer, param,
               spacetime, threshold, type, upper, varest,
               vartype, winconst, winstp)
}
\arguments{
  \item{coordx}{A numeric (\eqn{d \times 2}{d x 2})-matrix (where
    \code{d} is the number of points) assigning 2-dimensions of coordinates or a numeric vector assigning
    1-dimension of coordinates.}
  \item{coordy}{A numeric vector assigning 1-dimension of
    coordinates; \code{coordy} is interpreted only if \code{coordx} is a numeric
    vector otherwise it will be ignored.}
  \item{corrmodel}{Numeric; the id of the correlation model.}
  \item{data}{A numeric vector or a (\eqn{n \times d}{n x d})-matrix or (\eqn{d \times d
      \times n}{d x d x n})-matrix of observations.}
  \item{distance}{String; the name of the spatial distance. The default is \code{Eucl},
            the euclidean distance. See the Section \bold{Details}.}
  \item{flagcorr}{A numeric vector of binary values denoting which
    paramerters of the correlation function will be estimated.}
  \item{flagnuis}{A numeric vector of binary values denoting which
    nuisance paramerters will be estimated.}
  \item{fixed}{A numeric vector of parameters that
    will be considered as known values.}
  \item{grid}{Logical; if \code{FALSE} (the default) the data
    are interpreted as a vector or a (\eqn{n \times d}{n x d})-matrix,
    instead if \code{TRUE} then (\eqn{d \times d
      \times n}{d x d x n})-matrix is considered.}
  \item{likelihood}{String; the configuration of the
    compositelikelihood, see \code{\link{FitComposite}}.}
  \item{lower}{A numeric vector with the lower bounds of the parameters' ranges.}
  \item{model}{Numeric; the id value of the density associated to the likelihood
    objects.}
  \item{namescorr}{String; the names of the correlation parameters.}
  \item{namesnuis}{String; the names of the nuisance parameters.}
  \item{namesparam}{String; the names of the parameters to be maximised.}
  \item{numparam}{Numeric; the number of parameters to be maximised.}
  \item{numparamcorr}{Numeric; the number of correlation parameters.}
  \item{optimizer}{String; the optimization algorithm
    (see \code{\link{optim}} for details). 'Nelder-Mead' is the
    default.}
  \item{param}{A numeric vector of parameters' values.}
  \item{spacetime}{Logical; if \code{TRUE} the random field is
    spatial-temporal otherwise is a spatial field.}
  \item{threshold}{Numeric; a value indicating a threshold for the
    binary random field, see \code{\link{FitComposite}}.}
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
  \item{upper}{A numeric vector with the upper bounds of the parameters'
    ranges.}
  \item{varest}{Logical; if \code{TRUE} the estimate' variances and
    standard errors are returned.
    \code{FALSE} is the default.}
  \item{vartype}{String; the type of estimation method for computing
  the estimate variances, see \code{\link{FitComposite}}.}
  \item{winconst}{Numeric; a positive value -- if
    \code{vartype=SubSamp} -- determines the window size in the sub-sampling
    estimates of the variances, see \code{\link{FitComposite}}.}
  \item{winstp}{Numeric; a positive value -- if
    \code{vartype=SubSamp} -- determines the window step in the sub-sampling
    estimates of the variances, see \code{\link{FitComposite}}.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home}.}


\keyword{Composite}
