\name{FitGev}
\alias{FitGev}
\title{Maximum-likelihood Fitting of the Generalized Extreme Value
  Distribution}
\description{
  the function returns the parameters' estimates and the variances of
  the estimates (if required) of the generalized extreme value
  distribution for a given dataset of extreme values.
}
\usage{
FitGev(data, method='Nelder-Mead', start, varest=FALSE)
} 
\arguments{
  \item{data}{A vector of extreme values.}
  \item{method}{The optimization method (see \code{\link{optim}} for
    details). 'Nelder-Mead' is the default.}
  \item{start}{A named list with the initial values for the
    parameters over which the likelihood is to be maximized.}
  \item{varest}{Logical; if \code{TRUE} the estimate' variances and
    the standard errors are returned, instead if \code{FALSE} (the
    default) only the estimate are computed.}
}

\details{
If \code{start} is omitted the routine is computing the
    starting values using moment estimators.
}

\value{
  The returned object is a list with:
  
    \item{param}{The vector of  parameters' estimates.}
    \item{varcov}{The matrix of the variance-covariance of the
      estimates.}
    \item{stderr}{The vector of the standard errors.}
}

\references{
  Smith, R. L. (1985)
  Maximum likelihood estimation in a class of non-regular cases.
  \emph{Biometrika}, \bold{72}, 67--90.

  de Haan, L. and Ferreira, A. (2006) \emph{Extreme Value Theory An
    Introduction}. Springer Verlang, New York.
  
}

\seealso{\code{\link{GevLogLik}}, \code{\link{optim}}}


\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://people.epfl.ch/simone.padoan}}


\keyword{extremes}
