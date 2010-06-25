\name{Wls}
\alias{Wls}

\title{Computation of the Weighted Least Squares}
\description{
  Subroutine called by WLeastSquare. The procedure
  computes the weighted least squares for a given
  set of data and parameters.
}
\usage{
Wls(bins, corrmodel, fixed, lenbins, moments, numbins, param, weighted)
} 
\arguments{
  \item{bins}{A numeric vector with adjacent intervals of grouped distances}
  \item{corrmodel}{Numeric; the id of the correlation mode.}
  \item{fixed}{A named list giving the values of the parameters that
    will be considered as known values.}
  \item{lenbins}{A numeric vector with the number of observations that
    fall in each bin.}
  \item{moments}{A vector with the estimated means.}
  \item{numbins}{Numeric; the number of bins.}
  \item{param}{A numeric vector of parameters' values.}
  \item{weighted}{{Logical; if \code{TRUE} the least squares are
    weighted.}}
}

\seealso{\code{\link{WLeastSquare}}, \code{\link{FitComposite}}}


\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}.}


\keyword{Composite}
