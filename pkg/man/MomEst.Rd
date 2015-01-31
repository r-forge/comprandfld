\name{MomEst}
\alias{MomEst}
\title{Estimation of the GEV parameters by the Method of Moments}
\description{
  Using the moment estimator, the function returns the parameter estimates
  of the generalized extreme value distribution for a given dataset of extreme values.
}
\usage{
MomEst(data, n)
}
\arguments{
  \item{data}{A vector of extreme values.}
  \item{n}{The number of observations.}
}

\value{
  The returned object is a list with:

    \item{location}{The location estimate.}
    \item{scale}{The scale estimate.}
    \item{shape}{The shape estimate.}
}

\references{
  de Haan, L. and Ferreira, A. (2006) \emph{Extreme Value Theory An
    Introduction}. Springer Verlang, New York.

}

\seealso{\code{\link{GevLogLik}}, \code{\link{FitGev}}}

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home}.}



\keyword{extremes}
