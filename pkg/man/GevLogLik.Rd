\name{GevLogLik}
\alias{GevLogLik}
\title{Log-Likelihood of the GEV Distribution}
\description{
  The function returns the log-likelihood value of the Generalized
  Extreme Value Distribution for a given set of data and parameters.
}
\usage{
GevLogLik(data, numdata, param)
}
\arguments{
  \item{data}{A vector of extreme values.}
  \item{numdata}{The number of data observations.}
  \item{param}{The vector of GEV parameters (location, scale and shape).}
}
\value{
  The log-likelihood value is returned.}

\references{
  de Haan, L. and Ferreira, A. (2006) \emph{Extreme Value Theory An
    Introduction}. Springer Verlang, New York.
}


\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
   \url{http://staff.deuv.cl/bevilacqua/}.
  }

\keyword{extremes}
