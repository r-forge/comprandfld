\name{Dist2Dist}
\alias{Dist2Dist}
\title{Switches from an EV to Another EV Distribution}
\description{
  The function transforms observations belonging to the GEV class
  from one model to another.
}
\usage{
Dist2Dist(data, from='Gev', to='sFrechet', loc=NULL, scale=NULL, shape=NULL)
}
\arguments{
  \item{data}{A numeric vector or a matrix of extreme values.}
  \item{from}{The name of the original extreme value distribution,
    i.e. \code{Gev} (the default), see the \bold{Details} section.}
  \item{to}{The name of the desired extreme value distribution,
    i.e. \code{sFrechet} (the default),
    see the \bold{Details} section.}
  \item{loc}{A numeric value or vector of location parameters.}
  \item{scale}{A numeric value or vector of scale parameters.}
  \item{shape}{A numeric value or vector of shape parameters.}
}

\details{
  If \code{data} is a numeric vector of length \code{n} then the dataset is consider as a
  realisation from an univariate extreme value distribution. Instead, if
  \code{data} is a (\eqn{n \times d}{n x d})-matrix then the columns
  represent the different variables with extreme value distributions
  and the rows represent the iid replications. Finally,
  if \code{data} is a (\eqn{d \times d \times n}{d x d x n})-matrix then
  the columns and rows represent the different variables and the third
  dimension represents the iid replications.

  The parameters \code{from} and \code{to} indicate the original extreme
  value distribution(s) from which the observations are drawn and the
  target extreme value distribution(s) that the transformed data will
  follow. The options are:
  \enumerate{
    \item \code{from=Gev} (generalised extreme value distribution):
    \itemize{
      \item \code{to=Uniform}, which means uniform distribution;
      \item \code{to=sFrechet}, which means standard (or unit) Frechet distribution,
  that is GEV(1,1,1);
      \item \code{to=sGumbel}, which means standard Gumbel distribution, that is GEV(0,1,1);
      \item \code{to=sWeibull}, which means standard Weibull distribution, that is GEV(1,1,-1);
      \item \code{to=Gev}, which means generalised extreme value
  distribution. Note, that in this case, it is required to insert vectors of
  location, scale and shape parameters with dimension \code{n} in the
  univariate case, dimension \code{d} when \code{data} is (\eqn{n \times d}{n
    x d})-matrix and dimension \eqn{n \times d}{n x d} when \code{data} is
  (\eqn{d \times d \times n}{d x d x n})-matrix.
    }
    \item \code{from=sFrechet}
    \itemize{ \item \code{to=Gev}.}
    \item \code{from=sGumbel}
    \itemize{ \item \code{to=Gev}.}
    \item \code{from=sWeibull}
    \itemize{ \item \code{to=Gev}.}
  }
}

\value{
  A numeric vector or matrix of transformed values following
  the desired distribution.
}

\references{
  de Haan, L. and Ferreira, A. (2006) \emph{Extreme Value Theory An
    Introduction}. Springer Verlang, New York.

}

\seealso{\code{\link{FitGev}}}

\author{Simone Padoan, \email{simone.padoan@stat.unipd.it},
  \url{http://homes.stat.unipd.it/padoan}.}

\keyword{extremes}
