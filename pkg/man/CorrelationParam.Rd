\name{CorrelationParam}
\alias{CorrelationParam}

\title{Lists the Parameters of a Correlation Model}
\description{
  Subroutine called by InitParam and other procedures.
  The procedure returns a list with the parameters of a given
  correlation model.
}
\usage{
CorrelationParam(corrmodel)
}
\arguments{
  \item{corrmodel}{String; the name of a correlation model. See \code{\link{Covmatrix}}} for details}


\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{http://staff.deuv.cl/bevilacqua/}.}

\examples{
require(CompRandFld)
################################################################
###
### Example 1. Parameters of the  Matern correlation model
###
###############################################################

CorrelationParam("matern")

################################################################
###
### Example 1. Parameters of the  Gneiting correlation model
###
###############################################################

CorrelationParam("gneiting")
}
\keyword{Composite}
