\name{CheckLikelihood}
\alias{CheckLikelihood}

\title{Checking Composite-likelihood Type}
\description{
  Subroutine called by InitParam.
  The procedure controls the type of the composite-likelihood inserted
  by the users.
}
\usage{
CheckLikelihood(likelihood)
}
\arguments{
  \item{likelihood}{String; the configuration of the composite
    likelihood. \code{Marginal} is the default.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{http://staff.deuv.cl/bevilacqua/}.}

\keyword{Composite}
