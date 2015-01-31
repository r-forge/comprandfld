\name{CheckType}
\alias{CheckType}

\title{Checking Likelihood Objects}
\description{
  Subroutine called by InitParam.
  The procedure controls the type of likelihood objects inserted by the users.
}
\usage{
CheckType(type)
}
\arguments{
  \item{type}{String; the type of the likelihood objects. If \code{Pairwise} (the
    default) then the marginal composite likelihood is formed by
    pairwise marginal likelihoods.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home}.}

\keyword{Composite}
