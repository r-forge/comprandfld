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

\author{Simone Padoan, \email{simone.padoan@stat.unipd.it},
  \url{http://homes.stat.unipd.it/padoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@unibg.it},
  \url{www.unibg.it/pers/?moreno.bevilacqua}.}


\keyword{Composite}
