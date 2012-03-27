\name{CheckModel}
\alias{CheckModel}

\title{Checking Random Field type}
\description{
  Subroutine called by InitParam.
  The procedure controls the type of random field inserted by the users.
}
\usage{
CheckModel(model)
}
\arguments{
   \item{model}{String; the density associated to the likelihood
    objects. \code{Gaussian} is the default.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@stat.unipd.it},
  \url{http://homes.stat.unipd.it/padoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@unibg.it},
  \url{www.unibg.it/pers/?moreno.bevilacqua}.}

\keyword{Composite}
