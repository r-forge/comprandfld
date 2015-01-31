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

\author{Simone Padoan, \email{simone.padoan@unibocconi.it},
  \url{http://faculty.unibocconi.it/simonepadoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@uv.cl},
  \url{https://sites.google.com/a/uv.cl/moreno-bevilacqua/home}.}

\keyword{Composite}
