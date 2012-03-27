\name{CheckCorrModel}
\alias{CheckCorrModel}

\title{Checking Correlation Model}
\description{
  Subroutine called by InitParam.
  The procedure controls if the correlation model
  inserted is correct.
}
\usage{
CheckCorrModel(corrmodel)
}
\arguments{
  \item{corrmodel}{String; the name of a correlation model, for the
    description see \code{\link[RandomFields]{CovarianceFct}}.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@stat.unipd.it},
  \url{http://homes.stat.unipd.it/padoan};
  Moreno Bevilacqua, \email{moreno.bevilacqua@unibg.it},
   \url{www.unibg.it/pers/?moreno.bevilacqua}.
  }

\keyword{Composite}
