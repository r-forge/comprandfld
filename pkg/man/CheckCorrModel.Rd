\name{CheckCorrModel}
\alias{CheckCorrModel}

\title{Check of the Correlation Model}
\description{
  Subroutine called by InitParam.
  The procedure controls if the correlation model
  inserted has been implemented.
}
\usage{
CheckCorrModel(corrmodel)
} 
\arguments{
  \item{corrmodel}{String; the name of a correlation model, for the
    description see \code{\link[RandomFields]{CovarianceFct}}.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://people.epfl.ch/simone.padoan}.}

\keyword{Composite}
