\name{SetRangeParam}
\alias{SetRangeParam}

\title{Identification of the Parameters' Range}
\description{
  Subroutine called by InitParam and the fitting procedures.
  The procedure returns the range of the parameters for a given
  vector of parameters.
}
\usage{
SetRangeParam(namesparam, numparam)
} 
\arguments{
  \item{namesparam}{String; the names of the parameters.}
  \item{numparam}{Numeric; the numer of parameters.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}.}

\keyword{Composite}
