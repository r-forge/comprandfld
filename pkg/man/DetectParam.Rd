\name{DetectParam}
\alias{DetectParam}

\title{Identification  of the Parameters of the Correlation function}
\description{
  Subroutine called by Covariogram and others
  The procedure returns a list with the correlation model and
  the list of parameters.
}
\usage{
DetectParam(corrmodel, fixed, param)
} 
\arguments{
  \item{corrmodel}{String; the name of a correlation model.}
  \item{fixed}{A numeric vector with the fixed parameters.}
  \item{param}{A numeric vector with the parameters.}
}

\seealso{\code{\link{FitComposite}}}

\author{Simone Padoan, \email{simone.padoan@epfl.ch},
  \url{http://eflum.epfl.ch/people/padoan.en.php}.}

\keyword{Composite}
