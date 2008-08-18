\name{as.function.fd}
\alias{as.function.fd}
\alias{as.function.basisfd}
\alias{as.function.fdSmooth}
\title{
  Convert a functional data or functional basis object to a function
}

\description{
  Create a function to evaluate a functional data or basis object for
  specific argument values.  
}
\usage{
\method{as.function}{basisfd}(x, variable='x', \dots)
\method{as.function}{fd}(x, variable='x', \dots)
\method{as.function}{fdSmooth}(x, variable='x', \dots)
}
\arguments{
  \item{x}{
    an object of class \code{basisfd}, \code{fd}, or \code{fdSmooth}.
  }
  \item{variable}{
    A character string giving the name to be used for the formal
    argument of the resulting function.
  }
  \item{\dots}{Not currently used. }
}
\details{
  Create a function calling eval.basisfd(variable, x), eval.fd(variable,
  x) or eval.fdSmooth(variable, x).  
}
\value{
  An \code{R} function implementing the evaluation.
}
\author{Spencer Graves}
\examples{
bspl1.1 <- create.bspline.basis(norder=1, breaks=0:1)
bspl1.1f <- as.function(bspl1.1)

}
\keyword{symbolmath }

