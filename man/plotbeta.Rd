\name{plotbeta}
\alias{plotbeta}
\title{
  Plot a functional parameter object with confidence limits
}
\description{
  Plot a functional parameter object with confidence limits
}
\usage{
plotbeta(betaestlist, betastderrlist, argvals=NULL, xlab="", ...)
}
\arguments{
  \item{betaestlist}{
    a list containing one or more functional parameter objects
    (class = fdPar) or functional data objects (class = fd).
  }
  \item{betastderrlist}{
    a list containing functional data objects for the standard errors of
    the objects in \code{betaestlist}.
  }
  \item{argvals}{
    a sequence of values at which to evaluate \code{betaestlist} and
    \code{betastderrlist}.
  }
  \item{xlab}{ x axis label }
  \item{\dots }{
    additional plotting parameters passed to \code{plot}.
  }
}
\value{
  none
}
\section{Side Effects}{
  a plot of the basis functions
}
\references{
  Ramsay, James O., Hooker, Giles, and Graves, Spencer (2009),
    \emph{Functional data analysis with R and Matlab}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2005), 
    \emph{Functional Data Analysis, 2nd ed.}, Springer, New York.
  
  Ramsay, James O., and Silverman, Bernard W. (2002), 
    \emph{Applied Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{plot.fd}}
}
\keyword{smooth}
