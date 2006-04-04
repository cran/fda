\name{plotfit.fd}
\alias{plotfit.fd}
\title{
  Plot a Functional Data Object With Data
}
\description{
Functional observations are plotted one by one along with the data
used to estimate them, or residuals can be plotted.  This function
is useful for assessing how well a functional data object fits
the actual discrete data.
}
\usage{
plotfit.fd(x, argvals, fdobj, rng = rangeval,
           index = 1:nrep, nfine = 101, residual = FALSE,
           sortwrd = FALSE, titles=NULL, ...)
}
\arguments{
\item{x}{
a vector, matrix or array containing the discrete observations used to
estimate the functional data object.
}
\item{argvals}{
a vector containing the argument values corresponding to the
first dimension of \code{y}.
}
\item{fdobj}{
a functional data object estimated from the data.
}
\item{rng}{
a vector of length 2 specifying a reduced range contained within
the range defined by the functional data object or the arguments.
This is used to examine the fit over a limited segment of the
whole range of it.
}
\item{index}{
a set of indices of functions if only a subset of the
observations are to be plotted.
}
\item{nfine}{
the number of argument values used to define the plot of the
functional data object.  This may need to be increased if
the functions have a great deal of fine detail.
}
\item{residual}{
a logical variable:
if \code{TRUE}, the residuals are plotted rather than the
data and functional data object.
}
\item{sortwrd}{
a logical variable:
if \code{TRUE}, the observations are sorted for
plotting by the size of the sum of squared residuals.
}
\item{titles}{
a vector containing strings that are titles for each observation.
}
\item{\dots }{
additional arguments such as axis labels that may be used
with other \code{plot} functions.
}
}
\seealso{
plot.fd, lines.fd
}
% docclass is function
\keyword{smooth}
