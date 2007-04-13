\name{plot.basisfd}
\alias{plot.basisfd}
\title{
  Plot a Basis Object
}
\description{
  Plots all the basis functions.
}
\usage{
plot.basisfd(x, ...)
}
\arguments{
  \item{x}{
    a basis object
  }
  \item{\dots }{
    additional plotting parameters that can be used
    in \code{plot}
  }
}
\value{
  none
}
\section{Side Effects}{
  a plot of the basis functions
}
\seealso{
\code{\link{plot.fd}}
}
\examples{

# set up the b-spline basis for the lip data, using 23 basis functions,
#   order 4 (cubic), and equally spaced knots.
#  There will be 23 - 4 = 19 interior knots at 0.05, ..., 0.95
lipbasis <- create.bspline.basis(c(0,1), 23)
# plot the basis functions
plot(lipbasis)

}
\keyword{}
% docclass is function
\keyword{smooth}
