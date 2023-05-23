\name{create.polygonal.basis}
\alias{create.polygonal.basis}
\title{
  Create a Polygonal Basis
}
\description{
  A basis is set up for constructing polygonal lines, consisting of
  straight line segments that join together.
}
\usage{
create.polygonal.basis(rangeval=NULL, argvals=NULL, dropind=NULL,
        quadvals=NULL, values=NULL, basisvalues=NULL, names='polygon',
        axes=NULL)
}
\arguments{
  \item{rangeval}{
    a numeric vector of length 2 defining the interval over which the
    functional data object can be evaluated;  default value is
    \code{if(is.null(argvals)) 0:1 else range(argvals)}.

    If \code{length(rangeval) == 1} and \code{rangeval <= 0}, this is an
    error.  Otherwise, if \code{length(rangeval) == 1}, \code{rangeval}
    is replaced by \code{c(0,rangeval)}.

    If length(rangeval)>2 and \code{argvals} is not provided, this extra
    long \code{rangeval} argument is assigned to \code{argvals}, and
    then \code{rangeval = range(argvale)}.
  }
  \item{argvals}{
    a strictly increasing vector of argument values at which line
    segments join to form a polygonal line.
  }
  \item{dropind}{
    a vector of integers specifiying the basis functions to
    be dropped, if any.  For example, if it is required that
    a function be zero at the left boundary, this is achieved
    by dropping the first basis function, the only one that
    is nonzero at that point.
  }
  \item{quadvals}{
    a matrix with two columns and a number of rows equal to the number
    of quadrature points for numerical evaluation of the penalty
    integral.  The first column of \code{quadvals} contains the
    quadrature points, and the second column the quadrature weights.  A
    minimum of 5 values are required for each inter-knot interval, and
    that is often enough.  For Simpson's rule, these points are equally
    spaced, and the weights are proportional to These are proportional
    to 1, 4, 2, 4, ..., 2, 4, 1.
  }
  \item{values}{
    a list containing the basis functions and their derivatives
    evaluated at the quadrature points contained in the first
    column of \code{ quadvals }.
  }
  \item{basisvalues}{
    A list of lists, allocated by code such as vector("list",1).  This
    is designed to avoid evaluation of a basis system repeatedly
    at a set of argument values.  Each sublist corresponds to a specific
    set of argument values, and must have at least two components, which
    may be named as you wish.  The first component of a sublist contains
    the argument values.  The second component contains a matrix of
    values of the basis functions evaluated at the arguments in the
    first component.  The third and subsequent components, if present,
    contain matrices of values their derivatives up to a maximum
    derivative order.  Whenever function \code{getbasismatrix} is
    called, it checks the first list in each row to see, first, if the
    number of argument values corresponds to the size of the first
    dimension, and if this test succeeds, checks that all of the
    argument values match.  This takes time, of course, but is much
    faster than re-evaluation of the basis system.  Even this time can
    be avoided by direct retrieval of the desired array. For example,
    you might set up a vector of argument values called "evalargs" along
    with a matrix of basis function values for these argument values
    called "basismat".  You might want too use tags like "args" and
    "values", respectively for these.  You would then assign them to
    \code{basisvalues} with code such as the following:

    basisobj$basisvalues <- vector("list",1)

    basisobj$basisvalues[[1]] <- list(args=evalargs, values=basismat)
  }
  \item{names}{
    either a character vector of the same length as the number of basis
    functions or a single character string to which \code{1:nbasis} are
    appended as \code{paste(names, 1:nbasis, sep=''}.  For example, if
    \code{nbasis = 4}, this defaults to \code{c('polygon1', 'polygon2',
      'polygon3', 'polygon4')}.
  }
  \item{axes}{
    an optional list used by selected \code{plot} functions to create
    custom \code{axes}.  If this \code{axes} argument is not
    \code{NULL}, functions \code{plot.basisfd}, \code{plot.fd},
    \code{plot.fdSmooth} \code{plotfit.fd}, \code{plotfit.fdSmooth}, and
    \code{plot.Lfd} will create axes via \code{x$axes[[1]]} and
      \code{x$axes[-1]}.  The primary example of this uses
    \code{list("axesIntervals", ...)}, e.g., with \code{Fourier} bases
    to create \code{CanadianWeather} plots
  }
}
\value{
  a basis object with the type \code{polyg}.
}
\details{
  The actual basis functions consist of triangles, each with its apex
  over an argument value. Note that in effect the polygonal basis is
  identical to a B-spline basis of order 2 and a knot or break value at
  each argument value.  The range of the polygonal basis is set to the
  interval defined by the smallest and largest argument values.
}
\seealso{
\code{\link{basisfd}},
\code{\link{create.bspline.basis}},
\code{\link{create.constant.basis}},
\code{\link{create.exponential.basis}},
\code{\link{create.fourier.basis}},
\code{\link{create.monomial.basis}},
\code{\link{create.power.basis}}
}
\examples{
#  Create a polygonal basis over the interval [0,1]
#  with break points at 0, 0.1, ..., 0.95, 1
(basisobj <- create.polygonal.basis(seq(0,1,0.1)))
#  plot the basis
oldpar <- par(no.readonly=TRUE)
plot(basisobj)
par(oldpar)
}
% docclass is function
\keyword{smooth}
