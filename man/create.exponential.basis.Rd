\name{create.exponential.basis}
\alias{create.exponential.basis}
\title{
  Create an Exponential Basis
}
\description{
  Create an exponential basis object defining a set of exponential
  functions with rate constants in argument ratevec.
}
\usage{
create.exponential.basis(rangeval=c(0,1), nbasis=NULL, ratevec=NULL,
                         dropind=NULL, quadvals=NULL, values=NULL,
                         basisvalues=NULL, names='exp', axes=NULL)
}
\arguments{
  \item{rangeval}{
    a vector of length 2 containing the initial and final values of the
    interval over which the functional data object can be evaluated.
  }
  \item{nbasis}{
    the number of basis functions.  Default = \code{if(is.null(ratevec))
      2 else length(ratevec)}.
  }
  \item{ratevec}{
    a vector of length \code{nbasis} of rate constants defining basis
    functions of the form \code{exp(rate*x)}.  Default = 0:(nbasis-1).
  }
  \item{dropind}{
    a vector of integers specifiying the basis functions to be dropped,
    if any.  For example, if it is required that a function be zero at
    the left boundary, this is achieved by dropping the first basis
    function, the only one that is nonzero at that point.
  }
  \item{quadvals}{
    a matrix with two columns and a number of rows equal to the number
    of quadrature points for numerical evaluation of the penalty
    integral.  The first column of \code{quadvals} contains the
    quadrature points, and the second column the quadrature weights.  A
    minimum of 5 values are required for each inter-knot interval, and
    that is often enough.  For Simpson's rule, these points are equally
    spaced, and the weights are proportional to 1, 4, 2, 4, ..., 2, 4,
    1.
  }
  \item{values}{
    a list of matrices with one row for each row of \code{quadvals} and
    one column for each basis function.  The elements of the list
    correspond to the basis functions and their derivatives evaluated at
    the quadrature points contained in the first column of
    \code{quadvals}.
  }
  \item{basisvalues}{
    A list of lists, allocated by code such as vector("list",1).  This
    field is designed to avoid evaluation of a basis system repeatedly
    at a set of argument values.  Each list within the vector
    corresponds to a specific set of argument values, and must have at
    least two components, which may be tagged as you wish.  `The first
    component in an element of the list vector contains the argument
    values.  The second component in an element of the list vector
    contains a matrix of values of the basis functions evaluated at the
    arguments in the first component.  The third and subsequent
    components, if present, contain matrices of values their derivatives
    up to a maximum derivative order.  Whenever function getbasismatrix
    is called, it checks the first list in each row to see, first, if
    the number of argument values corresponds to the size of the first
    dimension, and if this test succeeds, checks that all of the
    argument values match.  This takes time, of course, but is much
    faster than re-evaluation of the basis system.  Even this time can
    be avoided by direct retrieval of the desired array.  For example,
    you might set up a vector of argument values called "evalargs" along
    with a matrix of basis function values for these argument values
    called "basismat".  You might want too use names like "args" and
    "values", respectively for these.  You would then assign them to
    \code{basisvalues} with code such as the following:

    basisobj$basisvalues <- vector("list",1)

    basisobj$basisvalues[[1]] <- list(args=evalargs,
    values=basismat)
  }
  \item{names}{
    either a character vector of the same length as the number of basis
    functions or a simple stem used to construct such a vector.

    For \code{exponential} bases, this defaults to paste('exp',
    0:(nbasis-1), sep='').
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
  }}
\value{
  a basis object with the type \code{expon}.
}
\details{
  Exponential functions are of the type $exp(bx)$ where $b$
  is the rate constant.  If $b = 0$, the constant function is
  defined.
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
  \code{\link{basisfd}},
  \code{\link{create.bspline.basis}},
  \code{\link{create.constant.basis}},
  \code{\link{create.fourier.basis}},
  \code{\link{create.monomial.basis}},
  \code{\link{create.polygonal.basis}},
  \code{\link{create.power.basis}}
}
\examples{

#  Create an exponential basis over interval [0,5]
#  with basis functions 1, exp(-t) and exp(-5t)
basisobj <- create.exponential.basis(c(0,5),3,c(0,-1,-5))
#  plot the basis
oldpar <- par(no.readonly=TRUE)
plot(basisobj)
par(oldpar)
}
% docclass is function
\keyword{smooth}
