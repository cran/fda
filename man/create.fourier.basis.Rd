\name{create.fourier.basis}
\alias{create.fourier.basis}
\title{
  Create a Fourier Basis
}
\description{
  Create an Fourier basis object defining a set of Fourier
  functions with specified period.
}
\usage{
create.fourier.basis(rangeval=c(0, 1), nbasis=3,
              period=diff(rangeval), dropind=NULL, quadvals=NULL,
              values=NULL, basisvalues=NULL, names=NULL,
              axes=NULL)
}
\arguments{
  \item{rangeval}{
    a vector of length 2 containing the initial and final values of the
    interval over which the functional data object can be evaluated.
  }
  \item{nbasis}{
    positive odd integer:  If an even number is specified, it is rounded
    up to the nearest odd integer to preserve the pairing of sine and
    cosine functions.  An even number of basis functions  only makes
    sense when there are always only an even number of observations at
    equally spaced points;  that case can be accommodated using dropind =
    nbasis-1 (because the bases are \code{const}, \code{sin},
    \code{cos}, ...).
  }
  \item{period}{
    the width of any interval over which the Fourier functions repeat
    themselves or are periodic.
  }
  \item{dropind}{
    an optional vector of integers specifiying basis functions to be
    dropped.
  }
  \item{quadvals}{
    an optional matrix with two columns and a number of rows equal to
    the number of quadrature points for numerical evaluation of the
    penalty integral.  The first column of \code{quadvals} contains the
    quadrature points, and the second column the quadrature weights.  A
    minimum of 5 values are required for each inter-knot interval, and
    that is often enough.  For Simpson's rule, these points are equally
    spaced, and the weights are proportional to 1, 4, 2, 4, ..., 2, 4,
    1.
  }
  \item{values}{
    an optional list of matrices with one row for each row of
    \code{quadvals} and one column for each basis function.  The
    elements of the list correspond to the basis functions and their
    derivatives evaluated at the quadrature points contained in the
    first column of \code{quadvals}.
  }
  \item{basisvalues}{
    an optional list of lists, allocated by code such as
    vector("list",1).  This field is designed to avoid evaluation of a
    basis system repeatedly at a set of argument values.  Each sublist
    corresponds to a specific set of argument values, and must have at
    least two components:  a vector of argument values and a matrix of
    the values the basis functions evaluated at the arguments in the
    first component.  Third and subsequent components, if present,
    contain matrices of values their derivatives.  Whenever function
    getbasismatrix is called, it checks the first list in each row to
    see, first, if the number of argument values corresponds to the size
    of the first dimension, and if this test succeeds, checks that all
    of the argument values match.  This takes time, of course, but is
    much faster than re-evaluation of the basis system.  Even this time
    can be avoided by direct retrieval of the desired array.  For
    example, you might set up a vector of argument values called
    "evalargs" along with a matrix of basis function values for these
    argument values called "basismat".  You might want too use tags like
    "args" and "values", respectively for these.  You would then assign
    them to \code{basisvalues} with code such as the following:

	basisobj\$basisvalues <- vector("list",1)

	basisobj\$basisvalues[[1]] <- list(args=evalargs,
	values=basismat)
  }
  \item{names}{
    either a character vector of the same length as the number of basis
    functions or a simple stem used to construct such a vector.

    If \code{nbasis} = 3, \code{names} defaults to c('const', 'cos',
    'sin').  If \code{nbasis} > 3, \code{names} defaults to c('const',
    outer(c('cos', 'sin'), 1:((nbasis-1)/2), paste, sep='')).

    If names = NA, no names are used.
  }
  \item{axes}{
    an optional list used by selected \code{plot} functions to create
    custom \code{axes}.  If this \code{axes} argument is not
    \code{NULL}, functions \code{plot.basisfd}, \code{plot.fd},
    \code{plot.fdSmooth} \code{plotfit.fd}, \code{plotfit.fdSmooth}, and
    \code{plot.Lfd} will create axes via \code{x$axes[[1]]} and
      \code{x$axes[-1]}.  The primary example of this is to create
    \code{CanadianWeather} plots using \code{list("axesIntervals")}
  }
}
\value{
  a basis object with the type \code{fourier}.
}
\details{
  Functional data objects are constructed by specifying a set of basis
  functions and a set of coefficients defining a linear combination of
  these basis functions.  The Fourier basis is a system
  that is usually used for periodic functions.  It has the advantages
  of very fast computation and great flexibility.   If the data are
  considered to be nonperiod, the Fourier basis is usually preferred.
  The first Fourier basis function is the constant function.  The
  remainder are sine and cosine pairs with integer multiples of the
  base period. The number of basis functions generated is always odd.
}
\seealso{
  \code{\link{basisfd}},
  \code{\link{create.bspline.basis}},
  \code{\link{create.constant.basis}},
  \code{\link{create.exponential.basis}},
  \code{\link{create.monomial.basis}},
  \code{\link{create.polygonal.basis}},
  \code{\link{create.power.basis}}
}
\examples{
# Create a minimal Fourier basis for annual data
#  using 3 basis functions
yearbasis3 <- create.fourier.basis(c(0,365),
                    axes=list("axesIntervals") )
#  plot the basis
plot(yearbasis3)

# Identify the months with letters
plot(yearbasis3, axes=list('axesIntervals', labels=monthLetters))

# The same labels as part of the basis object
yearbasis3. <- create.fourier.basis(c(0,365),
       axes=list("axesIntervals", labels=monthLetters) )
plot(yearbasis3.)

# set up the Fourier basis for the monthly temperature data,
#  using 9 basis functions with period 12 months.
monthbasis <- create.fourier.basis(c(0,12), 9, 12.0)

#  plot the basis
plot(monthbasis)

# Create a false Fourier basis using 1 basis function.
falseFourierBasis <- create.fourier.basis(nbasis=1)
#  plot the basis:  constant
plot(falseFourierBasis)

}
% docclass is function
\keyword{smooth}
