\name{basisfd}
\alias{basisfd}
\alias{[.basisfd}
\title{
  Define a Functional Basis Object
}
\description{
  This is the constructor function for objects of the \code{basisfd}
  class.  Each function that sets up an object of this class must call
  this function.  This includes functions \code{create.bspline.basis},
  \code{create.constant.basis}, \code{create.fourier.basis}, and so
  forth that set up basis objects of a specific type.  Ordinarily, user
  of the functional data analysis software will not need to call this
  function directly, but these notes are valuable to understanding what
  the "slots" or "members" of the \code{basisfd} class are.
}
\usage{
basisfd(type, rangeval, nbasis, params,
        dropind=vector('list', 0),
        quadvals=vector('list', 0),
        values=vector("list", 0),
        basisvalues=vector('list', 0))
}
\arguments{
  \item{type}{
    a character string indicating the type of basis.  Currently,
    there are eight possible types:
    \describe{
      \item{Bspline, bspline, Bsp, bsp}{b-spline basis}
      \item{const, con, constant}{constant basis}
      \item{exp, expon, exponen, exponential}{exponential basis}
      \item{Fourier, fourier, Fou, fou}{Fourier basis}
      \item{mon, monom, monomial}{monomial basis}
      \item{polyg, polygon, polygonal}{polygonal basis}
      \item{polynomial, poly}{polynomial basis}
      \item{power, pow}{power basis}
    }
  }
  \item{rangeval}{
    a vector of length 2 containing the lower and upper boundaries
    of the range over which the basis is defined
  }
  \item{nbasis}{
    the number of basis functions
  }
  \item{params}{
    a vector of parameter values defining the basis.

    If the basis is "fourier", this is a single number indicating the
    period.  That is, the basis functions are periodic on the interval
    (0,PARAMS) or any translation of it.

    If the basis is "bspline", the values are interior points at which
    the piecewise polynomials join.  Note that the number of basis
    functions NBASIS is equal to the order of the Bspline functions plus
    the number of interior knots, that is the length of PARAMS.  This
    means that NBASIS must be at least 1 larger than the length of PARAMS.
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
    of argument values used to approximate an integral using Simpson's
    rule.  The first column contains these argument values.  A minimum
    of 5 values are required for each inter-knot interval, and that is
    often enough. These are equally spaced between two adjacent knots.
    The second column contains the weights used for Simpson's rule.
    These are proportional to 1, 4, 2, 4, ..., 2, 4, 1.
  }
  \item{values}{
    a list, with entries containing the values of the basis function
    derivatives starting with 0 and going up to the highest derivative
    needed.  The values correspond to quadrature points in
    \code{quadvals} and it is up to the user to decide whether or not to
    multiply the derivative values by the square roots of the quadrature
    weights so as to make numerical integration a simple matrix
    multiplication.  Values are checked against \code{quadvals} to
    ensure the correct number of rows, and against \code{nbasis} to
    ensure the correct number of columns.

    \code{values} contains values of basis functions and derivatives at
    quadrature points weighted by square root of quadrature weights.
    These values are only generated as required, and only if the
    \code{quadvals} is not matrix("numeric",0,0).
    %    a list containing the basis functions and their derivatives
    %    evaluated at the quadrature points contained in the first
    %    column of \code{ quadvals }.
  }
  \item{basisvalues}{
    a list of lists.  This is designed to avoid evaluation of a basis
    system repeatedly at a set of argument values. Each sublist
    corresponds to a specific set of argument values, and must have at
    least two components, which may be named as you wish.  The first
    component in an element of the list vector contains the argument
    values.  The second component is a matrix of values of the basis
    functions evaluated at the arguments in the first component.
    Subsequent components, if present, are matrices of values their
    derivatives up to a maximum derivative order.  Whenever function
    \code{getbasismatrix} is called, it checks the first list in each
    row to see, first, if the number of argument values corresponds to
    the size of the first dimension, and if this test succeeds, checks
    that all of the argument values match.  This takes time, of course,
    but is much faster than re-evaluation of the basis system.  Even
    this time can be avoided by direct retrieval of the desired array.
    For example, you might set up a vector of argument values called
    "evalargs" along with a matrix of basis function values for these
    argument values called "basismat".  You might want too use tags like
    "args" and "values", respectively for these.  You would then assign
    them to BASISVALUES with code such as basisobj$basisvalues <-
    vector("list",1);  basisobj$basisvalues[[1]] <- list(args=evalargs,
    values=basismat).
  }
}
\details{

  Previous versions of the 'fda' software used the name \code{basis}
  for this class, and the code in Matlab still does.  However, this
  class name was already used elsewhere in the S languages, and there
  was a potential for a clash that might produce mysterious and perhaps
  disastrous consequences.

  To check that an object is of this class, use function
  \code{\link{is.basis}}.

  It is comparatively simple to add new basis types.  The code in
  the following functions needs to be estended to allow for the new
  type: \code{\link{basisfd}},
  \code{\link{getbasismatrix}} and \code{\link{getbasispenalty}}.
  In addition, a new "create" function should be written for the
  new type, as well as functions analogous to \code{\link{fourier}} and
  \code{\link{fourierpen}} for evaluating basis functions for basis
  penalty matrices.

  The "create" function names are rather long, and users who mind
  all that typing might be advised to modify these to versions with
  shorter names, such as "splbas", "conbas", and etc.  However, a
  principle of good programming practice is to keep the code readable,
  preferably by somebody other than the programmer.

  Normally only developers of new basis types will actually need
  to use this function, so no examples are provided.
}
\value{
  an object of class \code{basisfd}, being a list with the following
  components:

  \item{type}{type of basis}
  \item{rangeval}{acceptable range for the argument}
  \item{nbasis}{number of bases}
  \item{params}{
    a vector of parameter values defining the basis.
  }
  \item{dropind}{input argument dropind}
  \item{quadvals}{quadrature values ...}
  \item{values}{
    a list of basis functions and derivatives
  }
  \item{basisvalues}{input argument basisvalues}
}
\references{
  Ramsay, James O., Hooker, Giles and Graves, Spencer (2009)
  \emph{Functional Data Analysis with R and Matlab}, Springer, New
  York.

  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\keyword{smooth}
\seealso{
\code{\link{is.basis}},
\code{\link{is.eqbasis}},
\code{\link{plot.basisfd}},
\code{\link{getbasismatrix}},
\code{\link{getbasispenalty}},
\code{\link{create.bspline.basis}},
\code{\link{create.constant.basis}},
\code{\link{create.exponential.basis}},
\code{\link{create.fourier.basis}},
\code{\link{create.monomial.basis}},
\code{\link{create.polygonal.basis}},
\code{\link{create.power.basis}}
}
% docclass is function
\keyword{smooth}
\keyword{internal}

