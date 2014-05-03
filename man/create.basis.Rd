\name{create.basis}
\alias{create.basis}
\docType{package}
\title{Create Basis Set for Functional Data Analysis}
\description{
  Functional data analysis proceeds by selecting a finite basis set and
  fitting data to it.  The current \code{fda} package supports fitting
  via least squares penalized with lambda times the integral over the
  (finite) support of the basis set of the squared deviations from a
  linear differential operator.
}
\details{
  The most commonly used basis in \code{fda} is probably B-splines.  For
  periodic phenomena, Fourier bases are quite useful.  A constant basis
  is provided to facilitation arithmetic with functional data objects.
  To restrict attention to solutions of certain differential equations,
  it may be useful to use a corresponding basis set such as exponential,
  monomial or power basis sets.

  Power bases support the use of negative and fractional powers, while
  monomial bases are restricted only to nonnegative integer exponents.

  The polygonal basis is essentialy a B-spline of order 2, degree 1.

  The following summarizes arguments used by some or all of the current
  \code{create.basis} functions:

  \itemize{
    \item{rangeval}{
      a vector of length 2 giving the lower and upper limits of the
      range of permissible values for the function argument.

      For \code{bspline} bases, this can be inferred from
      range(breaks).  For \code{polygonal} bases, this can be inferred
      from range(argvals).  In all other cases, this defaults to 0:1.
    }
    \item{nbasis}{
      an integer giving the number of basis functions.

      This is not used for two of the \code{create.basis} functions:
      For \code{constant} this is 1, so there is no need to specify it.
      For \code{polygonal} bases, it is length(argvals), and again there
      is no need to specify it.

      For \code{bspline} bases, if \code{nbasis} is not specified, it
      defaults to (length(breaks) + norder - 2) if \code{breaks} is
      provided.  Otherwise, \code{nbasis} defaults to 20 for
      \code{bspline} bases.

      For \code{exponential} bases, if \code{nbasis} is not specified,
      it defaults to length(ratevec) if \code{ratevec} is provided.
      Otherwise, in \code{fda_2.0.2}, \code{ratevec} defaults to 1,
      which makes \code{nbasis} = 1;  in \code{fda_2.0.4},
      \code{ratevec} will default to 0:1, so \code{nbasis} will then
      default to 2.

      For \code{monomial} and \code{power} bases, if \code{nbasis} is
      not specified, it defaults to length(exponents) if
      \code{exponents} is provided.  Otherwise, \code{nbasis} defaults
      to 2 for \code{monomial} and \code{power} bases.  (Temporary
      exception:  In \code{fda_2.0.2}, the default \code{nbasis} for
      \code{power} bases is 1.  This will be increased to 2 in
      \code{fda_2.0.4}.)

    }

    In addition to \code{rangeval} and \code{nbasis}, all but
    \code{constant} bases have one or two parameters unique to that
    basis type or shared with one other:

    \itemize{
      \item{bspline}{
	Argument \code{norder} = the order of the spline, which is one
	more than the degree of the polynomials used.  This defaults to
	4, which gives cubic splines.

	Argument \code{breaks} = the locations of the break or join
	points;  also called \code{knots}.  This defaults to
	seq(rangeval[1], rangeval[2], nbasis-norder+2).
      }
      \item{polygonal}{
	Argument \code{argvals} = the locations of the break or join
	points;  also called \code{knots}.  This defaults to
	seq(rangeval[1], rangeval[2], nbasis).
      }
      \item{fourier}{
	Argument \code{period} defaults to diff(rangeval).
      }
      \item{exponential}{
	Argument \code{ratevec}.  In \code{fda_2.0.2}, this defaulted to
	1.  In \code{fda_2.0.3}, it will default to 0:1.
      }
      \item{monomial, power}{
	Argument \code{exponents}.  Default = 0:(nbasis-1).  For
	\code{monomial} bases, \code{exponents} must be distinct
	nonnegative integers.  For \code{power} bases, they must be
	distinct real numbers.
      }
    }

    Beginning with \code{fda_2.1.0}, the last 6 arguments for all the
    \code{create.basis} functions will be as follows;  some but not all
    are available in the previous versions of \code{fda}:

    \itemize{
      \item{dropind}{
	a vector of integers specifiying the basis functions to be
        dropped, if any.
      }
      \item{quadvals}{
	a matrix with two columns and a number of rows equal to the
	number of quadrature points for numerical evaluation of the
	penalty integral.  The first column of \code{quadvals} contains
	the quadrature points, and the second column the quadrature
	weights.  A minimum of 5 values are required for each inter-knot
	interval, and that is often enough.  For Simpson's rule, these
	points are equally spaced, and the weights are proportional to
	1, 4, 2, 4, ..., 2, 4, 1.
      }
      \item{values}{
	a list of matrices with one row for each row of \code{quadvals}
	and one column for each basis function.  The elements of the
	list correspond to the basis functions and their derivatives
	evaluated at the quadrature points contained in the first column
	of \code{quadvals}.
      }
      \item{basisvalues}{
	A list of lists, allocated by code such as vector("list",1).
	This field is designed to avoid evaluation of a basis system
	repeatedly at a set of argument values.  Each list within the
	vector corresponds to a specific set of argument values, and
	must have at least two components, which may be tagged as you
	wish.  `The first component in an element of the list vector
	contains the argument values.  The second component in an
	element of the list vector contains a matrix of values of the
	basis functions evaluated at the arguments in the first
	component.  The third and subsequent components, if present,
	contain matrices of values their derivatives up to a maximum
	derivative order.  Whenever function getbasismatrix is called,
	it checks the first list in each row to see, first, if the
	number of argument values corresponds to the size of the first
	dimension, and if this test succeeds, checks that all of the
	argument values match.  This takes time, of course, but is much
	faster than re-evaluation of the basis system.  Even this time
	can be avoided by direct retrieval of the desired array.  For
	example, you might set up a vector of argument values called
	"evalargs" along with a matrix of basis function values for
	these argument values called "basismat".  You might want too use
	tags like "args" and "values", respectively for these.  You
	would then assign them to \code{basisvalues} with code such as
	the following:

	basisobj\$basisvalues <- vector("list",1)

	basisobj\$basisvalues[[1]] <- list(args=evalargs,
	values=basismat)
      }
      \item{names}{
	either a character vector of the same length as the number of
	basis functions or a simple stem used to construct such a
	vector.

	For \code{bspline} bases, this defaults to paste('bspl', norder,
	'.', 1:nbreaks, sep='').

	For other bases, there are crudely similar defaults.
      }
      \item{axes}{
	an optional list used by selected \code{plot} functions to
	create custom \code{axes}.  If this \code{axes} argument is not
	NULL, functions \code{plot.basisfd}, \code{plot.fd},
	\code{plot.fdSmooth} \code{plotfit.fd}, \code{plotfit.fdSmooth},
	and \code{plot.Lfd} will create axes via
	\code{do.call(x$axes[[1]], x$axes[-1])}.  The primary example of
	this is to create \code{CanadianWeather} plots using
	\code{list("axesIntervals")}
      }
    }
  }
}

\author{
  J. O. Ramsay and Spencer Graves
}
\references{
  Ramsay, James O., and Silverman, Bernard W. (2006), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002), \emph{Applied
    Functional Data Analysis}, Springer, New York.
}
\seealso{
  \code{\link{create.bspline.basis}}
  \code{\link{create.constant.basis}}
  \code{\link{create.exponential.basis}}
  \code{\link{create.fourier.basis}}
  \code{\link{create.monomial.basis}}
  \code{\link{create.polygonal.basis}}
  \code{\link{create.power.basis}}
}
%\examples{}
\keyword{smooth}
