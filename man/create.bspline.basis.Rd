\name{create.bspline.basis}
\alias{create.bspline.basis}
\title{
  Create a B-spline Basis
}
\description{
  Functional data objects are constructed by specifying a set of basis 
  functions and a set of coefficients defining a linear combination of
  these basis functions.  The B-spline basis is a system that is usually
  used for non-periodic functions.  It has the advantages of very fast
  computation and great flexibility.  On the other hand, if the data are
  considered to be periodic, then the Fourier basis is usually
  preferred. 
}
\usage{
create.bspline.basis(rangeval=c(0, 1), nbasis=NULL, norder=4,
      breaks=NULL, dropind=NULL, quadvals=NULL, values=NULL,
      names="bspl") 
}
\arguments{
  \item{rangeval}{
    a numeric vector of length 2 defining the interval over
    which the functional data object can be evaluated;  default value
    \code{c(0,1)}.

    If length(rangeval) == 1 and rangeval <= 0, this is an error.
    Otherwise, if length(rangeval) == 1, rangeval is replaced by c(0,
    rangeval).  

    If length(rangeval)>2 and neither \code{breaks} nor \code{nbasis}
    are provided, this extra long \code{rangeval} argument is assigned
    to \code{breaks}, and then rangeval = range(breaks).  
  }
  \item{nbasis}{
    an integer variable specifying the number of basis functions.
    Default value NULL.

    This 'nbasis' argument is ignored if 'breaks' is supplied, in which
    case nbasis = nbreaks + norder - 2, where nbreaks = length(breaks).
  }
  \item{norder}{
    an integer specifying the order of b-splines, which is one higher
    than their degree. The default of 4 gives cubic splines.
  }
  \item{breaks}{
    a vector specifying the break points defining the b-spline.
    Also called knots, these are a strictly increasing sequence
    of junction points between piecewise polynomial segments.
    They must satisfy \code{breaks[1] = rangeval[1]} and
    \code{breaks[nbreaks] = rangeval[2]}, where \code{nbreaks} is the 
    length of \code{breaks}.  There must be at least 2 values in
    \code{breaks}.  There is a potential for inconsistency among
    arguments \code{nbasis}, \code{norder}, and \code{breaks}.  It is
    resolved as follows:

    If \code{breaks} is supplied, \code{nbreaks = length(breaks)}, and
    \code{nbasis = nbreaks + norder - 2}, no matter what value for
    \code{nbasis} is supplied. 

    If \code{breaks} is not supplied, but \code{nbasis} is,
    \code{nbreaks = nbasis - norder + 2}, and if this turns out to be
    less than 2, an error message results.
    
    If neither \code{breaks} nor \code{nbasis} is supplied,
    \code{nbreaks} is set to 21.

    If \code{breaks} is not provided, it is set to seq(rangeval[1],
    rangeval[2], length=nbreaks).  
  }
  \item{dropind}{
    a vector of integers specifiying the basis functions to
    be dropped, if any.  For example, if it is required that
    a function be zero at the left boundary, this is achieved
    by dropping the first basis function, the only one that
    is nonzero at that point. Default value NULL.
  }
  \item{quadvals}{
    a matrix with two columns and a number of rows equal to the number 
    of argument values used to approximate an integral using Simpson's
    rule.  The first column contains these argument values.
    A minimum of 5 values are required for
    each inter-knot interval, and that is often enough. These
    are equally spaced between two adjacent knots.
    The second column contains the weights used for Simpson's
    rule.  These are proportional to 1, 4, 2, 4, ..., 2, 4, 1.
  }
  \item{values}{
    a list containing the basis functions and their derivatives
    evaluated at the quadrature points contained in the first
    column of \code{ quadvals }.
  }
  \item{names}{
    either a character vector of the same length as the number of basis 
    functions or a single character string to which norder, "." and
    1:nbasis are appended as paste(names, norder, ".", 1:nbreaks,
    sep="").  For example, if norder = 4, this defaults to 'bspl4.1',
    'bspl4.2', ... .   
  }
}
\value{
  a basis object of the type \code{bspline}
}
\details{
  Spline functions are constructed by joining polynomials end-to-end at
  argument values called \emph{break points} or \emph{knots}. First, the
  interval is subdivided into a set of adjoining intervals
  separated the knots.  Then a polynomial of order $m$ and
  degree $m-1$ is defined for each interval.  In order to make the
  resulting piece-wise polynomial smooth, two adjoining polynomials are
  constrained to have their values and all their derivatives up to order
  $m-2$ match at the point where they join.

  Consider as an illustration the very common case where the order is 4 
  for all polynomials, so that degree of each polynomials is 3.  That
  is, the polynomials are \emph{cubic}.  Then at each break point or
  knot, the values of adjacent polynomials must match, and so also for
  their first and second derivatives.  Only their third derivatives will
  differ at the point of junction.

  The number of degrees of freedom of a cubic spline function of this
  nature is calculated as follows.  First, for the first interval, there
  are four degrees of freedom.  Then, for each additional interval, the
  polynomial over that interval now has only one degree of freedom
  because of the requirement for matching values and derivatives.  This
  works out to the following equation:  The number of degrees of freedom
  is equal to the number of interior knots (that is, not counting the
  lower and upper limits) plus the order of the polynomials.

  \emph{B-splines} are a set of special spline functions that can be
  used to construct any such piece-wise polynomial by computing the
  appropriate linear combination.  They derive their computational
  convience from the fact that any B-spline basis function is nonzero
  over at most m adjacent intervals.  The number of basis functions is
  given by the rule above for the number of degrees of freedom.
  
  The number of intervals controls the flexibility of the spline;  the
  more knots, the more flexible the resulting spline will be. But the 
  position of the knots also plays a role.  Where do we position the
  knots?  There is room for judgment here, but two considerations must
  be kept in mind:  (1) you usually want at least one argument value
  between two adjacent knots, and (2)  there should be more knots where
  the curve needs to have sharp curvatures such as a sharp peak or
  valley or an abrupt change of level, but only a few knots are required
  where the curve is changing very slowly.

  This function automatically includes norder replicates of the end
  points rangeval.  By contrast, the analogous splines' package
  functions \link[splines]{splineDesign} and \link[splines]{spline.des}
  do NOT automatically replicate the end points.  To compare answers,
  the end knots must be replicated manually when using
  \link[splines]{splineDesign} or \link[splines]{spline.des}.  
}
\seealso{
  \code{\link{basisfd}}, 
  \code{\link{create.constant.basis}}, 
  \code{\link{create.exponential.basis}}, 
  \code{\link{create.fourier.basis}}, 
  \code{\link{create.monomial.basis}}, 
  \code{\link{create.polygonal.basis}}, 
  \code{\link{create.polynomial.basis}}, 
  \code{\link{create.power.basis}}
  \code{\link[stats]{splineDesign}}
  \code{\link[stats]{spline.des}}
}
\examples{
# The simplest basis currently available with this function:
str(bspl1.1 <- create.bspline.basis(norder=1, breaks=0:1))
plot(bspl1.1)
# 1 basis function, order 1 = degree 0 = step function:  
# constant 1 between 0 and 1.  

str(bspl1.2 <- create.bspline.basis(norder=1, breaks=c(0,.5, 1)))
plot(bspl1.2)
# 2 bases, order 1 = degree 0 = step functions:  
# (1) constant 1 between 0 and 0.5 and 0 otherwise
# (2) constant 1 between 0.5 and 1 and 0 otherwise.

str(bspl2.3 <- create.bspline.basis(norder=2, breaks=c(0,.5, 1)))
plot(bspl2.3)
# 3 bases:  order 2 = degree 1 = linear 
# (1) line from (0,1) down to (0.5, 0), 0 after
# (2) line from (0,0) up to (0.5, 1), then down to (1,0)
# (3) 0 to (0.5, 0) then up to (1,1).

str(bspl3.4 <- create.bspline.basis(norder=3, breaks=c(0,.5, 1)))
plot(bspl3.4)
# 4 bases:  order 3 = degree 2 = parabolas.  
# (1) (x-.5)^2 from 0 to .5, 0 after
# (2) 2*(x-1)^2 from .5 to 1, and a parabola
#     from (0,0 to (.5, .5) to match
# (3 & 4) = complements to (2 & 1).  

# Default b-spline basis
(bSpl4.23 <- create.bspline.basis())
# Cubic bspline (norder=4) with nbasis=23,
# so nbreaks = nbasis-norder+2 = 21, 
# 2 of which are rangeval, leaving 19 Interior knots.

(bSpl4. <- create.bspline.basis(c(-1,1)))
# Same as bSpl4.23 but over (-1,1) rather than (0,1).

# set up the b-spline basis for the lip data, using 23 basis functions,
#   order 4 (cubic), and equally spaced knots.
#  There will be 23 - 4 = 19 interior knots at 0.05, ..., 0.95
lipbasis <- create.bspline.basis(c(0,1), 23)
all.equal(bSpl4.23, lipbasis)
# TRUE 
# plot the basis functions
plot(lipbasis)

(bSpl.growth <- create.bspline.basis(growth$age))

(bSpl.growth6 <- create.bspline.basis(growth$age,norder=6))

}
% docclass is function
\keyword{smooth}
