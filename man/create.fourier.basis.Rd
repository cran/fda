\name{create.fourier.basis}
\alias{create.fourier.basis}
\title{
  Create a Fourier Basis
}
\description{
Create an Fourier basis object defining a set of Fourier
functions with period in argument period.
}
\usage{
create.fourier.basis(rangeval=c(0, 1), nbasis=3,
                     period=width,  dropind=NULL,
                     quadvals=NULL, values=NULL)
}
\arguments{
\item{rangeval}{
a vector of length 2 containing the initial and final
values of the interval over which the functional
data object can be evaluated.  Default value \code{c(0,1)}
}
\item{nbasis}{
an integer variable specifying the number of basis functions.  The
number of basis functions is always odd, even when an even number
is specified, so as to preserve the pairing of sine and cosine
functions.  Default value 3.
}
\item{period}{
the width of any interval over which the Fourier functions repeat
themselves, or are periodic.  The default is the width of the
interval defined in rangeval.
}
\item{dropind}{
a vector of integers specifiying the basis functions to
be dropped, if any.  For example, if it is required that
a function be zero at the left boundary, this is achieved
by dropping the first basis function, the only one that
is nonzero at that point. Default value NULL.
}
\item{quadvals}{
a matrix with two columns and a number of rows equal to the number of
argument values used to approximate an integral using Simpson's rule.
The first column contains these argument values.
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
basisfd, create.bspline.basis, create.constant.basis, create.exponential.basis,
create.monomial.basis, create.polygonal.basis, create.polynomial.basis,
create.power.basis
}
\examples{

# set up the Fourier basis for the monthly temperature data,
#  using 9 basis functions with period 12 months.
monthbasis <- create.fourier.basis(c(0,12), 9, 12.0)
#  plot the basis
plot(monthbasis)

}
% docclass is function
\keyword{smooth}
