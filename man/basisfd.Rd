\name{basisfd}
\alias{basisfd}
\title{
  Define a Functional Basis Object
}
\description{
This is the constructor function for objects of the \code{basisfd} class.
Each function that sets up an object of this class must call this function.
This includes functions \code{create.bspline.basis},
\code{create.constant.basis}, \code{create.fourier.basis},  and
so forth that set up basis objects of a specific type.  Ordinarily, user
of the functional data analysis software will not need to call this function
directly, but these notes are valuable to understanding what the "slots"
or "members" of the \code{basisfd} class are.
}
\usage{
basisfd(type, rangeval, nbasis, params,
        dropind=NULL, quadvals=NULL,
        values=vector("list", 0))
}
\arguments{
\item{type}{
a character string indicating the type of basis.  Currently,
there are eight possible types:

\enumerate{
	\item bspline
	\item const
	\item expon
	\item fourier
	\item monom
	\item polyg
	\item polynom
	\item power
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
a vector of parameter values defining the basis
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
an object of class \code{basisfd}
}
\details{

Older versions of the software used the name \code{basis}
for this class, and the code in Matlab still does.  However, this
class name was already used elsewhere in the S languages, and there
was a potential for a clash that might produce mysterious and perhaps
disastrous consequences.

To check that an object is of this class, use function
\code{\link{is.basis}}.

 It is comparatively simple to add new basis types.  The code in
the following functions needs to be estended to allow for the new
type: \code{\link{basisfd}}, \code{\link{use.proper.basis}},
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
\code{\link{create.polynomial.basis}}, 
\code{\link{create.power.basis}}
}
% docclass is function
\keyword{smooth}
\keyword{internal}

