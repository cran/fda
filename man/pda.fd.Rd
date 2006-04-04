\name{pda.fd}
\alias{pda.fd}
\title{
  Principal Differential Analysis
}
\description{
Principal differential analysis (PDA) estimates a
system of linear differential equations that define functions
that fit the data and their derivatives.
}
\usage{
pda.fd(xfdlist, bwtlist=NULL,
       awtlist=NULL, ufdlist=NULL, nfine=501)
}
\arguments{
\item{xfdlist}{
a list whose members are functional data objects.  Each of these objects
contain one or more functions that are variables to be represented
by a differential equation.  The length of the list is the size of the
system of differential equations. The number of replications must be
the same for each member functional data object.
}
\item{bwtlist}{


a list array with the first two dimensions are equal to the number of
variables in the system (the length of list \code{xfdlist}) and
the last dimension equal to the order of each equation.  The order of
the equations is assumed to be the same for each equation.  If only a
single differential equation is involved, this argument can be an
ordinary list with length equal to the order of the equation.

 Each member of \code{bwtlist} defines a weight function
in the system of linear differential equation.  For any equation,
there can be a term in the equation for each order of derivative
from 0 to m-1, where m is the order of the equation.  These derivatives
can be those of each variable in the system, moreover.  Thus, a two-variable
system of order 3 can have, for each of the two equations, contributions
from derivatives of order 0, 1 and 2 for each of the variables, or
six terms in total.

 The first dimension of \code{bwtlist} corresponds to equations,
and the second dimension to contributions of variables for a fixed derivative
within the equation defined by the first index.  The third index indicates
the order of derivative of the contribution, which is one less than the
index value.

 Each member of \code{bwtlist} is a either a functional
parameter object or a functional data object.  Functional data objects
are converted to functional parameters objects with no smoothing.
Functional parameter objects permit selected weight functions to be
held fixed and not estimated.  These are often constant functions with
the value 0, corresponding to variables and derivatives that make no
contribution.

 For example, the harmonic acceleration differential equation
would be of the form $D^3 x(t) = -b Dx(t)$ so that
\code{bwtlist}
would be an ordinary list of length 3, each functional parameter object
being defined with a constant basis, but only the second object would
have a value of TRUE for the \code{estimate} slot, and the first
and third coefficient functions would have value 0.
}
\item{awtlist}{
a list containing weight functions for forcing functions.

 In addition to terms in each of the equations involving terms corresponding
to each derivative of each variable in the system, each equation can
also have a contribution from one or more exogenous variables, often
called forcing functions.

 This argument defines the weights multiplying these
forcing functions, and is a list array with first dimension equal to
the number of variables or equations, and second dimension equal to the
number of forcing functions. It is assumed that the number of forcing functions is
the same for all equations. If only one forcing function is involved,
or there is only one equation, then \code{awtlist} can be an
ordinary list.

 Each member is a functional parameter object
(or a functional data object) defining a weighting function for a
forcing function.  As with
\code{bwtlist}, each of these weighting functions may be estimated
or held fixed.  If the number of forcing functions actually varies from
equation to equation, one can still use all forcing functions for all
variables, but just weight unwanted ones in a particular equation by a
fixed zero function.

 Each member's functional data object has only a single replicate.
}
\item{ufdlist}{
a list containing forcing functions.
This is a list array of the same size as \code{awtlist} and each
member is a functional data object corresponding to a forcing function.
The number of replicates must be equal to that of the variables themselves,
which is assumed to be the same for all variables.
}
\item{nfine}{
a number of values for a fine mesh.
The estimation of the differential equation involves discrete
numerical quadrature estimates of integrals, and these require
that functions be evaluated at a fine mesh of values of the
argument.  This argument defines the number to use.  If not
constant, this number is set to five times the largest number of
basis functions used to represent any variable in the system.
}
}
\value{
a named list of length 3 with components:

\item{bwtlist}{
a list array of the same dimensions as the
corresponding argument, containing the estimated or fixed weight
functions defining the system of linear differential equations.
}
\item{resfdlist}{
a list of length equal to the number of variables
or equations.  Each members is a functional data object giving the
residual functions or forcing functions defined as the left side
of the equation (the derivative of order m of a variable) minus
the linear fit on the right side.  The number of replicates for each
residual functional data object is the same as that for the variables.
}
\item{awtlist}{
a list of the same dimensions as the corresponding
argument.  Each member is an estimated or fixed weighting function for
a forcing function.
}
}
\seealso{
pca.fd, cca.fd
}
\examples{
#See analyses of daily weather data for examples.
}
% docclass is function
\keyword{smooth}
