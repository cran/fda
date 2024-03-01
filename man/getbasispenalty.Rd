\name{getbasispenalty}
\alias{getbasispenalty}
\title{
  Evaluate a Roughness Penalty Matrix
}
\description{
A basis roughness penalty matrix is the matrix containing
the possible inner products of pairs of basis functions.
These inner products are typically defined in terms of
the value of a derivative or of a linear differential
operator applied to the basis function.  The basis penalty
matrix plays an important role in the computation of
functions whose roughness is controlled by a roughness
penalty.
}
\usage{
getbasispenalty(basisobj, Lfdobj=NULL)
}
\arguments{
\item{basisobj}{
a basis object.
}
\item{Lfdobj}{
an Lfd object
}
}
\value{
a symmetric matrix of order equal to the number of basis functions
defined by the B-spline basis object.  Each element is the inner product
of two B-spline basis functions after taking the derivative.
}
\details{
A roughness penalty for a function $x(t)$ is defined by
integrating the square of either the derivative of  $ x(t) $ or,
more generally, the result of applying a linear differential operator
$L$ to it.  The most common roughness penalty is the integral of
the square of the second derivative, and
this is the default. To apply this roughness penalty, the matrix of
inner products of the basis functions defining this function is
necessary. This function just calls the roughness penalty evaluation
function specific to the basis involved.
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
  \code{\link{eval.penalty}}
}
\examples{

#  set up a B-spline basis of order 4 with 13 basis functions
#  and knots at 0.0, 0.1,..., 0.9, 1.0.
basisobj <- create.bspline.basis(c(0,1),13)
#  compute the 13 by 13 matrix of inner products of second derivatives
penmat <- getbasispenalty(basisobj)
#  set up a Fourier basis with 13 basis functions
#  and and period 1.0.
basisobj <- create.fourier.basis(c(0,1),13)
#  compute the 13 by 13 matrix of inner products of second derivatives
penmat <- getbasispenalty(basisobj)

}
% docclass is function
\keyword{smooth}
